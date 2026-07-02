#include "scan_fastq.h"
#include "id_parser.h"
#include <iostream>
#include <array>
#include <algorithm> // std::sort, std::max
#include <string>
#include <utility>
#include <vector>
#include <seqan3/utility/views/chunk.hpp>

fastq_metadata set_up_metadata(sequence_file_input_phred94& fin, const int& id_index, const std::string& suffix) {
    fastq_metadata metadata{};
    auto rec = *fin.begin();
    if (id_index < 0) {
        metadata.id_index = set_id_index(rec.id(), suffix);
    } else {
        metadata.id_index = id_index;
    }
    metadata.n_ID_fields = get_n_ID_fields(rec.id(), metadata.id_index, suffix);
    metadata.illumina_second_id_style = is_valid_second_id_field_illumina(rec.id(), metadata.id_index);
    metadata.print_n_ID_fields();
    if (metadata.illumina_second_id_style) {
        std::cerr << "Assuming the Illumina ID format for obtaining read group IDs\n";
    }
    if (!metadata.valid_n_ID_fields()) {
        std::cerr << "The read group ID will be always 'A'\n";
    }
    return metadata;
}

fastq_metadata scan_fastq_iter(fastq_metadata metadata, sequence_file_input_phred94& fin, const std::string& suffix,
                               const bool& allow_early_termination, const bool& use_index_sequence, 
                               std::size_t& n_check_phred_after_determined) {
    std::size_t n_processed = 0u;
    const bool check = true;
    bool is_phred_determined = metadata.format.determined();
    std::string ID;
    std::size_t len;
    for (auto && rec : fin) {
        ++n_processed;
        ID = get_rg_id(parse_ID(rec.id(), metadata.id_index, suffix, metadata.illumina_second_id_style), metadata.n_ID_fields, check, use_index_sequence);
        auto [count_it, inserted] = metadata.rg_counts.try_emplace(ID, 0u);
        ++count_it->second;
        if (inserted) {
            metadata.rg_ids.push_back(ID); // keep the order of appearance
        }
        len = rec.sequence().size();
        if (std::find(metadata.lengths.begin(), metadata.lengths.end(), len) == metadata.lengths.end()) {
            metadata.lengths.push_back(len); // add the record if it was not already added
        }
        if (!is_phred_determined || n_check_phred_after_determined > 0) {
            metadata.format.update(check_phred(rec.sequence(), rec.base_qualities()));
            if (is_phred_determined) {
                --n_check_phred_after_determined;
            } else {
                is_phred_determined = metadata.format.determined();
            }
        } else if (allow_early_termination == true) {
            break;
        }
    }
    std::cerr << "Scanned " << n_processed << " records\n";
    return metadata;
}

// Merge read groups that share the same run/lane and have an index sequence within
// max_index_mismatch mismatches. Read groups are processed in descending order of read
// counts, so minor (mismatched) index sequences are assigned to the most frequent
// (canonical) index sequence.
void collapse_read_groups(fastq_metadata& metadata, const std::size_t& max_index_mismatch) {
    // (rg_id, count) sorted by count descending, then rg_id ascending (for determinism)
    std::vector<std::pair<std::string, std::size_t>> groups(metadata.rg_counts.begin(), metadata.rg_counts.end());
    std::sort(groups.begin(), groups.end(), [](const auto& lhs, const auto& rhs) {
        if (lhs.second != rhs.second) {
            return lhs.second > rhs.second;
        }
        return lhs.first < rhs.first;
    });

    // a confirmed (canonical) read group: its rg_id, run/lane prefix and index sequence
    struct canonical_group {
        std::string rg_id;
        std::string prefix;
        std::string index;
    };
    std::vector<canonical_group> confirmed;
    std::vector<std::string> collapsed_rg_ids; // canonical rg_ids in count-descending order

    for (const auto& [rg_id, count] : groups) {
        // split rg_id ("run.lane.index") into the "run.lane" prefix and the index sequence
        std::size_t last_dot = rg_id.find_last_of('.');
        std::string prefix = rg_id.substr(0, last_dot);
        std::string index = rg_id.substr(last_dot + 1);

        bool merged = false;
        for (const auto& c : confirmed) {
            if (c.prefix == prefix && index_hamming_distance(index, c.index) <= max_index_mismatch) {
                metadata.rg_id_map[rg_id] = c.rg_id; // link to the most frequent canonical group
                merged = true;
                break;
            }
        }
        if (!merged) {
            metadata.rg_id_map[rg_id] = rg_id; // this group is a new canonical group
            confirmed.push_back(canonical_group{rg_id, prefix, index});
            collapsed_rg_ids.push_back(rg_id);
        }
    }

    std::cerr << "Collapsed " << metadata.rg_ids.size() << " read groups into " << collapsed_rg_ids.size()
              << " (max index mismatch = " << max_index_mismatch << ")\n";
    metadata.rg_ids = collapsed_rg_ids;
}

fastq_metadata scan_fastq(std::filesystem::path & fastq1, std::filesystem::path & fastq2,
                          const int& id_index, const std::string& suffix1, const std::string& suffix2,
                          const phred& prespecified_phred, const bool& allow_early_termination, const bool& use_index_sequence,
                          const std::size_t& max_index_mismatch,
                          std::size_t n_check_phred_after_determined /*= 100000u*/) {
    // iterate over fastq1 to obtain read group IDs (in the {run.lane} format)
    std::cerr << "Scanning the first fastq\n";
    sequence_file_input_phred94 fin{fastq1};
    fastq_metadata metadata = set_up_metadata(fin, id_index, suffix1);
    metadata.format = prespecified_phred;
    if (metadata.format.determined()) {
        n_check_phred_after_determined = 0; //overwrite
    } else if (metadata.n_ID_fields == 0) {
        std::cerr << "Assuming Phred-33 because the read names indicate that the sequencer was DNBSEQ\n";
        metadata.format.update(phred{0B1101});
        n_check_phred_after_determined = 0; //overwrite
    }
    // iterate over the rest of fastq1 to obtain read group IDs (in the {run.lane} format)
    metadata = scan_fastq_iter(metadata, fin, suffix1, allow_early_termination, use_index_sequence, n_check_phred_after_determined);
    try {
        metadata.format.finalize();
    } catch (...) {
        // unfinalizable
        if (!metadata.format.valid()) {
            throw;
        } else {
            std::cerr << "Scanning the second fastq as the format of quality score has not been determined yet.\n";
            sequence_file_input_phred94 fin{fastq2};
            metadata = scan_fastq_iter(metadata, fin, suffix2, allow_early_termination, use_index_sequence, n_check_phred_after_determined);
            metadata.format.finalize();
        }
    }
    // merge read groups whose index sequences differ within the allowed mismatch
    if (use_index_sequence && metadata.n_ID_fields == 7u) {
        collapse_read_groups(metadata, max_index_mismatch);
    }
    return metadata;
}