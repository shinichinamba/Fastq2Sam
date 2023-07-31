#include "scan_fastq.h"
#include "id_parser.h"
#include <iostream>
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
    metadata.print_n_ID_fields();
    if (!metadata.valid_n_ID_fields()) {
        std::cerr << "The read group ID will be always 'A'\n";
    }
    return metadata;
}

fastq_metadata scan_fastq_iter(fastq_metadata metadata, sequence_file_input_phred94& fin, const std::size_t& batch_size, const std::string& suffix,
                               std::size_t& n_check_phred_after_determined) {
    std::size_t n_processed = 0u;
    const bool check = true;
    bool is_phred_determined = false;
    std::string ID;
    std::size_t len;
    for (auto && records : fin | seqan3::views::chunk(batch_size)) { // `&&` is important because seqan3::views::chunk returns temporaries!
        for (auto rec : records) { // scan all records
            ++n_processed;
            ID = get_rg_id(parse_ID(rec.id(), metadata.id_index, suffix), metadata.n_ID_fields, check);
            if (std::find(metadata.rg_ids.begin(), metadata.rg_ids.end(), ID) == metadata.rg_ids.end()) {
                metadata.rg_ids.push_back(ID); // add the record if it was not already added
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
            }

        }
    }

    std::cerr << "Scanned " << n_processed << " records\n";
    return metadata;
}

fastq_metadata scan_fastq(std::filesystem::path & fastq1, std::filesystem::path & fastq2, 
                          const std::size_t& batch_size, const int& id_index, const std::string& suffix1, const std::string& suffix2,
                          std::size_t n_check_phred_after_determined /*= 10000u*/) {
    // iterate over fastq1 to obtain read group IDs (in the {run.lane} format)
    std::cerr << "Scanning the first fastq\n";
    sequence_file_input_phred94 fin{fastq1};
    fastq_metadata metadata = set_up_metadata(fin, id_index, suffix1);
    // iterate over the rest of fastq1 to obtain read group IDs (in the {run.lane} format)
    metadata = scan_fastq_iter(metadata, fin, batch_size, suffix1, n_check_phred_after_determined);
    try {
        metadata.format.finalize();
    } catch (...) {
        // unfinalizable
        if (!metadata.format.valid()) {
            throw;
        } else {
            std::cerr << "Scanning the second fastq as the format of quality score has not been determined yet.\n";
            sequence_file_input_phred94 fin{fastq2};
            metadata = scan_fastq_iter(metadata, fin, batch_size, suffix2, n_check_phred_after_determined);
            metadata.format.finalize();
        }
    }
    return metadata;
}