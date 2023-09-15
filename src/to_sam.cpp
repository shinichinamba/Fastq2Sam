// This file is deprecated
#include "to_sam.h"
#include <cstddef> // std::size_t
#include <stdint.h> // uint64_t
#include <iostream>
#include <string> // for 0-padding
#include <sstream> // for 0-padding
#include <fstream> // write out hash
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/utility/views/chunk.hpp>
#include <seqan3/utility/views/zip.hpp>
#include <seqan3/alphabet/quality/phred94.hpp> // seqan3::phred94
#include "convert_phred_score.h"
#include "sequence_file_input.h"
#include "id_parser.h"
#include "bamhash.h"
using namespace seqan3::literals; 

// sam_flag
auto first_flag = seqan3::sam_flag::paired | seqan3::sam_flag::unmapped | seqan3::sam_flag::mate_unmapped | seqan3::sam_flag::first_in_pair;
auto second_flag = seqan3::sam_flag::paired | seqan3::sam_flag::unmapped | seqan3::sam_flag::mate_unmapped | seqan3::sam_flag::second_in_pair;

// convert fastq files to a sam file
template <class T>
void to_sam(fastq_metadata metadata, T& fin1, T& fin2, sam_file_output_type& fout, const std::vector<seqan3::sam_tag_dictionary>& dicts,
            std::filesystem::path & hash, std::filesystem::path & hash_no_quality, 
            const std::size_t& batch_size, const std::size_t& min_length, const std::string& suffix1, const std::string& suffix2) {
    if (!metadata.enough_info()) {
        throw std::runtime_error("The metadata does not contain enough information.");
    }
    const bool check = false;
    bool error_flag;
    std::string ID1;
    std::string ID2;
    std::string rg_id;
    seqan3::sam_tag_dictionary dict;
    std::vector<sam_record_type> outs(batch_size * 2, sam_record_type{});
    std::size_t n_processed = 0u;
    std::size_t n_skipped = 0u;
    uint64_t sum = 0;
    uint64_t sum_no_quality = 0;
    // iterate over records
    for (auto && [record1, record2] : seqan3::views::zip(fin1 | seqan3::views::chunk(batch_size), fin2 | seqan3::views::chunk(batch_size))) {// && is important!
        std::size_t i = 0u;
        outs.resize(batch_size * 2);
        for (auto && [rec1, rec2] : seqan3::views::zip(record1, record2)) {
            ++n_processed;
            // parse id
            ID1 = parse_ID(rec1.id(), metadata.id_index, suffix1);
            ID2 = parse_ID(rec2.id(), metadata.id_index, suffix2);
            if (ID1 != ID2) {
                throw std::runtime_error("Your pairs don't match. ID in the file 1, " + rec1.id() + "; ID in the file 2, " + rec2.id());
            }
            
            rg_id = get_rg_id(ID1, metadata.n_ID_fields, check);
            error_flag = true;
            for (auto d : dicts) {
                if (rg_id == d.get<"RG"_tag>()) {
                    dict = d;
                    error_flag = false;
                    break;
                }
            }
            if (error_flag == true) {
                throw std::runtime_error("In the " + std::to_string(n_processed) + " th records, the read group ID is not included in the dict: " + rg_id +
                                         " (fastq ID: " + rec1.id() + ")");
            }
            if (rec1.sequence().size() < min_length || rec2.sequence().size() < min_length) {
                ++n_skipped;
            } else {
                // store records
                auto bq1 = convert_phred_score(rec1.base_qualities(), metadata.format);
                auto bq2 = convert_phred_score(rec2.base_qualities(), metadata.format);
                outs.at(i * 2) = sam_record_type{ID1, rec1.sequence(), bq1, first_flag, dict};
                outs.at(i * 2 + 1) = sam_record_type{ID2, rec2.sequence(), bq2, second_flag, dict};

                // hash
                if (!hash.empty()) {
                    hexSum(calc_hash_str(ID1 + "/1" + vec2string(rec1.sequence()) + vec2string(bq1)), sum);
                    hexSum(calc_hash_str(ID2 + "/2" + vec2string(rec2.sequence()) + vec2string(bq2)), sum);
                }
                if (!hash_no_quality.empty()) {
                    hexSum(calc_hash_str((ID1 + "/1" + vec2string(rec1.sequence()))), sum_no_quality);
                    hexSum(calc_hash_str((ID2 + "/2" + vec2string(rec2.sequence()))), sum_no_quality);
                }
                ++i;
            }
        }
        outs.resize(i * 2); // to shrink the vector for the final iteration and skipped records
        fout = outs; //write out
    }
    std::string zero_padding_hex(uint64_t x, int width = 16) {
        std::ostringstream ss;
        ss << std::hex << x;
        std::string ss_str = ss.str();
        ss_str = std::string(std::max(0, width - (int)S.size()), '0') + ss_str;
        return ss_str;
    }
    if (!hash.empty()) {
        std::ofstream writing_file;
        writing_file.open(hash);
        writing_file << zero_padding_hex(sum) << "\t" << std::dec << n_processed << "\n";
        writing_file.close();
    }
    if (!hash_no_quality.empty()) {
        std::ofstream writing_file;
        writing_file.open(hash_no_quality);
        writing_file << zero_padding_hex(sum_no_quality) << "\t" << std::dec << n_processed << "\n";
        writing_file.close();
    }
   std::cerr << "Done. Processed " << n_processed << " record pairs and skipped " << n_skipped << " pairs\n";
};



// explicit instantiation
template void to_sam<sequence_file_input_phred94> (fastq_metadata, sequence_file_input_phred94&, sequence_file_input_phred94&, 
            sam_file_output_type&, const std::vector<seqan3::sam_tag_dictionary>&,
            std::filesystem::path &, std::filesystem::path &, const std::size_t&, const std::size_t&, const std::string&, const std::string&);
template void to_sam<sequence_file_input_phred68solexa> (fastq_metadata, sequence_file_input_phred68solexa&, sequence_file_input_phred68solexa&,
            sam_file_output_type&, const std::vector<seqan3::sam_tag_dictionary>&,
            std::filesystem::path &, std::filesystem::path &, const std::size_t&, const std::size_t&, const std::string&, const std::string&);
