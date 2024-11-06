#ifndef SCAN_FASTQ_H_INCLUDED
#define SCAN_FASTQ_H_INCLUDED
#include "sequence_file_input.h"
#include "fastq_metadata.h"
#include <cstddef> // std::size_t

/* 
fastq_metadata set_up_metadata(sequence_file_input_phred94& fin, const int& id_index, const std::string& suffix);
fastq_metadata scan_fastq_iter(fastq_metadata metadata, sequence_file_input_phred94& fin, const std::size_t& batch_size, const std::string& suffix,
        std::size_t& n_check_phred_after_determined);
*/

fastq_metadata scan_fastq(std::filesystem::path & fastq1, std::filesystem::path & fastq2, 
                          const std::size_t& batch_size, const int& id_index, const std::string& suffix1, const std::string& suffix2,
                          const phred& prespecified_phred, std::size_t n_check_phred_after_determined = 10000u);
#endif