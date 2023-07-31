#ifndef TO_SAM_H_INCLUDED
#define TO_SAM_H_INCLUDED
#include <cstddef> // std::size_t
#include <string>
#include <vector>
#include <filesystem>
#include <seqan3/io/sam_file/output.hpp>
#include <seqan3/io/sam_file/record.hpp>
#include "fastq_metadata.h"
using namespace std::literals; 

// typedef sam_record_type
using types = seqan3::type_list<std::string, std::vector<seqan3::dna5>, std::vector<seqan3::phred94>, seqan3::sam_flag, seqan3::sam_tag_dictionary>;
using fields = seqan3::fields<seqan3::field::id, seqan3::field::seq, seqan3::field::qual, seqan3::field::flag, seqan3::field::tags>;
using sam_record_type = seqan3::sam_record<types, fields>;

using sam_file_fields = seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::ref_id, seqan3::field::ref_offset, seqan3::field::cigar, 
    seqan3::field::mapq, seqan3::field::qual, seqan3::field::flag, seqan3::field::mate, seqan3::field::tags, seqan3::field::header_ptr>;
using sam_file_types = seqan3::type_list<seqan3::format_bam, seqan3::format_sam>;
using sam_file_output_type = seqan3::sam_file_output<sam_file_fields, sam_file_types>; // same as default, but explicitly defined

template <class T>
void to_sam(fastq_metadata metadata, T& fin1, T& fin2, sam_file_output_type& fout, const std::vector<seqan3::sam_tag_dictionary>& dicts,
            std::filesystem::path & hash, std::filesystem::path & hash_no_quality, 
            const std::size_t& batch_size, const std::size_t& min_length, const std::string& suffix1, const std::string& suffix2);

#endif