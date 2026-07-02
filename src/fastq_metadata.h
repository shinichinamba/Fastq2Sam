#ifndef FASTQ_METADATA_H_INCLUDED
#define FASTQ_METADATA_H_INCLUDED
#include <cstddef> // std::size_t
#include <string>
#include <vector>
#include <unordered_map>
#include "phred.h"
struct fastq_metadata {
    phred format;
    std::size_t id_index;
    std::size_t n_ID_fields;
    bool illumina_second_id_style;
    std::vector<std::string> rg_ids;
    std::vector<std::size_t> lengths;
    std::unordered_map<std::string, std::size_t> rg_counts; // rg_id -> the number of reads
    std::unordered_map<std::string, std::string> rg_id_map; // raw rg_id -> canonical rg_id (empty means no collapsing)
    //constructor
    fastq_metadata();
    //member functions
    bool valid_n_ID_fields();
    void print_n_ID_fields();
    bool enough_info();
    void print();
};
#endif