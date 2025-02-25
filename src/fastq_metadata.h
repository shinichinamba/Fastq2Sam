#ifndef FASTQ_METADATA_H_INCLUDED
#define FASTQ_METADATA_H_INCLUDED
#include <cstddef> // std::size_t
#include "phred.h"
struct fastq_metadata {
    phred format; 
    std::size_t id_index;
    std::size_t n_ID_fields;
    std::vector<std::string> rg_ids;
    std::vector<std::size_t> lengths;
    //constructor
    fastq_metadata();
    //member functions
    bool valid_n_ID_fields();
    void print_n_ID_fields();
    bool enough_info();
    void print();
};
#endif