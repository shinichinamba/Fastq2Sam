#include "fastq_metadata.h"
#include <seqan3/core/debug_stream.hpp>

//constructor
fastq_metadata::fastq_metadata() {
    format = phred();
    id_index = 9999u; // implausibly large value
    n_ID_fields = 9999u; // implausibly large value
    rg_ids = std::vector<std::string>{};
    lengths = std::vector<std::size_t>{};
}

// n_ID_fields
bool fastq_metadata::valid_n_ID_fields() {
    return (n_ID_fields == 5u) || (n_ID_fields == 7u);
}

void fastq_metadata::print_n_ID_fields() {
    if (n_ID_fields == 7u) {
        seqan3::debug_stream << "The estimated ID format is CASAVA-1.8\n";
    } else if (n_ID_fields == 5u) {
        seqan3::debug_stream << "The estimated ID format is the previous format used before CASAVA-1.8\n";
    } else {
        seqan3::debug_stream << "The number of the fields in the fastq ID is " << n_ID_fields << ", not 5 or 7.\n";
    }
}

// whether we can use the metadata for to_sam()
bool fastq_metadata::enough_info() {
    return format.determined() && 
        id_index != 9999u &&
        n_ID_fields != 9999u &&
        rg_ids.size() > 0;
}

// for n_ID_fields, use print_n_ID_fields()
void fastq_metadata::print() {
    seqan3::debug_stream
        << "The format of quality score: " << format.to_string() << '\n'
        << "Read group IDs (N=" << rg_ids.size() << ") : " << rg_ids << '\n'
        << "Read lengths (N=" << lengths.size() << ") : " << lengths << '\n';
    /*
    std::size_t id_index;
    */
}
