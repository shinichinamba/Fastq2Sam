#ifndef SEQUENCE_FILE_INPUT_H_INCLUDED
#define SEQUENCE_FILE_INPUT_H_INCLUDED
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/alphabet/quality/phred94.hpp> // seqan3::phred94
#include <seqan3/alphabet/quality/phred68solexa.hpp> // seqan3::phred68solexa

struct sequence_file_input_traits_phred94 : seqan3::sequence_file_input_default_traits_dna
{
    using quality_alphabet = seqan3::phred94; 
};
using sequence_file_input_phred94 = seqan3::sequence_file_input<sequence_file_input_traits_phred94>;
// The phred94 covers entire codes for both phred33 and phred68solexa

// typedef
struct sequence_file_input_traits_phred68solexa : seqan3::sequence_file_input_default_traits_dna
{
    using quality_alphabet = seqan3::phred68solexa; // for old Illumina format (<1.8)
};
using sequence_file_input_phred68solexa = seqan3::sequence_file_input<sequence_file_input_traits_phred68solexa>;

#endif