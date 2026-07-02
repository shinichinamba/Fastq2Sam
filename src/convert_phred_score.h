#ifndef CONVERT_PHRED_SCORE_H_INCLUDED
#define CONVERT_PHRED_SCORE_H_INCLUDED
#include <vector>
#include <seqan3/alphabet/quality/phred94.hpp> // seqan3::phred94
#include <seqan3/alphabet/quality/phred68solexa.hpp> // seqan3::phred68solexa
#include "phred.h"
std::vector<seqan3::phred94> convert_phred_score(const std::vector<seqan3::phred94>& qual_vec, const phred& format);
std::vector<seqan3::phred94> convert_phred_score(const std::vector<seqan3::phred68solexa>& qual_vec, const phred& format);

#endif