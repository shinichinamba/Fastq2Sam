#ifndef PHRED_H_INCLUDED
#define PHRED_H_INCLUDED

#include <string>
#include <vector>
#include <bitset>
#include <seqan3/alphabet/quality/phred94.hpp> // seqan3::phred94
#include <seqan3/alphabet/nucleotide/dna5.hpp> // seqan3::dna5

struct phred {
    std::bitset<4> format; // <phred33 (max 41)><phred33 (max >41)><phred64><solexa64> : 3210
    std::bitset<4> suggestive; // <phred33 (max 41)><phred33 (max >41)><phred64><solexa64> : 3210
    //constructor
    phred();
    phred(const int& x, const int y = 15);
    phred(const std::bitset<4>& x, const std::bitset<4> y = {0B1111});
    //member functions
    std::vector<std::string> to_string_vec() const;
    std::string to_string() const;
    bool determined() const;
    bool none() const;
    bool valid() const;
    phred operator&(const phred&) const;
    void operator&=(const phred&);
    bool operator==(const phred&);
    void update(const phred&); // a fail-safe version of &=. Will throw runtime error if the result is invalid
    void finalize(); // assume one phred format from an ambiguous format
};

// check the phred format of quality score
phred check_phred(const std::vector<seqan3::dna5>&, const std::vector<seqan3::phred94>&);

// construct a phred object from a string
phred string_to_phred(const std::string&);
#endif