#include "bamhash.h"
#include <seqan3/alphabet/quality/phred94.hpp> // seqan3::phred94
#include <seqan3/alphabet/nucleotide/dna5.hpp> // seqan3::dna5

// codes derived from BamHash (https://github.com/DecodeGenetics/BamHash)
hash_t str2md5(const char *str, int length) {
  hash_t out;
  MD5((unsigned char *)str, length, (unsigned char *)(out.c));
  return out;
}

void hexSum(hash_t out, uint64_t& sum) {
  sum += out.p.low;
}
// end of functions derived from BamHash 

// add hash values
hash_t calc_hash_str(const std::string & x) {
    hash_t hex = str2md5(x.c_str(), x.size());
    return hex;
}

template<class T>
std::string vec2string(const std::vector<T> & x) {
    std::string res{};
    for (T i : x)
    {
        res += i.to_char();
    }
    return res;
};

// explicit instantiation
template std::string vec2string<seqan3::phred94> (const std::vector<seqan3::phred94> &);
template std::string vec2string<seqan3::dna5> (const std::vector<seqan3::dna5> &);
