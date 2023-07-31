#ifndef BAMHASH_H_INCLUDED
#define BAMHASH_H_INCLUDED

#include <openssl/md5.h> // bamhash
#include <vector>
#include <string>
// codes derived from BamHash (https://github.com/DecodeGenetics/BamHash)
union hash_t {
  unsigned char c[16];
  struct {
    uint64_t low;
    uint64_t high;
  } p;
};

hash_t str2md5(const char *str, int length);
void hexSum(hash_t out, uint64_t& sum);
// end of functions derived from BamHash 

hash_t calc_hash_str(const std::string & x);
template<class T>
std::string vec2string(const std::vector<T> & x);

#endif
