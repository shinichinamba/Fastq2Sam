#ifndef SPLIT_H_INCLUDED
#define SPLIT_H_INCLUDED
#include <vector>
#include <string>

// split function for std::string. see https://marycore.jp/prog/cpp/std-string-split/
template<class T> 
std::vector<std::string> split(const std::string& s, const T& separator, bool ignore_empty = 0, bool split_empty = 0);

#endif