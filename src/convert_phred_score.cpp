#include "convert_phred_score.h"
#include <cstddef> // std::size_t
#include <string>
using namespace std::literals; 

// check phred and return the input as is
std::vector<seqan3::phred94> convert_phred_score(const std::vector<seqan3::phred94>& qual_vec, const phred& format) {
    if (format == phred{0B1000} || format == phred{0B0100}) {
        return qual_vec;
    } else {
        throw "quality scores were read as seqan3::phred94, but the phred format was incompatible: "s + format.to_string();
    }
}

int solexa_to_phred(const seqan3::phred68solexa& qual) {
    // solexa format
    // see https://cell-innovation.nig.ac.jp/SurfWiki/FASTQ.html
    // see also https://cell-innovation.nig.ac.jp/SurfWiki/solexaQuality.html
    int res;
    switch (qual.to_phred()) {
    case -5:
        res = 1;
        break;
    case -4:
        res = 1;
        break;
    case -3:
        res = 2;
        break;
    case -2:
        res = 2;
        break;
    case -1:
        res = 3;
        break;
    case 0:
        res = 3;
        break;
    case 1:
        res = 4;
        break;
    case 2:
        res = 4;
        break;
    case 3:
        res = 5;
        break;
    case 4:
        res = 5;
        break;
    case 5:
        res = 6;
        break;
    case 6:
        res = 7;
        break;
    case 7:
        res = 8;
        break;
    case 8:
        res = 9;
        break;
    case 9:
        res = 10;
        break;
    default:
        res = qual.to_phred();
        break;
    }
    return res;
}

// convert phred68solexa to phred94
std::vector<seqan3::phred94> convert_phred_score(const std::vector<seqan3::phred68solexa>& qual_vec, const phred& format) {
    std::vector<seqan3::phred94> res(qual_vec.size());
    seqan3::phred94 q;
    if (format == phred{0B0010}) {
        for (std::size_t i=0 ; i<qual_vec.size() ; i++) {
            q.assign_phred(qual_vec[i].to_phred());
            res[i] = q;
        }
    } else if (format == phred{0B0001}) {
        for (std::size_t i=0 ; i<qual_vec.size() ; i++) {
            q.assign_phred(solexa_to_phred(qual_vec[i]));
            res[i] = q;
        }
    } else {
        throw "quality scores were read as seqan3::phred68solexa, but the phred format was incompatible: "s + format.to_string();
    }
    return res;
}
