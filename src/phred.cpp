#include "phred.h"
#include <iostream>
#include <cstddef> // std::size_t
using namespace std::literals; 
using namespace seqan3::literals; 

// struct phred
// constructor
phred::phred() {
    format = std::bitset<4>{0B1111};
}

phred::phred(const int& x) {
    format = x;
}
phred::phred(const std::bitset<4>& x) {
    format = x;
}

std::vector<std::string> phred::to_string_vec() const {
    std::vector<std::string> res;
    if (format[2] && format[3]) {
        res.push_back("Phred33 (undetermined max quality)"s);
    } else if (format[2]) {
        res.push_back("Phred33 (max quality > 41)"s);
    } else if (format[3]) {
        res.push_back("Phred33 (max quality <= 41)"s);
    }
    if (format[1]) {
        res.push_back("Phred64"s);
    }
    if (format[0]) {
        res.push_back("Solexa64"s);
    }
    return res;
}
std::string phred::to_string() const {
    if (format.all()) {
        return "undetermined"s;
    } else if (this->none()) {
        return "none of possible formats"s;
    }
    std::vector<std::string> strvec = this->to_string_vec();
    std::string res{""s};
    for (auto str : strvec) {
        if (res == ""s) {
            res = str;
        } else {
            res = res + " / "s + str;
        }
    }
    return res;
}

bool phred::determined() const {
    return (format.count() == 1) || (format == std::bitset<4>{0B1100});
}

bool phred::none() const {
    return format.none();
}

bool phred::valid() const {
    if (this->none()) {
        return false;
    } else if (this->determined()) {
        return true;
    } else {
        if (format[1] && !format[0]) {
            return false;
        } else if (format[3] && !format[2]) {
            return false;
        } else {
            return true;
        }
    }
}

phred phred::operator&(const phred& x) const {
    phred res(format & x.format);
    return res;
}

void phred::operator&=(const phred& x) {
    *this = *this & x;
}

bool phred::operator==(const phred& x) {
    return format == x.format; 
}


void phred::update(const phred& x) {
    phred phred_ = *this & x;
    if (phred_.none()) {
        throw std::runtime_error("No phred format got possible. Original: ["s + this->to_string() + "] and ["s + x.to_string() + "]"s);
    } else if (!phred_.valid()) {
        throw std::runtime_error("The phred format got invalid. Original: ["s + this->to_string() + "] and ["s + x.to_string() + 
                                 "]; New: ["s + phred_.to_string() + "]"s);
    }
    *this = phred_;
}

void phred::finalize() {
    if (this->determined()) {
        if (format == std::bitset<4>{0B1100}) {
            format = std::bitset<4>{0B1000}; 
        }
        return;
    }
    if (!this->valid()) {
        throw std::runtime_error("Cannot finalize the invalid phred format: ["s + this->to_string() + "]"s);
    }

    if (format == std::bitset<4>{0B0011} | format == std::bitset<4>{0B0111}) {
        phred phred_(0B0010);
        std::cerr << "The phred format is [" << this->to_string() << "] and is set as [" << phred_.to_string() << "]\n";
        *this = phred_;
    } else {
        throw std::runtime_error("Cannot finalize the phred format: ["s + this->to_string() + "]"s);
    }
}
// end of struct phred

// see https://bi.biopapyrus.jp/rnaseq/qc/fastq-quality-score.html
int phred33_specific_chars_max = (':'_phred94).to_phred();
int phred64_suggestive_chars_min = ('K'_phred94).to_phred();
std::vector<seqan3::phred94> phred64_0_chars = {'@'_phred94, 'B'_phred94};
seqan3::phred94 solexa_0_char = ';'_phred94;
std::vector<seqan3::phred94> solexa_suggestive_chars = {'<'_phred94, '='_phred94, '>'_phred94, '?'_phred94}; // except for solexa_0_char

// check the phred format of quality score
phred check_phred(const std::vector<seqan3::dna5>& seq_vec, const std::vector<seqan3::phred94>& qual_vec) {
    phred res(0B1111);
    int qp;
    for (std::size_t i = 0; i < seq_vec.size(); i++) {
        if (seq_vec[i] == 'N'_dna5) {
            if (std::find(phred64_0_chars.begin(), phred64_0_chars.end(), qual_vec[i]) != phred64_0_chars.end()) {
                return phred{0B0010}; // determined
            } else if (qual_vec[i] == solexa_0_char) {
                return phred{0B0001}; // determined
            }
        }
        qp = qual_vec[i].to_phred();
        if (qp <= phred33_specific_chars_max) {
            return phred{0B1100}; // determined
        } else if (qp >= phred64_suggestive_chars_min) {
            res &= phred{0B0111}; // Phred+64 or Solexa+64. we cannot differentiate Phred+64 from Solexa+64
        } else if (std::find(solexa_suggestive_chars.begin(), solexa_suggestive_chars.end(), qual_vec[i]) != solexa_suggestive_chars.end()) {
            res &= phred{0B1101}; // Phred+33 or Solexa+64
        }
    }
    return res;
}

phred string_to_phred(const std::string& phred_string) {
    if (phred_string == "auto"s) {
        return phred{0B1111};
    } else if (phred_string == "phred33"s) {
        return phred{0B0100};
    } else if (phred_string == "phred64"s) {
        return phred{0B0010};
    } else if (phred_string == "solexa64"s) {
        return phred{0B0001};
    } else {
        throw std::runtime_error("Cannot determine the phred format: "s + phred_string);
    }
}