#include "id_parser.h"
#include "split.h"
#include <iostream>
#include <stdexcept> // std::runtime_error
#include <algorithm> //std::replace
#include <regex>
#include <array>
using namespace std::literals; 

std::regex DNBSEQ_regex(R"(^(\w+)L(\d+)C\d{3}R\d+$)");

// get read group ID from fastq ID
std::string get_rg_id(const std::array<std::string, 2>& id, const std::size_t& n_fields, const bool& check /*=false"*/) {
    std::vector<std::string> split_id = split(id[0], ':');
    std::string rg_id;
    if (check == true && n_fields != 0 && n_fields != split_id.size())
    {
        throw std::runtime_error("Inconsistent ID formats; expected " + std::to_string(n_fields) + " fields, but found " + std::to_string(split_id.size()));
    }
    if (n_fields == 7)
    {
        if (id[1] != "")
        {
            std::vector<std::string> second_split_id = split(id[1], ':');
            rg_id = split_id[2] + "." + split_id[3] + "." + second_split_id[3];
        }
        else
        {
            rg_id = split_id[2] + "." + split_id[3];
        }
    }
    else if (n_fields == 5)
    {
        rg_id = split_id[0] + "." + split_id[1]; 
        // Although the first field is the machine ID in the previous Illumina fastq ID format,
        // we use this field as the run ID is not included in the fastq ID.
    }
    else if (n_fields == 0) // DNBSEQ
    {
        // DNBSEQ
        std::smatch matched;
        bool is_matched = std::regex_match(id[0], matched, DNBSEQ_regex);
        if (check == true && is_matched == false) {
            throw std::runtime_error("Inconsistent ID formats; expected the DNBSEQ format, but got " + id[0]);
        }
        rg_id = matched[1].str() + "." + matched[2].str(); 
    }
    else
    {
        rg_id = "A"s;
    }
    return rg_id;
}

// remove a suffix
std::string remove_suffix(std::string& query, const std::string& suffix) {
    if (query.length() <= suffix.length())
    {
        return query; // query.substr raises error due to integer overflow if query.length() - suffix.length() < 0
    }
    if (query.substr(query.length() - suffix.length(), query.length()) == suffix)
    {
        return query.substr(0, query.length() - suffix.length());
    }
    else
    {
        return query;
    }
}

// split by either tab or spaces
std::vector<std::string> split_by_spaces(std::string x) {
    std::replace(x.begin(), x.end(), '\t', ' ');
    return split(x, ' ', true); 
}

// parse ID
std::array<std::string, 2> parse_ID(const std::string& ID, const std::size_t& id_index, const std::string& suffix, const bool& also_return_second_id /*= false*/) {
    auto split_IDs = split_by_spaces(ID);
    std::string first_id = remove_suffix(split_IDs[id_index], suffix);
    std::string second_id;
    if (also_return_second_id)
    {
        second_id = split_IDs[id_index + 1];
    }
    else
    {
        second_id = "";
    }
    std::array<std::string, 2> parsed{first_id, second_id};
    return parsed;
}

// obtain the number of fields in the fastq ID
std::size_t get_n_ID_fields(const std::string& ID, const std::size_t& index, const std::string& suffix) {
    const bool also_return_second_id = false;
    std::string parsed_ID = parse_ID(ID, index, suffix, also_return_second_id)[0];
    std::size_t n_ID_fields = (split(parsed_ID, ':')).size();
    if (n_ID_fields == 1 && std::regex_search(parsed_ID, DNBSEQ_regex)) {
        n_ID_fields = 0; // DNBSEQ -> set n_ID_fields as 0
    }
    return n_ID_fields;
}

// determine the id_index
std::size_t set_id_index(const std::string& ID, const std::string& suffix) {
    auto split_IDs = split_by_spaces(ID);
    std::size_t n_ID_fields;
    std::size_t id_index = 0;
    bool is_id_found = false;
    for (auto split_ID : split_IDs) {
        n_ID_fields = get_n_ID_fields(split_ID, 0, suffix);
        if (n_ID_fields == 7u || n_ID_fields == 5u || n_ID_fields == 0u) {
            is_id_found = true;
            break;
        }
        ++id_index;
    }
    if (is_id_found) {
        std::cerr << "The " << id_index + 1 << " th part of the ID row will be used as a read ID.\n";
    }
    else
    {
        std::cerr << "No part of the ID row was in the illumina/DNBSEQ formats. The first part of the ID row will be used as a read ID.\n";
        id_index = 0;
    }
    return id_index;
}

// Illumina BaseSpace-style
bool is_valid_second_id_field_illumina(const std::string& ID, const std::size_t& idx) {
    auto split_IDs = split_by_spaces(ID);
    std::string first_id = split_IDs[idx];
    bool is_valid = false;
    if ((split(first_id, ':')).size() == 7u && split_IDs.size() > idx) {
        std::string second_id = split_IDs[idx + 1];
        std::vector<std::string> split_second_id = split(second_id, ':');
        is_valid = split_second_id.size() == 4u;
    }
    return is_valid;
}
