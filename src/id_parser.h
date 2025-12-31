#ifndef ID_PARSER_H_INCLUDED
#define ID_PARSER_H_INCLUDED
#include <vector>
#include <string>
#include <cstddef> // std::size_t
std::string get_rg_id(const std::array<std::string, 2>& id, const std::size_t& n_fields, const bool& check = false);
std::string remove_suffix(std::string& query, const std::string& suffix);
std::vector<std::string> split_by_spaces(std::string x);
std::array<std::string, 2> parse_ID(const std::string& ID, const std::size_t& id_index, const std::string& suffix, const bool& also_return_second_id = false);
std::size_t get_n_ID_fields(const std::string& ID, const std::size_t& index, const std::string& suffix);
std::size_t set_id_index(const std::string& ID, const std::string& suffix);
bool is_valid_second_id_field_illumina(const std::string& ID, const std::size_t& idx);
#endif