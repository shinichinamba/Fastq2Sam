#include <cstdlib> //EXIT_SUCCESS, EXIT_FAILURE
#include <filesystem>
#include <string>
// #include <algorithm>  // std::unique
#include <vector>
#include <utility>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sam_file/output.hpp>
#include <seqan3/io/sam_file/record.hpp>
#include <seqan3/utility/views/chunk.hpp>
#include <seqan3/utility/views/zip.hpp>
#include <seqan3/contrib/stream/bgzf.hpp> // seqan3::contrib::bgzf_thread_count
#include <seqan3/alphabet/quality/phred68solexa.hpp> // seqan3::phred68solexa
#include <seqan3/alphabet/quality/phred94.hpp> // seqan3::phred94
#include <sharg/all.hpp> // argparser for seqan3.3

#include <stdint.h> // uint64_t
#include <openssl/md5.h> // bamhash
#include <fstream> // write out hash

using namespace std::literals; 
using namespace seqan3::literals; 

// program name and version
const std::string PG{"fastq2sam"s};
const std::string VER{"0.0.1"s};


// default values of the arguments
struct cmd_arguments
{
    std::filesystem::path fastq1{};
    std::filesystem::path fastq2{};
    std::filesystem::path out{};
    std::string sample_name{};
    std::filesystem::path hash{};
    std::filesystem::path hash_no_quality{};
    std::string platform{"ILLUMINA"s};
    unsigned int min_length{0u};
    int id_index{-1};
    std::string suffix1{"/1"s};
    std::string suffix2{"/2"s};
    unsigned int batch_size{1000000u};
    unsigned int bam_writers{4u};
    bool no_pg{false};
};

void initialise_parser(sharg::parser & parser, cmd_arguments & args)
{
    parser.info.author = "Shinichi Namba";
    parser.info.short_description = "Converting paired-end illumina fastq files into an unmapped sam/bam file while properly handling the RG tags.";
    parser.info.version = VER;
    parser.info.short_copyright = "GPL v3.0";

    parser.add_option(args.fastq1,
                      sharg::config{.short_id = '1',
                                    .long_id = "fastq1",
                                    .description = "The first file of the paired-end fastq files.", 
                                    .required = true,
                                    .validator = sharg::input_file_validator{{"fq", "fastq", "gz"}}});
    parser.add_option(args.fastq2,
                      sharg::config{.short_id = '2',
                                    .long_id = "fastq2",
                                    .description = "The second file of the paired-end fastq files.", 
                                    .required = true,
                                    .validator = sharg::input_file_validator{{"fq", "fastq", "gz"}}});
    parser.add_option(args.out,
                      sharg::config{.short_id = 'o',
                                    .long_id = "out",
                                    .description = "The output sam/bam file. If not set, this program will only scan the first fastq file.", 
                                    .required = false,
                                    .validator = sharg::output_file_validator{sharg::output_file_open_options::open_or_create, {"sam", "bam"}}});
    parser.add_option(args.sample_name,
                      sharg::config{.short_id = 'n',
                                    .long_id = "sample-name",
                                    .description = "Sample ID.", 
                                    .required = true});
    parser.add_option(args.hash,
                      sharg::config{.long_id = "hash",
                                    .description = "The output file to write the read hash, which is compatible to BamHash."});
    parser.add_option(args.hash_no_quality,
                      sharg::config{.long_id = "hash-no-quality",
                                    .description = "The output file to write the read hash, calculated without read quality."});
    parser.add_option(args.platform,
                      sharg::config{.short_id = 'p',
                                    .long_id = "platform",
                                    .description = "The platform/technology used to produce the reads. This value will be used for the PL tag.", 
                                    .advanced = true});
    parser.add_option(args.min_length,
                      sharg::config{.short_id = 'm',
                                    .long_id = "min-length",
                                    .description = "The minimum read length to be written to the output file.", 
                                    .advanced = true});
    parser.add_option(args.id_index,
                      sharg::config{.short_id = 'i',
                                    .long_id = "id-index",
                                    .description = "The column index of ID to be parsed (0-based). With the default, the index will be automatically determined.", 
                                    .advanced = true,
                                    .validator = sharg::arithmetic_range_validator{-1, 10}});
    parser.add_option(args.suffix1,
                      sharg::config{.long_id = "suffix1",
                                    .description = "The suffix of the read ID in the first fastq to be removed.", 
                                    .advanced = true});
    parser.add_option(args.suffix2,
                      sharg::config{.long_id = "suffix2",
                                    .description = "The suffix of the read ID in the second fastq to be removed.", 
                                    .advanced = true});
    parser.add_option(args.batch_size,
                      sharg::config{.short_id = 'b',
                                    .long_id = "batch-size",
                                    .description = "The batch size to read fastq files at a time.", 
                                    .advanced = true,
                                    .validator = sharg::arithmetic_range_validator{1, 100000000}});
    parser.add_option(args.bam_writers,
                      sharg::config{.long_id = "bam-writers",
                                    .description = "The number of threads used to write a bam file.", 
                                    .advanced = true,
                                    .validator = sharg::arithmetic_range_validator{1, 100}});
    parser.add_flag(args.no_pg,
                      sharg::config{.long_id = "no-pg",
                                    .description = "Do not write the PG tag (for testing purpose).", 
                                    .advanced = true});                                    
}

// typedef
struct my_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4; // instead of dna5
 
    template <typename alph>
    using sequence_container = std::vector<alph>; // must be defined as a template!
};
using sequence_file_input_my_traits = seqan3::sequence_file_input<my_traits>;

struct sequence_file_input_traits_phred64 : seqan3::sequence_file_input_default_traits_dna
{
    using quality_alphabet = seqan3::phred68solexa; // for old Illumina format (<1.8)
};
using sequence_file_input_phred64 = seqan3::sequence_file_input<sequence_file_input_traits_phred64>;

struct sequence_file_input_traits_phred94 : seqan3::sequence_file_input_default_traits_dna
{
    using quality_alphabet = seqan3::phred94; 
};
using sequence_file_input_phred94 = seqan3::sequence_file_input<sequence_file_input_traits_phred94>;
// The phred94 covers entire codes for both phred33 and phred68solexa



// derived from BamHash (https://github.com/DecodeGenetics/BamHash)
union hash_t {
  unsigned char c[16];
  struct {
    uint64_t low;
    uint64_t high;
  } p;
};

hash_t str2md5(const char *str, int length) {
  hash_t out;
  MD5((unsigned char *)str, length, (unsigned char *)(out.c));
  return out;
}

void hexSum(hash_t out, uint64_t& sum) {
  sum += out.p.low;
}

// add hash values
auto calc_hash_str = [](const std::string & x)
{
    auto hex = str2md5(x.c_str(), x.size());
    // seqan3::debug_stream << "Str: " << x << ". Size: " << x.size() << ". Hash: " << std::hex << hex.p.low << "\n";
    return hex;
};

auto vec2string = [](const auto & x)
{
    std::string res{};
    for (auto i : x)
    {
        res += i.to_char();
    }
    return res;
};

// split function for std::string. see https://marycore.jp/prog/cpp/std-string-split/
template<class T> std::vector<std::string> split(const std::string& s, const T& separator, bool ignore_empty = 0, bool split_empty = 0) {
  struct {
    auto len(const std::string&             s) { return s.length(); }
    auto len(const std::string::value_type* p) { return p ? std::char_traits<std::string::value_type>::length(p) : 0; }
    auto len(const std::string::value_type  c) { return c == std::string::value_type() ? 0 : 1; /*return 1;*/ }
  } util;
  
  if (s.empty()) { /// empty string ///
    if (!split_empty || util.len(separator)) return {""};
    return {};
  }
  
  auto v = std::vector<std::string>();
  auto n = static_cast<std::string::size_type>(util.len(separator));
  if (n == 0) {    /// empty separator ///
    if (!split_empty) return {s};
    for (auto&& c : s) v.emplace_back(1, c);
    return v;
  }
  
  auto p = std::string::size_type(0);
  while (1) {      /// split with separator ///
    auto pos = s.find(separator, p);
    if (pos == std::string::npos) {
      if (ignore_empty && p - n + 1 == s.size()) break;
      v.emplace_back(s.begin() + p, s.end());
      break;
    }
    if (!ignore_empty || p != pos)
      v.emplace_back(s.begin() + p, s.begin() + pos);
    p = pos + n;
  }
  return v;
}

// get read group ID from fastq ID
auto get_rg_id = [](const auto& id, const auto& n_fields, const auto& check)
{
    std::vector<std::string> split_id = split(id, ":");
    std::string rg_id;
    if (check == true && n_fields != split_id.size())
    {
        throw std::runtime_error("Inconsistent ID formats; expected " + std::to_string(n_fields) + " fields, but found " + std::to_string(split_id.size()));
    }
    if (n_fields == 7)
    {
        rg_id = split_id[2] + "." + split_id[3];
    }
    else if (n_fields == 5)
    {
        rg_id = split_id[0] + "." + split_id[1]; 
        // Although the first field is the machine ID in the previous Illumina fastq ID format,
        // we use this field as the run ID is not included in the fastq ID.
    }
    else
    {
        rg_id = "A"s;
    }
    return rg_id;
};

// remove a suffix
auto remove_suffix = [](auto query, const auto& suffix)
{
    if (query.substr(query.length() - suffix.length(), query.length()) == suffix)
    {
        return query.substr(0, query.length() - suffix.length());
    }
    else
    {
        return query;
    }
};

// split by either tab or spaces
auto split_by_spaces = [](auto x)
{
    std::replace(x.begin(), x.end(), '\t', ' ');
    return split(x, ' ', true); 
};

// parse ID
auto parse_ID = [](const auto& ID, const auto& id_index, const auto& suffix)
{
    auto parsed = split_by_spaces(ID)[id_index];
    return remove_suffix(parsed, suffix);
};

// obtain the number of fields in the fastq ID
auto get_n_ID_fields = [](const auto& ID, const auto& index, const auto& suffix)
{
    return (split(parse_ID(ID, index, suffix), ":")).size();
};

// determine the id_index
auto set_id_index = [](const auto& ID, const auto& suffix)
{
    auto split_IDs = split_by_spaces(ID);
    unsigned int n_ID_fields;
    int id_index = 0;
    bool is_id_found = false;
    for (auto split_ID : split_IDs)
    {
        n_ID_fields = split(remove_suffix(split_ID, suffix), ":").size();
        if (n_ID_fields == 7u || n_ID_fields == 5u)
        {
            is_id_found = true;
            break;
        }
        ++id_index;
    }
    if (is_id_found)
    {
        seqan3::debug_stream << "The " << id_index + 1 << " th part of the ID row will be used as a read ID.\n";
    }
    else
    {
        seqan3::debug_stream << "No part of the ID row has 5 or 7 fieids with the colon delimiter. The first part of the ID row will be used as a read ID.\n";
        id_index = 0;
    }
    return id_index;
};

// check whether the format of quality score is phred33 or phred64
auto check_quality_format = [](unsigned int& phred, const std::vector<seqan3::phred94>& qual_char_vec)
{
  std::vector<seqan3::phred94> phred33_specific_chars = {
    '!'_phred94, '"'_phred94, '#'_phred94, '$'_phred94, '%'_phred94, '&'_phred94, '\''_phred94, '('_phred94, ')'_phred94,
    '*'_phred94, '+'_phred94, ','_phred94, '-'_phred94, '.'_phred94, '/'_phred94, '0'_phred94, '1'_phred94, '2'_phred94, 
    '3'_phred94, '4'_phred94, '5'_phred94, '6'_phred94, '7'_phred94, '8'_phred94, '9'_phred94, ':'_phred94};
  std::vector<seqan3::phred94> phred64_suggestive_chars = {
    'K'_phred94, 'L'_phred94, 'M'_phred94, 'N'_phred94, 'O'_phred94, 'P'_phred94, 'Q'_phred94, 'R'_phred94, 'S'_phred94,
    'T'_phred94, 'U'_phred94, 'V'_phred94, 'W'_phred94, 'X'_phred94, 'Y'_phred94, 'Z'_phred94, '['_phred94, '\\'_phred94,
    ']'_phred94, '^'_phred94, '_'_phred94, '`'_phred94, 'a'_phred94, 'b'_phred94, 'c'_phred94, 'd'_phred94, 'e'_phred94,
    'f'_phred94, 'g'_phred94, 'h'_phred94};
  std::vector<seqan3::phred94> solexa_suggestive_chars = {';'_phred94, '<'_phred94, '='_phred94, '>'_phred94, '?'_phred94};
  if (phred == 1 || phred == 2 || phred == 4) {
    throw std::runtime_error("Internal error: Attempt to update the determined phred variable")
  }
  for (auto c : phred33_specific_chars) {
     if (std::find(qual_char_vec.begin(), qual_char_vec.end(), c) != qual_char_vec.end()) {
        phred = 1; // determined
        return;
     };
  }
  for (auto c : phred64_suggestive_chars) {
     if (std::find(qual_char_vec.begin(), qual_char_vec.end(), c) != qual_char_vec.end()) {
        phred = 6; // Phred+64 or Solexa+64. we cannot differentiate Phred+64 from Solexa+64
        break;
     };
  }
  for (auto c : solexa_suggestive_chars) {
     if (std::find(qual_char_vec.begin(), qual_char_vec.end(), c) != qual_char_vec.end()) {
        if (phred == 6) {
            phred = 4; // determined
            return;
        } else {
            phred = 5; // Phred+33 or Solexa+64
            break;
        }
     };
  }
};

//0, undetermined; 1, phred33; 2, phred64; 4, solexa
auto translate_phred_indicator = [](const auto& phred)
{
    std::string phred_char;
    switch (phred) {
		case 0:
			phred_char = "undetermined"s;
			break;
		case 1:
			phred_char = "Phred+33"s;
			break;
		case 2:
			phred_char = "Phred+64"s;
			break;
		case 4:
			phred_char = "Solexa"s;
			break;
		case 5:
			phred_char = "Phred+33 or Solexa+64"s;
			break;
		case 6:
			phred_char = "Phred+64 or Solexa+64"s;
			break;
		default:
			throw std::runtime_error("Internal error: 'phred' must be one of 0,1,2,4,5,6");
			break;
	}
    return phred_char;
};

// scan fastq
// update rg_ids and phred, which indicates whether the fastq format is Illumina 1.8 (Phred+33) or not (Phred+64)
auto scan_fastq = [](auto&& fin, auto& rg_ids, auto& phred, const auto& batch_size, const auto& index, const auto& suffix, const auto& n_ID_fields)
{
    unsigned int n_processed = 0u;
    const bool check = true;
    std::string ID;
    unsigned int scanned_phred = 0; //0, undetermined; 1, phred33; 2, phred64
    std::vector<unsigned int> lengths{};
    for (auto && records : fin | seqan3::views::chunk(batch_size)) // `&&` is important because seqan3::views::chunk returns temporaries!
    {
        for (auto rec : records) // scan all records
        {
            ++n_processed;
            ID = get_rg_id(parse_ID(rec.id(), index, suffix), n_ID_fields, check);
            if (std::find(rg_ids.begin(), rg_ids.end(), ID) == rg_ids.end())
            {
                rg_ids.push_back(ID); // add the record if it was not already added
            }
            if (std::find(lengths.begin(), lengths.end(), rec.sequence().size()) == lengths.end())
            {
                lengths.push_back(rec.sequence().size()); // add the record if it was not already added
            }
            if (scanned_phred != 1 && scanned_phred != 2 && scanned_phred != 4)
            {
                check_quality_format(scanned_phred, rec.base_qualities());
            }
        }
    }

    // report
    seqan3::debug_stream << "Scanned " << n_processed << " records\n";
    seqan3::debug_stream << "Found " << rg_ids.size() << " read group IDs: " << rg_ids << '\n';
    seqan3::debug_stream << "Found " << lengths.size() << " types of the read length: " << lengths << '\n';

    // check consistency between scanned_phred and phred (inferred from the seq ID)
    if (phred == 0)
    {
        if (scanned_phred == 0)
        {
            throw std::runtime_error("The format of quality score couldn't be determined from the fastq. Please contact the developer, as the option to specify the format from the command line has not been implemented yet.");
        }
        phred = scanned_phred;
        seqan3::debug_stream << "The format of quality score was " << translate_phred_indicator(phred) << '\n';
    }
    else if (phred != scanned_phred && scanned_phred != 0)
    {
        seqan3::debug_stream << "The format of quality score was inconsistent between IDs (" << translate_phred_indicator(phred) 
                             << ") and scores themselves (" << translate_phred_indicator(scanned_phred) << ")\nWill use that from scores themselves\n";
        phred = scanned_phred; 
    }
    else 
    {
        seqan3::debug_stream << "The format of quality score was " << translate_phred_indicator(phred) << '\n';
    }
};

// convert phred68solexa to phred42
// TODO
std::vector<seqan3::phred42> convert_phred_score(const std::vector<seqan3::phred42>& qual_vec, const unsigned int& phred)
{
    std::vector<seqan3::phred42> res(qual_vec.size());
    seqan3::phred42 q;
    for (unsigned int i=0 ; i<qual_vec.size() ; i++)
    {
        q.assign_phred(qual_vec.at(i).to_phred());
        res.at(i) = q;
    }
    return res;
}
std::vector<seqan3::phred42> convert_phred_score(const std::vector<seqan3::phred68solexa>& qual_vec, const unsigned int& phred)
{
    std::vector<seqan3::phred42> res(qual_vec.size());
    seqan3::phred42 q;
    for (unsigned int i=0 ; i<qual_vec.size() ; i++)
    {
        q.assign_phred(qual_vec.at(i).to_phred());
        res.at(i) = q;
    }
    return res;
}

// typedef sam_record_type
using types = seqan3::type_list<std::string, std::vector<seqan3::dna5>, std::vector<seqan3::phred42>, seqan3::sam_flag, seqan3::sam_tag_dictionary>;
using fields = seqan3::fields<seqan3::field::id, seqan3::field::seq, seqan3::field::qual, seqan3::field::flag, seqan3::field::tags>;
using sam_record_type = seqan3::sam_record<types, fields>;


sam_record_type make_sam_record(const auto& ID, const auto& seq, const auto& flag, const auto& dict, const std::vector<seqan3::phred42>& qual)
{
    return sam_record_type{ID, seq, qual, flag, dict};
}

sam_record_type make_sam_record(const auto& ID, const auto& seq, const auto& flag, const auto& dict, const std::vector<seqan3::phred68solexa>& qual)
{
    return sam_record_type{ID, seq, convert_phred_score(qual), flag, dict};
}

// main function for converting fastq files to a sam file
auto fastq2sam_iter = [] (auto&& fin1, auto&& fin2, auto&& fout, const auto& hash, const auto& hash_no_quality,
                          const auto& batch_size, const auto& min_length,
                          const auto& id_index, const auto& suffix1, const auto& suffix2, const std::vector<seqan3::sam_tag_dictionary>& dicts, 
                          const auto& n_ID_fields, const auto& phred)
{
    const bool check = false;
    // sam_flag
    auto first_flag = seqan3::sam_flag::paired | seqan3::sam_flag::unmapped | seqan3::sam_flag::mate_unmapped | seqan3::sam_flag::first_in_pair;
    auto second_flag = seqan3::sam_flag::paired | seqan3::sam_flag::unmapped | seqan3::sam_flag::mate_unmapped | seqan3::sam_flag::second_in_pair;

    bool error_flag;
    std::string ID1;
    std::string ID2;
    std::string rg_id;
    seqan3::sam_tag_dictionary dict;
    std::vector<sam_record_type> outs(batch_size * 2, sam_record_type{});
    unsigned int n_processed = 0u;
    unsigned int n_skipped = 0u;
    uint64_t sum = 0;
    uint64_t sum_no_quality = 0;
    // iterate over records
    for (auto && [record1, record2] : seqan3::views::zip(fin1 | seqan3::views::chunk(batch_size), fin2 | seqan3::views::chunk(batch_size))) // && is important!
    {                                                           // because seqan3::views::zip returns temporaries
        unsigned int i = 0u;
        outs.resize(batch_size * 2);
        for (auto && [rec1, rec2] : seqan3::views::zip(record1, record2)) 
        {
            ++n_processed;
            // seqan3::debug_stream << "ID:  " << rec1.id() << '\n';
            // parse id
            ID1 = parse_ID(rec1.id(), id_index, suffix1);
            ID2 = parse_ID(rec2.id(), id_index, suffix2);
            if (ID1 != ID2) {
                throw std::runtime_error("Your pairs don't match. ID in the file 1, " + rec1.id() + "; ID in the file 2, " + rec2.id());
            }
            
            rg_id = get_rg_id(ID1, n_ID_fields, check);
            error_flag = true;
            for (auto d : dicts) {
                if (rg_id == d.get<"RG"_tag>()) {
                    dict = d;
                    error_flag = false;
                    break;
                }
            }
            if (error_flag == true) {
                throw std::runtime_error("In the " + std::to_string(n_processed) + " th records, the read group ID is not included in the dict: " + rg_id +
                                         " (fastq ID: " + rec1.id() + ")");
            }
            if (rec1.sequence().size() < min_length || rec2.sequence().size() < min_length) {
                ++n_skipped;
            } else {
                // store records
                auto bq1 = convert_phred_score(rec1.base_qualities(), phred);
                auto bq2 = convert_phred_score(rec2.base_qualities(), phred);
                outs.at(i * 2) = sam_record_type{ID1, rec1.sequence(), bq1, first_flag, dict};
                outs.at(i * 2 + 1) = sam_record_type{ID2, rec2.sequence(), bq2, second_flag, dict};

                // hash
                if (!hash.empty())
                {
                    hexSum(calc_hash_str(ID1 + "/1" + vec2string(rec1.sequence()) + vec2string(bq1)), sum);
                    hexSum(calc_hash_str(ID2 + "/2" + vec2string(rec2.sequence()) + vec2string(bq2)), sum);
                }
                if (!hash_no_quality.empty())
                {
                    hexSum(calc_hash_str((ID1 + "/1" + vec2string(rec1.sequence()))), sum_no_quality);
                    hexSum(calc_hash_str((ID2 + "/2" + vec2string(rec2.sequence()))), sum_no_quality);
                }
                ++i;
            }

        }
        outs.resize(i * 2); // to shrink the vector for the final iteration and skipped records
        fout = outs; //write out
    }
    if (!hash.empty())
    {
        std::ofstream writing_file;
        writing_file.open(hash);
        writing_file << std::hex << sum << "\t";
        writing_file << std::dec << n_processed << "\n";
        writing_file.close();
    }
    if (!hash_no_quality.empty())
    {
        std::ofstream writing_file;
        writing_file.open(hash_no_quality);
        writing_file << std::hex << sum_no_quality << "\t";
        writing_file << std::dec << n_processed << "\n";
        writing_file.close();
    }
    seqan3::debug_stream << "Done. Processed " << n_processed << " record pairs and skipped " << n_skipped << " pairs\n";
};

// fastq2sam
void fastq2sam(std::filesystem::path & fastq1, std::filesystem::path & fastq2, std::filesystem::path & out, std::string & sample_name, 
               std::filesystem::path & hash, std::filesystem::path & hash_no_quality,
               std::string & platform, unsigned int & min_length, int & id_index, std::string & suffix1, std::string & suffix2, 
               unsigned int & batch_size, bool & no_pg)
{
    // iterate over fastq1 to obtain read group IDs (in the {run.lane} format)
    // get_n_ID_fields
    std::vector<std::string> rg_ids;
    unsigned int phred; 
    seqan3::debug_stream << "Scanning the first fastq\n";
    sequence_file_input_phred94 fin{fastq1};
    auto rec = *fin.begin();
    if (id_index < 0)
    {
        id_index = set_id_index(rec.id(), suffix1);
    }
    unsigned int n_ID_fields = get_n_ID_fields(rec.id(), id_index, suffix1);
    // msg and get rg from the first record 
    if (n_ID_fields == 7u)
    {
        seqan3::debug_stream << "The estimated ID format is CASAVA-1.8\n";
        rg_ids = {get_rg_id(parse_ID(rec.id(), id_index, suffix1), n_ID_fields, true)};
        phred = 1;
    }
    else if (n_ID_fields == 5u)
    {
        seqan3::debug_stream << "The estimated ID format is the previous format used before CASAVA-1.8\n";
        rg_ids = {get_rg_id(parse_ID(rec.id(), id_index, suffix1), n_ID_fields, true)};
        phred = 6;
    }
    else
    {
        seqan3::debug_stream << "The number of the fields in the fastq ID is " << n_ID_fields << ", not 5 or 7.\n";
        seqan3::debug_stream << "The read group ID will be always 'A'\n";
        rg_ids = {"A"s};
        phred = 0;
    }
    // iterate over the rest of fastq1 to obtain read group IDs (in the {run.lane} format)
    scan_fastq(fin, rg_ids, phred, batch_size, id_index, suffix1, n_ID_fields); //update rg_ids and phred
   
    if (out.empty())
    {
        return; // scan only
    }
    
    // Prep fout and dicts
    std::vector<std::string> ref_ids{"dummy_ref"};
    std::vector<size_t> ref_lengths{1};
    seqan3::sam_file_output fout{out, ref_ids, ref_lengths};

    // header: program_infos
    if (no_pg == false)
    {
        seqan3::sam_file_program_info_t pg{};
        pg.command_line_call = 
            PG + " --fastq1 "s + fastq1.string() + " --fastq2 "s + fastq2.string() + " --out "s + out.string() + " --sample-name "s + sample_name + 
            " --platform "s + platform + " --min-length "s + std::to_string(min_length) + " --id-index "s + std::to_string(id_index) + 
            " --suffix1 "s + suffix1 + " --suffix2 "s + suffix2 + " --batch-size "s + std::to_string(batch_size); 
        pg.id = PG;
        pg.name = PG;
        pg.version = VER;
        fout.header().program_infos.push_back(pg);
    }
    
    // header: order
    fout.header().grouping = "none"s;
    fout.header().sorting = "unsorted"s;
    fout.header().subsorting = "queryname"s;

    // read_group & sam_tag
    // typedef read_group
    using read_group = std::pair<std::string, std::string>; //for header
    read_group rg; 
    seqan3::sam_tag_dictionary dict{}; // initialise empty dictionary
    std::vector<seqan3::sam_tag_dictionary> dicts{};
    for (auto i : rg_ids)
    {
        rg = {i, "SM:"s + sample_name + "\tLB:"s + sample_name + "\tPL:"s + platform};
        fout.header().read_groups.push_back(rg);        
        dict.get<"RG"_tag>() = i;
        dicts.push_back(dict);
    }

    // iterate over records
    seqan3::debug_stream << "Writing the sam/bam file\n";
    if (phred == 1) {
        seqan3::sequence_file_input fin1{fastq1};
        seqan3::sequence_file_input fin2{fastq2};
        fastq2sam_iter(fin1, fin2, fout, hash, hash_no_quality, batch_size, min_length, id_index, suffix1, suffix2, dicts, n_ID_fields, phred);
    } else {
        sequence_file_input_phred64 fin1{fastq1};
        sequence_file_input_phred64 fin2{fastq2};
        fastq2sam_iter(fin1, fin2, fout, hash, hash_no_quality, batch_size, min_length, id_index, suffix1, suffix2, dicts, n_ID_fields, phred);
    }
}

int main(int argc, char ** argv)
{
    sharg::parser myparser{"fastq2sam", argc, argv}; // initialise myparser
    cmd_arguments args{};
    initialise_parser(myparser, args); 
    try
    {
        myparser.parse(); // trigger command line parsing
    }
    catch (sharg::parser_error const & ext) // catch user errors
    {
        std::cerr << "[ERROR] " << ext.what() << "\n"; // error message
        return EXIT_FAILURE;
    }

    seqan3::contrib::bgzf_thread_count = args.bam_writers;
    fastq2sam(args.fastq1, args.fastq2, args.out, args.sample_name, args.hash, args.hash_no_quality,
              args.platform, args.min_length, args.id_index, args.suffix1, args.suffix2, args.batch_size, args.no_pg);
    return EXIT_SUCCESS;
}


