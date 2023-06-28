#include <cstdlib> //EXIT_SUCCESS, EXIT_FAILURE
#include <filesystem>
#include <string>
#include <algorithm>  // std::unique
#include <vector>
#include <utility>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/utility/views/all.hpp>
#include <seqan3/io/sam_file/all.hpp>
#include <sharg/all.hpp> // argparser
using namespace std::literals; 
using namespace seqan3::literals; 

// program name and version
const std::string PG{"fastq2sam"s};
const std::string VER{"0.0.1"s};

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
auto get_rg_id_illumina1_8 = [](auto id)
{
    std::vector<std::string> split_id = split(id, ":");
    std::string rg_id = split_id[2] + "." + split_id[3];
    return rg_id;
};

// default values of the arguments
struct cmd_arguments
{
    std::filesystem::path fastq1{};
    std::filesystem::path fastq2{};
    std::filesystem::path out{};
    std::string sample_name{};
    int batch_size{2}; // TODO 10000
};

void initialise_parser(sharg::parser & parser, cmd_arguments & args)
{
    parser.info.author = "Shinichi Namba";
    parser.info.short_description = "Converting paired-end fastq files with illumina-format quality scores into an unmapped sam/bam file while properly handling the RG tags.";
    parser.info.version = VER;

    parser.add_option(args.fastq1,
                      sharg::config{.short_id = '1',
                                    .long_id = "fastq1",
                                    .description = "The first file of the paired-end fastq files.", 
                                    .required = true,
                                    .validator = sharg::input_file_validator{{"fq", "fastq", "fq.gz", "fastq.gz"}}});
    parser.add_option(args.fastq2,
                      sharg::config{.short_id = '2',
                                    .long_id = "fastq2",
                                    .description = "The second file of the paired-end fastq files.", 
                                    .required = true,
                                    .validator = sharg::input_file_validator{{"fq", "fastq", "fq.gz", "fastq.gz"}}});
    parser.add_option(args.out,
                      sharg::config{.short_id = 'o',
                                    .long_id = "out",
                                    .description = "The output sam/bam file.", 
                                    .required = true,
                                    .validator = sharg::output_file_validator{sharg::output_file_open_options::open_or_create, {"sam", "bam"}}});
    parser.add_option(args.sample_name,
                      sharg::config{.short_id = 'n',
                                    .long_id = "sample-name",
                                    .description = "Sample ID.", 
                                    .required = true});
    parser.add_option(args.batch_size,
                      sharg::config{.short_id = 'b',
                                    .long_id = "batch-size",
                                    .description = "The batch size to read fastq files at a time.", 
                                    .required = false,
                                    .validator = sharg::arithmetic_range_validator{1, 1000000}});
}

void fastq2sam(std::filesystem::path & fastq1, std::filesystem::path & fastq2, std::filesystem::path & out, std::string sample_name, int batch_size)
{
    // iterate over fastq1 to obtain read group IDs (in the {run.lane} format)
    seqan3::debug_stream << "Obtaining read group IDs from the first fastq" << '\n';
    std::vector<std::string> rg_ids{};
    std::string rg_id;
    seqan3::sequence_file_input fin{fastq1};
    for (auto && records : fin | seqan3::views::chunk(batch_size)) // `&&` is important because seqan3::views::chunk returns temporaries!
    {
        // `records` contains batch_size elements (or less at the end)
        rg_id = get_rg_id_illumina1_8((*records.begin()).id());
        rg_ids.push_back(rg_id); // store the first ID in batch
    }
    // sort&unique
    std::sort(rg_ids.begin(), rg_ids.end());
    rg_ids.erase(std::unique(rg_ids.begin(), rg_ids.end()), rg_ids.end());
    seqan3::debug_stream << "Read group IDs: " << rg_ids << '\n';
    
    // def sam_record_type
    using types = seqan3::type_list<std::string &, std::vector<seqan3::dna5> &, std::vector<seqan3::phred42> &, seqan3::sam_flag, seqan3::sam_tag_dictionary>;
    using fields = seqan3::fields<seqan3::field::id, seqan3::field::seq, seqan3::field::qual, seqan3::field::flag, seqan3::field::tags>;
    using sam_record_type = seqan3::sam_record<types, fields>;
    using read_group = std::pair<std::string, std::string>; //for header

    // sam_flag
    auto first_flag = seqan3::sam_flag::paired | seqan3::sam_flag::unmapped | seqan3::sam_flag::mate_unmapped | seqan3::sam_flag::first_in_pair;
    auto second_flag = seqan3::sam_flag::paired | seqan3::sam_flag::unmapped | seqan3::sam_flag::mate_unmapped | seqan3::sam_flag::second_in_pair;
    
    // Output
    std::vector<std::string> ref_ids{"dummy_ref"};
    std::vector<size_t> ref_lengths{1};
    seqan3::sam_file_output fout{out, ref_ids, ref_lengths};

    // header: program_infos
    seqan3::sam_file_program_info_t pg{};
    pg.command_line_call = PG + " -1 " + fastq1.string() + " -2 " + fastq2.string() + " -o " + out.string() + " -n " + sample_name + " -b " + std::to_string(batch_size);
    pg.id = PG;
    pg.name = PG;
    pg.version = VER;
    fout.header().program_infos.push_back(pg);
    
    // header: order
    fout.header().grouping = "none"s;
    fout.header().sorting = "unsorted"s;
    fout.header().subsorting = "queryname"s;

    // read_group & sam_tag
    read_group rg; 
    seqan3::sam_tag_dictionary dict{}; // initialise empty dictionary
    std::vector<seqan3::sam_tag_dictionary> dicts{};
    for (auto i : rg_ids)
    {
        rg = {i, "SM:" + sample_name + "\tLB:" + sample_name + "\tPL:ILLUMINA"};
        fout.header().read_groups.push_back(rg);        
        dict.get<"RG"_tag>() = i;
        dicts.push_back(dict);
    }

    // iterate over records
    // Input for simplicity we take the same file
    seqan3::sequence_file_input fin1{fastq1};
    seqan3::sequence_file_input fin2{fastq2};
    bool error_flag;
    seqan3::debug_stream << "Writing the sam/bam file: " << '\n';
    for (auto && [record1, record2] : seqan3::views::zip(fin1 | seqan3::views::chunk(batch_size), fin2 | seqan3::views::chunk(batch_size))) // && is important!
    {                                                           // because seqan3::views::zip returns temporaries
        for (auto && [rec1, rec2] : seqan3::views::zip(record1, record2)) 
        {
            seqan3::debug_stream << "ID:  " << rec1.id() << '\n';
            // TODO parse id
            if (rec1.id() != rec2.id()) { // TODO
                throw std::runtime_error("Your pairs don't match."); // TODO more details
            }
            
            rg_id = get_rg_id_illumina1_8(rec1.id());
            error_flag = true;
            for (auto d : dicts)
            {
                if (rg_id == d.get<"RG"_tag>())
                {
                    dict = d;
                    error_flag = false;
                    break;
                }
            }
            if (error_flag == true)
            {
                throw std::runtime_error("The read group ID is not included in the dict."); // TODO more details
            }
            std::vector<sam_record_type> range{
                sam_record_type{rec1.id(), rec1.sequence(), rec1.base_qualities(), first_flag, dict}, 
                sam_record_type{rec2.id(), rec2.sequence(), rec2.base_qualities(), second_flag, dict}
            };
            fout = range; //write out 2 lines
        }
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
        std::cerr << "[ERROR] " << ext.what() << "\n"; // customise your error message
        return EXIT_FAILURE;
    }

    fastq2sam(args.fastq1, args.fastq2, args.out, args.sample_name, args.batch_size);
    return EXIT_SUCCESS;
}