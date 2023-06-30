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
#include <sharg/all.hpp> // argparser for seqan3.3
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
    unsigned int min_length{0u};
    unsigned short id_index{0u};
    std::string suffix1{"/1"};
    std::string suffix2{"/2"};
    unsigned int batch_size{1000000u};
    unsigned int bam_writers{4u};
    bool no_pg{false};
};

void initialise_parser(sharg::parser & parser, cmd_arguments & args)
{
    parser.info.author = "Shinichi Namba";
    parser.info.short_description = "Converting paired-end fastq files with illumina-format quality scores into an unmapped sam/bam file while properly handling the RG tags.";
    parser.info.version = VER;
    parser.info.short_copyright = "Internal only";
    parser.info.long_copyright = "Users can freely use and modify this program. However, please always ask Shinichi Namba for the re-distribution.";

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
                                    .description = "The output sam/bam file. If not set, this program will only scan the first fastq file.", 
                                    .required = false,
                                    .validator = sharg::output_file_validator{sharg::output_file_open_options::open_or_create, {"sam", "bam"}}});
    parser.add_option(args.sample_name,
                      sharg::config{.short_id = 'n',
                                    .long_id = "sample-name",
                                    .description = "Sample ID.", 
                                    .required = true});
    parser.add_option(args.min_length,
                      sharg::config{.short_id = 'm',
                                    .long_id = "min-length",
                                    .description = "The minimum read length to be written to the output file.", 
                                    .advanced = true});
    parser.add_option(args.id_index,
                      sharg::config{.short_id = 'i',
                                    .long_id = "id-index",
                                    .description = "The column index of ID to be parsed (0-based).", 
                                    .advanced = true,
                                    .validator = sharg::arithmetic_range_validator{0, 10}});
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

// parse ID
auto parse_ID = [](auto ID, const auto& id_index, const auto& suffix)
{
    std::replace(ID.begin(), ID.end(), '\t', ' ');
    ID = split(ID, ' ', true)[id_index]; 
    ID = remove_suffix(ID, suffix);
    return ID;
};

// obtain the number of fields in the fastq ID
auto check_n_ID_fields = [](const auto& ID, const auto& index, const auto& suffix)
{
    auto n_ID_fields = (split(parse_ID(ID, index, suffix), ":")).size();
    if (n_ID_fields == 7u)
    {
        seqan3::debug_stream << "The estimated ID format is CASAVA-1.8\n";
    }
    else if (n_ID_fields == 5u)
    {
        seqan3::debug_stream << "The estimated ID format is the previous format used before CASAVA-1.8\n";
    }
    else
    {
        seqan3::debug_stream << "The number of the fields in the fastq ID is " << n_ID_fields << ", not 5 or 7.\n";
    }
    return n_ID_fields;
};

// scan fastq and obtain rg_ids
auto scan_fastq = [](auto&& fin, const auto& batch_size, const auto& index, const auto& suffix, const auto& n_ID_fields)
{
    unsigned int n_processed = 0u;
    bool check = true;
    std::string ID;
    std::vector<std::string> rg_ids{};
    std::vector<unsigned int> lengths{};
    for (auto && records : fin | seqan3::views::chunk(batch_size)) // `&&` is important because seqan3::views::chunk returns temporaries!
    {
        for (auto rec : records) // scan all records
        {
            ++n_processed;
            ID = parse_ID(rec.id(), index, suffix);
            ID = get_rg_id(ID, n_ID_fields, check);
            if (std::find(rg_ids.begin(), rg_ids.end(), ID) == rg_ids.end())
            {
                rg_ids.push_back(ID); // add the record if it was not already added
            }
            if (std::find(lengths.begin(), lengths.end(), rec.sequence().size()) == lengths.end())
            {
                lengths.push_back(rec.sequence().size()); // add the record if it was not already added
            }
            
        }
    }    
    // report
    seqan3::debug_stream << "Scanned " << n_processed << " records\n";
    seqan3::debug_stream << "Found " << rg_ids.size() << " read group IDs: " << rg_ids << '\n';
    seqan3::debug_stream << "Found " << lengths.size() << " types of the read length: " << lengths << '\n';
    return rg_ids;
};

void fastq2sam(std::filesystem::path & fastq1, std::filesystem::path & fastq2, std::filesystem::path & out, std::string sample_name, 
               unsigned int min_length, unsigned short id_index, std::string suffix1, std::string suffix2, unsigned int batch_size, bool no_pg)
{
    // iterate over fastq1 to obtain read group IDs (in the {run.lane} format)
    std::vector<std::string> rg_ids;
    seqan3::debug_stream << "Scanning the first fastq\n";
    seqan3::sequence_file_input fin{fastq1};
    
    auto rec = *fin.begin();
    unsigned int n_ID_fields = check_n_ID_fields(rec.id(), id_index, suffix1);
    if (n_ID_fields != 7u && n_ID_fields != 5u && !out.empty())
    {
        seqan3::debug_stream << "The read group ID will be always 'A'\n";
        seqan3::debug_stream << "The read length: " << rec.sequence().size() << "\n";
        rg_ids = {"A"s};
    }
    else
    {
        rg_ids = scan_fastq(fin, batch_size, id_index, suffix1, n_ID_fields);
    }
   
    if (out.empty())
    {
        return; // scan only
    }
    // typedef sam_record_type and read_group
    using types = seqan3::type_list<std::string, std::vector<seqan3::dna5>, std::vector<seqan3::phred42>, seqan3::sam_flag, seqan3::sam_tag_dictionary>;
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
    if (no_pg == false)
    {
        seqan3::sam_file_program_info_t pg{};
        pg.command_line_call = 
            PG + " --fastq1 " + fastq1.string() + " --fastq2 " + fastq2.string() + " --out " + out.string() + " --sample-name " + sample_name + 
            " --min-length " + std::to_string(min_length) + " --id-index " + std::to_string(id_index) + " --suffix1 " + suffix1 + " --suffix2 " + suffix2 +
            " --batch-size " + std::to_string(batch_size); 
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
    std::string ID1;
    std::string ID2;
    std::string rg_id;
    seqan3::debug_stream << "Writing the sam/bam file\n";
    std::vector<sam_record_type> outs(batch_size * 2, sam_record_type{});
    unsigned int n_processed = 0u;
    unsigned int n_skipped = 0u;
    bool check = false;
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
                throw std::runtime_error("In the " + std::to_string(n_processed) + " th records, the read group ID is not included in the dict: " + rg_id +
                                         " (fastq ID: " + rec1.id() + ")");
            }
            if (rec1.sequence().size() < min_length || rec2.sequence().size() < min_length)
            {
                ++n_skipped;
            }
            else
            {
                // store records
                outs.at(i * 2) = sam_record_type{ID1, rec1.sequence(), rec1.base_qualities(), first_flag, dict};
                outs.at(i * 2 + 1) = sam_record_type{ID2, rec2.sequence(), rec2.base_qualities(), second_flag, dict};   
                ++i;
            }
        }
        outs.resize(i * 2); // to shrink the vector for the final iteration and skipped records
        fout = outs; //write out
    }
    seqan3::debug_stream << "Done. Processed " << n_processed << " record pairs and skipped " << n_skipped << " pairs\n";
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
    fastq2sam(args.fastq1, args.fastq2, args.out, args.sample_name, args.min_length, args.id_index, args.suffix1, args.suffix2, args.batch_size, args.no_pg);
    return EXIT_SUCCESS;
}
