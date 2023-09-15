#include <cstdlib> //EXIT_SUCCESS, EXIT_FAILURE
#include <cstddef> // std::size_t
#include <filesystem>
#include <string>
#include <vector>
#include <utility>
#include <sharg/all.hpp> // argparser for seqan3.3
#include <seqan3/contrib/stream/bgzf.hpp> // seqan3::contrib::bgzf_thread_count

#include <seqan3/core/debug_stream.hpp>
#include "sequence_file_input.h"
#include "phred.h"
#include "scan_fastq.h"
// #include "to_sam.h"
using namespace std::literals; 
using namespace seqan3::literals; 

/*
phred test
    phred a(0B1110);
    phred b(0B0011);
    phred c(0B1101);
    std::cout << "Hello World! " << a.to_string() << "\n";
    a &= b;
    std::cout << "Hello World! " << a.to_string() << "\n";
    a.update(c);
    std::cout << "Hello World! " << a.to_string() << "\n";
*/

// we cannot separate the declaration of to_sam, because we define the class of the output object in the main program
#include <seqan3/io/sam_file/output.hpp>
#include <seqan3/io/sam_file/record.hpp>
#include "fastq_metadata.h"

#include <stdint.h> // uint64_t
#include <iostream>
#include <fstream> // write out hash
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/utility/views/chunk.hpp>
#include <seqan3/utility/views/zip.hpp>
#include <seqan3/alphabet/quality/phred94.hpp> // seqan3::phred94
#include "convert_phred_score.h"
#include "id_parser.h"
#include "bamhash.h"

// sam_flag
auto first_flag = seqan3::sam_flag::paired | seqan3::sam_flag::unmapped | seqan3::sam_flag::mate_unmapped | seqan3::sam_flag::first_in_pair;
auto second_flag = seqan3::sam_flag::paired | seqan3::sam_flag::unmapped | seqan3::sam_flag::mate_unmapped | seqan3::sam_flag::second_in_pair;

// typedef sam_record_type
using types = seqan3::type_list<std::string, std::vector<seqan3::dna5>, std::vector<seqan3::phred94>, seqan3::sam_flag, seqan3::sam_tag_dictionary>;
using fields = seqan3::fields<seqan3::field::id, seqan3::field::seq, seqan3::field::qual, seqan3::field::flag, seqan3::field::tags>;
using sam_record_type = seqan3::sam_record<types, fields>;

// convert fastq files to a sam file
template <class T>
void to_sam(fastq_metadata metadata, T& fin1, T& fin2, auto& fout, const std::vector<seqan3::sam_tag_dictionary>& dicts,
            std::filesystem::path hash, std::filesystem::path hash_no_quality, 
            const std::size_t& batch_size, const std::size_t& min_length, const std::string& suffix1, const std::string& suffix2) {
    if (!metadata.enough_info()) {
        throw std::runtime_error("The metadata does not contain enough information.");
    }
    const bool check = false;
    bool error_flag;
    std::string ID1;
    std::string ID2;
    std::string rg_id;
    seqan3::sam_tag_dictionary dict;
    std::vector<sam_record_type> outs(batch_size * 2, sam_record_type{});
    std::size_t n_processed = 0u;
    std::size_t n_skipped = 0u;
    uint64_t sum = 0;
    uint64_t sum_no_quality = 0;
    // iterate over records
    for (auto && [record1, record2] : seqan3::views::zip(fin1 | seqan3::views::chunk(batch_size), fin2 | seqan3::views::chunk(batch_size))) {// && is important!
        std::size_t i = 0u;
        outs.resize(batch_size * 2);
        for (auto && [rec1, rec2] : seqan3::views::zip(record1, record2)) {
            ++n_processed;
            // parse id
            ID1 = parse_ID(rec1.id(), metadata.id_index, suffix1);
            ID2 = parse_ID(rec2.id(), metadata.id_index, suffix2);
            if (ID1 != ID2) {
                throw std::runtime_error("Your pairs don't match. ID in the file 1, " + rec1.id() + "; ID in the file 2, " + rec2.id());
            }
            
            rg_id = get_rg_id(ID1, metadata.n_ID_fields, check);
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
                auto bq1 = convert_phred_score(rec1.base_qualities(), metadata.format);
                auto bq2 = convert_phred_score(rec2.base_qualities(), metadata.format);
                outs.at(i * 2) = sam_record_type{ID1, rec1.sequence(), bq1, first_flag, dict};
                outs.at(i * 2 + 1) = sam_record_type{ID2, rec2.sequence(), bq2, second_flag, dict};

                // hash
                if (!hash.empty()) {
                    hexSum(calc_hash_str(ID1 + "/1" + vec2string(rec1.sequence()) + vec2string(bq1)), sum);
                    hexSum(calc_hash_str(ID2 + "/2" + vec2string(rec2.sequence()) + vec2string(bq2)), sum);
                }
                if (!hash_no_quality.empty()) {
                    hexSum(calc_hash_str((ID1 + "/1" + vec2string(rec1.sequence()))), sum_no_quality);
                    hexSum(calc_hash_str((ID2 + "/2" + vec2string(rec2.sequence()))), sum_no_quality);
                }
                ++i;
            }
        }
        outs.resize(i * 2); // to shrink the vector for the final iteration and skipped records
        fout = outs; //write out
    }
    if (!hash.empty()) {
        std::ofstream writing_file;
        writing_file.open(hash);
        writing_file << std::setw(16) << std::setfill('0') << std::hex << sum;
        writing_file << "\t" << std::dec << n_processed << "\n";
        writing_file.close();
    }
    if (!hash_no_quality.empty()) {
        std::ofstream writing_file;
        writing_file.open(hash_no_quality);
        writing_file << std::setw(16) << std::setfill('0') << std::hex << sum_no_quality;
        writing_file << "\t" << std::dec << n_processed << "\n";
        writing_file.close();
    }
   std::cerr << "Done. Processed " << n_processed << " record pairs and skipped " << n_skipped << " pairs\n";
};
// end of to_sam()


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
    std::string library{};
    std::filesystem::path hash{};
    std::filesystem::path hash_no_quality{};
    std::string platform{"ILLUMINA"s};
    std::size_t min_length{0u};
    int id_index{-1};
    std::string suffix1{"/1"s};
    std::string suffix2{"/2"s};
    std::size_t batch_size{1000000u};
    std::size_t bam_writers{4u};
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
                                    .validator = sharg::output_file_validator{sharg::output_file_open_options::open_or_create, {"sam", "bam"}}});
    parser.add_option(args.sample_name,
                      sharg::config{.short_id = 'n',
                                    .long_id = "sample-name",
                                    .description = "Sample ID. Required if the output file is set."});
    parser.add_option(args.library,
                      sharg::config{.short_id = 'l',
                                    .long_id = "library",
                                    .description = "Library ID. If not set, the sample ID will be used."});
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

// fastq2sam
void fastq2sam(const cmd_arguments& args, const fastq_metadata metadata) {
    // Prep fout and dicts
    std::vector<std::string> ref_ids{"dummy_ref"};
    std::vector<std::size_t> ref_lengths{1};
    seqan3::sam_file_output fout(args.out, ref_ids, ref_lengths);

    // header: program_infos
    if (args.no_pg == false) {
        seqan3::sam_file_program_info_t pg{};
        pg.command_line_call = 
            PG + " --fastq1 "s + args.fastq1.string() + " --fastq2 "s + args.fastq2.string() + " --out "s + args.out.string() + 
            " --sample-name "s + args.sample_name + " --library "s + args.library + " --platform "s + args.platform + 
            " --min-length "s + std::to_string(args.min_length) + " --id-index "s + std::to_string(args.id_index) + 
            " --suffix1 "s + args.suffix1 + " --suffix2 "s + args.suffix2 + " --batch-size "s + std::to_string(args.batch_size); 
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
    for (auto i : metadata.rg_ids) {
        rg = {i, "SM:"s + args.sample_name + "\tLB:"s + args.library + "\tPL:"s + args.platform};
        fout.header().read_groups.push_back(rg);        
        dict.get<"RG"_tag>() = i;
        dicts.push_back(dict);
    }

    // iterate over records
    seqan3::debug_stream << "Writing the sam/bam file\n";
    if (metadata.format == phred{0B0100} || metadata.format == phred{0B1000}) {
        sequence_file_input_phred94 fin1{args.fastq1};
        sequence_file_input_phred94 fin2{args.fastq2};
        to_sam(metadata, fin1, fin2, fout, dicts, args.hash, args.hash_no_quality, args.batch_size, args.min_length, args.suffix1, args.suffix2);
    } else {
        sequence_file_input_phred68solexa fin1{args.fastq1};
        sequence_file_input_phred68solexa fin2{args.fastq2};
        to_sam(metadata, fin1, fin2, fout, dicts, args.hash, args.hash_no_quality, args.batch_size, args.min_length, args.suffix1, args.suffix2);
    }
}

int main(int argc, char ** argv)
{
    // parse args
    sharg::parser myparser{"fastq2sam", argc, argv, sharg::update_notifications::off}; // initialise myparser
    cmd_arguments args{};
    initialise_parser(myparser, args); 
    try {
        myparser.parse(); // trigger command line parsing
    } catch (sharg::parser_error const & ext) { // catch user errors
        std::cerr << "[ERROR] " << ext.what() << "\n"; // error message
        return EXIT_FAILURE;
    }
    if (args.sample_name.empty() && !args.out.empty()) {
        std::runtime_error("`--sample-name` must be provided if `--out` is not empty."s);
    }
    if (args.library.empty()) {
        args.library = args.sample_name;
    }

    // scan
    fastq_metadata metadata = scan_fastq(args.fastq1, args.fastq2, args.batch_size, args.id_index, args.suffix1, args.suffix2);
    metadata.print();
    if (args.out.empty()) {
        return EXIT_SUCCESS; // scan only
    }

    // fastq2sam
    seqan3::contrib::bgzf_thread_count = args.bam_writers;
    fastq2sam(args, metadata);

    return EXIT_SUCCESS;
}
