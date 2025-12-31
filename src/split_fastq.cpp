#include <cstdlib> //EXIT_SUCCESS, EXIT_FAILURE
#include <cstddef> // std::size_t
#include <filesystem>
#include <string>
#include <vector>
#include <array>
#include <utility>
#include <sharg/all.hpp> // argparser for seqan3.3
// #include <seqan3/contrib/stream/bgzf.hpp> // seqan3::contrib::bgzf_thread_count
#include <format> // for zero padding (requiring c++20)

#include <seqan3/core/debug_stream.hpp>
#include "sequence_file_input.h"
#include "phred.h"
#include "scan_fastq.h"
using namespace std::literals; 
using namespace seqan3::literals; 

#include "fastq_metadata.h"

#include <stdint.h> // uint64_t
#include <iostream>
#include <fstream> // write out hash
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/utility/views/zip.hpp>
#include <seqan3/alphabet/quality/phred94.hpp> // seqan3::phred94
#include "convert_phred_score.h"
#include "id_parser.h"
#include "bamhash.h"
#include "version.h"
// typedef fq_record_type
using types = seqan3::type_list<std::string, std::vector<seqan3::dna5>, std::vector<seqan3::phred94>>;
using fields = seqan3::fields<seqan3::field::id, seqan3::field::seq, seqan3::field::qual>;
using fq_record_type = seqan3::sequence_record<types, fields>;

using out_t = seqan3::sequence_file_output<fields>;
using out_ptr = std::unique_ptr<out_t>;
struct FastqPair
{
    out_ptr r1;
    out_ptr r2;
};

// convert fastq files to a list of fastq files
template <class T>
void to_fastq_list(fastq_metadata metadata, T& fin1, T& fin2, std::filesystem::path out,
            std::filesystem::path hash, std::filesystem::path hash_no_quality, 
            const std::size_t& min_length, const std::string& suffix1, const std::string& suffix2) {
    if (!metadata.enough_info()) {
        throw std::runtime_error("The metadata does not contain enough information.");
    }
    const bool check = false;
    std::array<std::string, 2> ID1, ID2;
    std::string rg_id;
    std::size_t n_processed = 0u;
    std::size_t n_skipped = 0u;
    uint64_t sum = 0;
    uint64_t sum_no_quality = 0;
    int nout = 0; // n output files
    std::unordered_map<std::string, FastqPair> fastq_outputs; // for storing output files
    // create the output directory
    std::filesystem::create_directories(out);
    std::filesystem::path out_flie_list = out / "file_list.txt";
    std::ofstream ofs(out_flie_list);

    // iterate over records
    for (auto && [rec1, rec2] : seqan3::views::zip(fin1, fin2)) {// && is important!
        ++n_processed;
        // parse id
        ID1 = parse_ID(rec1.id(), metadata.id_index, suffix1, metadata.illumina_second_id_style);
        ID2 = parse_ID(rec2.id(), metadata.id_index, suffix2, metadata.illumina_second_id_style);
        if (ID1[0] != ID2[0] || ID1[1] != ID2[1]) {
            throw std::runtime_error("Your pairs don't match. ID in the file 1, " + rec1.id() + "; ID in the file 2, " + rec2.id());
        }
        
        rg_id = get_rg_id(ID1, metadata.n_ID_fields, check);
        if (rec1.sequence().size() < min_length || rec2.sequence().size() < min_length) {
            ++n_skipped;
        } else {
            // create output files if necessary
            if (!fastq_outputs.contains(rg_id)) {
                nout++;
                std::string file_idx = std::format("{:03d}", nout);
                ofs << file_idx << ".R1.fastq.gz\t" 
                    << file_idx << ".R2.fastq.gz\t"
                    << rg_id << '\n';
                fastq_outputs.emplace(
                    rg_id,
                    FastqPair{
                        std::make_unique<out_t>(out / (file_idx + ".R1.fastq.gz")),
                        std::make_unique<out_t>(out / (file_idx + ".R2.fastq.gz"))
                    }
                );
            }
            // store records
            auto bq1 = convert_phred_score(rec1.base_qualities(), metadata.format);
            auto bq2 = convert_phred_score(rec2.base_qualities(), metadata.format);
            fq_record_type rec1o{ID1[0] + " " + ID1[1], rec1.sequence(), bq1};
            fq_record_type rec2o{ID2[0] + " " + ID2[1], rec2.sequence(), bq2};
            auto & fout = fastq_outputs.at(rg_id);
            fout.r1->push_back(rec1o);
            fout.r2->push_back(rec2o);
    
            // hash
            if (!hash.empty()) {
                hexSum(calc_hash_str(ID1[0] + "/1" + vec2string(rec1.sequence()) + vec2string(bq1)), sum);
                hexSum(calc_hash_str(ID2[0] + "/2" + vec2string(rec2.sequence()) + vec2string(bq2)), sum);
            }
            if (!hash_no_quality.empty()) {
                hexSum(calc_hash_str((ID1[0] + "/1" + vec2string(rec1.sequence()))), sum_no_quality);
                hexSum(calc_hash_str((ID2[0] + "/2" + vec2string(rec2.sequence()))), sum_no_quality);
            }
        }
    }
    ofs.close();
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
    std::cerr << "Done. Writing out " << n_processed << " record pairs into " << nout << " file pairs and skipped " << n_skipped << " pairs\n";
};
// end of to_fastq_list()

// default values of the arguments
struct cmd_arguments
{
    std::filesystem::path fastq1{};
    std::filesystem::path fastq2{};
    std::filesystem::path out{};
    std::filesystem::path hash{};
    std::filesystem::path hash_no_quality{};
    std::size_t min_length{0u};
    int id_index{-1};
    std::string suffix1{"/1"s};
    std::string suffix2{"/2"s};
    std::string phred{"auto"s};
};

void initialise_parser(sharg::parser & parser, cmd_arguments & args)
{
    parser.info.author = "Shinichi Namba";
    parser.info.short_description = "Splitting paired-end illumina/DNBSEQ fastq files by read groups.";
    parser.info.version = FASTQ2SAM_VERSION_STRING;
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
                                    .description = "The output directory for storing the fastq files per RG. If not set, this program will only scan the first fastq file.", 
                                    .validator = sharg::output_directory_validator{}});
    parser.add_option(args.hash,
                      sharg::config{.long_id = "hash",
                                    .description = "The output file to write the read hash, which is compatible to BamHash."});
    parser.add_option(args.hash_no_quality,
                      sharg::config{.long_id = "hash-no-quality",
                                    .description = "The output file to write the read hash, calculated without read quality."});
    parser.add_option(args.min_length,
                      sharg::config{.short_id = 'm',
                                    .long_id = "min-length",
                                    .description = "The minimum read length to be written to the output file.", 
                                    .advanced = true});
    parser.add_option(args.id_index,
                      sharg::config{.short_id = 'i',
                                    .long_id = "id-index",
                                    .description = "The column index of ID to be parsed (0-based). By default, the index will be automatically estimated.", 
                                    .advanced = true,
                                    .validator = sharg::arithmetic_range_validator{-1, 10}});
    parser.add_option(args.suffix1,
                      sharg::config{.long_id = "suffix1",
                                    .description = "The suffix of the read ID in the first fastq.", 
                                    .advanced = true});
    parser.add_option(args.suffix2,
                      sharg::config{.long_id = "suffix2",
                                    .description = "The suffix of the read ID in the second fastq.", 
                                    .advanced = true});
    parser.add_option(args.phred,
                      sharg::config{.long_id = "phred",
                                    .description = "The phred format.", 
                                    .advanced = true,
                                    .validator = sharg::value_list_validator{"auto", "phred33", "phred64", "solexa64"}});
}


// split_fastq
void split_fastq(const cmd_arguments& args, const fastq_metadata metadata) {
    // iterate over records
    seqan3::debug_stream << "Writing the paired fastq files\n";
    if (metadata.format == phred{0B0100} || metadata.format == phred{0B1000}) {
        sequence_file_input_phred94 fin1{args.fastq1};
        sequence_file_input_phred94 fin2{args.fastq2};
        to_fastq_list(metadata, fin1, fin2, args.out, args.hash, args.hash_no_quality, args.min_length, args.suffix1, args.suffix2);
    } else {
        sequence_file_input_phred68solexa fin1{args.fastq1};
        sequence_file_input_phred68solexa fin2{args.fastq2};
        to_fastq_list(metadata, fin1, fin2, args.out, args.hash, args.hash_no_quality, args.min_length, args.suffix1, args.suffix2);
    }
}

int main(int argc, char ** argv)
{
    // parse args
    sharg::parser myparser{"split_fastq", argc, argv, sharg::update_notifications::off}; // initialise myparser
    cmd_arguments args{};
    initialise_parser(myparser, args); 
    try {
        myparser.parse(); // trigger command line parsing
    } catch (sharg::parser_error const & ext) { // catch user errors
        std::cerr << "[ERROR] " << ext.what() << "\n"; // error message
        return EXIT_FAILURE;
    }
    /*
    if (args.sample_name.empty() && !args.out.empty()) {
        std::runtime_error("`--sample-name` must be provided if `--out` is not empty."s);
    }
    if (args.library.empty()) {
        args.library = args.sample_name;
    }
    */

    // scan
    fastq_metadata metadata = scan_fastq(args.fastq1, args.fastq2, args.id_index, args.suffix1, args.suffix2, string_to_phred(args.phred), !args.out.empty());
    metadata.print();
    if (args.out.empty()) {
        return EXIT_SUCCESS; // scan only
    }

    // split_fastq
    split_fastq(args, metadata);

    return EXIT_SUCCESS;
}
