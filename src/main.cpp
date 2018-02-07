#include "sparsepp.h"
using spp::sparse_hash_map;

#include "optionsParser.hpp"
#include "hash.hpp"
#include "degenerate.hpp"
#include "node.hpp"
#include "mutual_information.hpp"
#include "fisher_test.hpp"
#include "graph_simplification.hpp"
#include "meme_format.hpp"
#include "find_redundant_motif.hpp"

#include <utility>
using std::pair;

#include <chrono>

#include <string>
using std::stoi;
using std::stod;


int main (int argc, char **argv) {

    auto start_chrono_parsing_options = std::chrono::high_resolution_clock::now();

    std::initializer_list<const string> positive_file_options = {"-pf", "--positive-file"};
    std::initializer_list<const string> negative_file_options = {"-nf", "--negative-file"};
    std::initializer_list<const string> motif_length_options = {"-l", "--motif-length"};
    std::initializer_list<const string> degeneration_level_options = {"-d", "--degeneration-level"};
    std::initializer_list<const string> output_file_options = {"-o", "--output-file"};
    std::initializer_list<const string> position_options = {"-p", "--position"};
    std::initializer_list<const string> no_reverse_complement_options = {"--norc"};
    std::initializer_list<const string> pvalue_threshold_options = {"-t", "--threshold"};
    std::initializer_list<const string> help_options = {"-h", "--help"};
    std::initializer_list<const string> no_log_options = {"--no-log"};

    InputParser input(argc, argv);

    if (input.cmdOptionExists(no_log_options)) {
        std::clog.setstate(std::ios_base::badbit);
    }

    /* Parsing filename */

    if (input.cmdOptionExists(help_options)) {
        std::cout << "Usage :" << endl <<
                     "dinamo (-pf|--positive-file) path/to/positive (-nf|--negative-file) path/to/negative (-l|--motif-length) k" << endl;
        std::cout << "Available options :" << endl;
        std::cout << "\t(-d|--degeneration-level) k         : Limits the degeneration to at most k positions" << endl;
        std::cout << "\t(-o|--output-file) path/to/output   : Output the meme file to the desired path (has no effect when -p option is used)" << endl;
        std::cout << "\t(-p|--position) k                   : Only process motif that end at position k in the sequence." << endl <<
                     "\t\t(Important note : position 0 corresponds to the last motif of each sequence)" << endl;
        std::cout << "\t--norc                              : When -p is not used, prevents dinamo to link motif to their reverse complement" << endl <<
                     "\t\t(Please be warned : not linking the motif to their reverse complement usually doubles memory usage)" << endl;
        std::cout << "\t(-t|--threshold) r                  : Sets the pvalue threshold to r (0 <= r <= 1)" << endl;
        std::cout << "\t(-h|--help)                         : Displays this help" << endl;
        std::cout << "\t--no-log                            : Prevents the log ouput from being displayed" << endl;
        exit(EXIT_SUCCESS);
    }

    const string &positive_filename = input.getCmdOption(positive_file_options);
    if (positive_filename.empty()){
        cerr << "You did not input the positive sequence file, try adding \"-pf path/to/positive\" to the end of your command line" << endl;
        exit(EXIT_FAILURE);
    }

    const string &negative_filename = input.getCmdOption(negative_file_options);
    if (negative_filename.empty()){
        cerr << "You did not input the negative sequence file, try adding \"-nf path/to/negative\" to the end of your command line" << endl;
        exit(EXIT_FAILURE);
    }

    /* Parsing kmere's size */

    const string &length_string = input.getCmdOption(motif_length_options);
    if (length_string.empty()){
        cerr << "Please input a motif length putting \"-l length\" to the end of your command line" << endl;
        exit(EXIT_FAILURE);
    }

    int signed_l;
    try {
      signed_l = stoi(length_string);
    } catch (std::invalid_argument& e) {
        cerr << "The motif length you input (" << length_string << ") could not be converted to an int." << endl;
        exit(EXIT_FAILURE);
    } catch (std::out_of_range& e) {
        cerr << "The motif length you input (" << length_string << ") is too big !" << endl;
        cerr << "Note that most of the k-mer will be unique when the size is >= 15" << endl;
        exit(EXIT_FAILURE);
    }
    if (signed_l < 2) {
        cerr << "The motif length must be greater than 1 !" << endl;
        exit(EXIT_FAILURE);
    }

    unsigned int l = signed_l;
    unsigned int d = l;

    if(input.cmdOptionExists(degeneration_level_options)) {

        const string &degeneration_level = input.getCmdOption(degeneration_level_options);
        int signed_d;
        if (degeneration_level.empty()) {
            cerr << "You specified the option -d (limit degeneration) but did not input a value." << endl;
            cerr << "Option usage : append \"-d X\" to the end of your command line." << endl;
            cerr << "Description : Limits the degeneration to at most X positions per motif. Please note that X cannot be greater than l." << endl;
            exit(EXIT_FAILURE);
        }
        else {
            try {
                signed_d = stoi(degeneration_level);
            } catch (std::invalid_argument &e) {
                cerr << "The degeneration limit you entered (" << degeneration_level << ") could not be converted to an int." << endl;
                exit(EXIT_FAILURE);
            } catch (std::out_of_range &e) {
                cerr << "The degeneration limit you entered (" << degeneration_level << ") is way too big !" << endl;
                cerr << "Please note that at most l positions can be degenerated (l being the motif size you chose)." << endl;
                exit(EXIT_FAILURE);
            }
            if (signed_d < 0 || signed_d > signed_l) {
                cerr << "The degeneration limit must be positive, and <= l ! (l being the motif size you chose)" << endl;
                exit(EXIT_FAILURE);
            }
            d = signed_d;
        }
    }

    unsigned int p = 0;
    sparse_hash_map<char, unsigned int> neg_nuc_count{
        {'A', 0},
        {'C', 0},
        {'G', 0},
        {'T', 0}
    };

    string meme_filename = input.getCmdOption(output_file_options);

    if(input.cmdOptionExists(position_options)) {
        const string &position = input.getCmdOption(position_options);
        int signed_p;
        if (position.empty()) {
            cerr << "You specified the option -p (count only motif that end at a certain pos) but did not input a value." << endl;
            cerr << "Option usage : append \"-p X\" to the end of your command line." << endl;
            cerr << "Description : count only motif that end at pos X in the sequence." << endl;
            cerr << "Please note that the counting is made from the end (ie. -p 0 will only count the last motif of each sequence)." << endl;
            exit(EXIT_FAILURE);
        }
        else {
            try {
                signed_p = stoi(position);
            } catch (std::invalid_argument &e) {
                cerr << "The position you entered (" << position << ") could not be converted to an int." << endl;
                exit(EXIT_FAILURE);
            } catch (std::out_of_range &e) {
                cerr << "The position you entered (" << position << ") is too big !" << endl;
                cerr << "If you believe you should be able to count a motif this far please consider creating an issue." << endl;
                exit(EXIT_FAILURE);
            }
            if (signed_p < 0) {
                cerr << "The degeneration limit must be positive !" << endl;
                exit(EXIT_FAILURE);
            }
            p = signed_p;
        }
    } else {
        if (meme_filename.empty()){
            std::clog << "warning : you did not input a meme output file. \"dinamo_results.meme\" will be created/overriden as default output." << endl;
            meme_filename = "dinamo_results.meme";
        }
    }

    bool rc = !input.cmdOptionExists(position_options);

    if(input.cmdOptionExists(position_options) && input.cmdOptionExists(no_reverse_complement_options)) {
        const string &norc_str = input.getCmdOption(no_reverse_complement_options);
        if (!norc_str.empty()) {
            std::cerr << "warning : you provided an argument (" << norc_str << ") to the -norc option but this option doesn't require one." << endl;
        }
        rc = false;
    }

    double alpha = 0.05;

    if(input.cmdOptionExists(pvalue_threshold_options)) {
        const string &threshold_str = input.getCmdOption(pvalue_threshold_options);
        if (threshold_str.empty()) {
            cerr << "You specified the option -t (set pvalue threshold) but did not input a value." << endl;
            cerr << "Option usage : append \"-t X\" to the end of your command line." << endl;
            cerr << "Description : sets the threshold of the Holm-Bonferroni to X." << endl;
            cerr << "Please note that X must be a value between 0 and 1." << endl;
            exit(EXIT_FAILURE);
        } else {
            try {
                alpha = stod(threshold_str);
            } catch (std::invalid_argument &e) {
                cerr << "The alpha you entered (" << threshold_str << ") could not be converted to a double." << endl;
                exit(EXIT_FAILURE);
            } catch (std::out_of_range &e) {
                cerr << "The alpha you entered (" << threshold_str << ") is too big !" << endl;
                cerr << "Please remember that 0 <= alpha <= 1" << endl;
                exit(EXIT_FAILURE);
            }
            if (alpha < 0) {
                cerr << "Alpha must be positive !" << endl;
                exit(EXIT_FAILURE);
            } else if (alpha > 1){
                cerr << "Alpha must be less than 1 !" << endl;
                exit(EXIT_FAILURE);
            }
        }
    }

    auto end_chrono_parsing_options = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> parsing_options_time = end_chrono_parsing_options - start_chrono_parsing_options;
    std::clog << "Parsed options in : " << parsing_options_time.count() << " seconds\n";

    //vector of pointers to hash_map
    vector<sparse_hash_map<string, Node *> *> hash_map_holder(d+1);

    //0_degree_motifs_hash_map
    hash_map_holder[0] = new sparse_hash_map<string, Node *>();

    unsigned int global_motif_count_positive;
    unsigned int global_motif_count_negative;

    //counting k-mers

    auto start_chrono_counting = std::chrono::high_resolution_clock::now();
    std::clog << "======== Counting ========" << endl << endl;

    std::clog << "\tBeginning to read positive file..." << endl;
    if(input.cmdOptionExists(position_options))
         global_motif_count_positive = fill_hash_map_from_pos(*hash_map_holder[0], positive_filename, l, p, true, rc);
    else global_motif_count_positive = fill_hash_map         (*hash_map_holder[0], positive_filename, l,    true, rc, neg_nuc_count);
    std::clog << "\tPositive file read." << endl;

    std::clog << "\tBeginning to read negative file..." << endl;
    if(input.cmdOptionExists(position_options))
         global_motif_count_negative = fill_hash_map_from_pos(*hash_map_holder[0], negative_filename, l, p, false, rc);
    else global_motif_count_negative = fill_hash_map         (*hash_map_holder[0], negative_filename, l,    false, rc, neg_nuc_count);
    std::clog << "\tNegative file read." << endl;

    auto end_chrono_counting = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> counting_time = end_chrono_counting - start_chrono_counting;
    std::clog << "Counting done in : " << counting_time.count() << " seconds" << endl;

    auto start_chrono_degeneration = std::chrono::high_resolution_clock::now();

    std::clog << endl << "======== Degeneration ========" << endl << endl;

    std::clog << "Before : Total motif count before degeneration" << endl;
    std::clog << "After  : Total motif count after degeneration" << endl;
    std::clog << "\t\t\tBefore\t\tAfter\t\tRatio" << endl;
    unsigned int total_motif_created = hash_map_holder[0]->size();

    for (unsigned int i=0; i < d; i++) {

        double count_before_deg = total_motif_created;
        std::clog << "\tLevel " << i+1 << "\t\t" << total_motif_created << "\t\t";

        hash_map_holder[i+1] = new sparse_hash_map<string, Node *>();
        degenerate(*(hash_map_holder[i]), *(hash_map_holder[i+1]), l, rc);

        total_motif_created += hash_map_holder[i+1]->size();
        std::clog << total_motif_created << "\t\t" << total_motif_created / count_before_deg << endl;
    }

    auto end_chrono_degeneration = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> degeneration_time = end_chrono_degeneration - start_chrono_degeneration;
    std::clog << "Degeneration done in : " << degeneration_time.count() << " seconds\n";

    // for (auto const &hash_map_reference : hash_map_holder)
    //     for (auto const &motif : *hash_map_reference)
    //         std::cerr << motif.first << "\t" << motif.second->get_positive_count() << "\t" << motif.second->get_negative_count() << "\t" << motif.second->get_predecessors().size() << "\t" << motif.second->get_successors().size() << endl;


    auto start_chrono_simplification = std::chrono::high_resolution_clock::now();
    std::clog << endl << "======== Simplification ========" << endl << endl;

    std::clog << "\tGenerating MIs..." << endl;
    //TODO get the number of nodes to define the vector as an array
    std::vector<pair<const string, Node *> *> mi_sorted_hash_map_entries;
    std::vector<pair<const string, Node *> *> pvalue_sorted_hash_map_entries;

    for (auto const &hash_map_reference : hash_map_holder) {
        for (auto &hash_map_entry_ref : *hash_map_reference ) {
            mi_sorted_hash_map_entries.push_back(&hash_map_entry_ref);
            hash_map_entry_ref.second->calculate_mi(global_motif_count_positive, global_motif_count_negative);
        }
    }

    std::clog << "\tDone generating MIs." << endl << endl;

    std::clog << "\tSorting entries by MI..." << endl << endl;
    std::sort(  mi_sorted_hash_map_entries.begin(),
                mi_sorted_hash_map_entries.end(),
                [] (const auto entry_one, const auto entry_two) {
                    return entry_one->second->get_mi() > entry_two->second->get_mi();
                }
             );


    std::clog << "\tSimplificating graph..." << endl;

    // graph_simplification(mi_sorted_hash_map_entries, input.cmdOptionExists(position_options));
    graph_simplification(mi_sorted_hash_map_entries, false);

    std::clog << "\tSimplification done !" << endl << endl;

    std::clog << "\tGenerating pvalue of the remaining entries..." << endl;
    unsigned int m = 0;
    for (auto &entry : mi_sorted_hash_map_entries) {
        if (entry->second->get_state() == validated) {
            entry->second->suppress();
            entry->second->calculate_pvalue(global_motif_count_positive, global_motif_count_negative);
            pvalue_sorted_hash_map_entries.push_back(entry);
            m++;
        }
    }
    std::clog << "\tDone generating pvalues." << endl << endl;

    std::clog << "\tSorting the remaining entries by pvalue..." << endl;
    std::sort(  pvalue_sorted_hash_map_entries.begin(),
                pvalue_sorted_hash_map_entries.end(),
                [] (const auto entry_one, const auto entry_two) {
                    return entry_one->second->get_pvalue() < entry_two->second->get_pvalue();
                }
             );

    auto end_chrono_simplification = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> simplification_time = end_chrono_simplification - start_chrono_simplification;
    std::clog << "Simplification made in : " << simplification_time.count() << " seconds\n";


    auto start_chrono_holm_test = std::chrono::high_resolution_clock::now();

    std::clog << endl << "======== Holm-Bonferroni ========" << endl << endl;

    std::clog << "\tApplying Holm-Bonferroni method to determine independance of observed values in both files" << endl;

    mi_sorted_hash_map_entries.clear();

    unsigned int k = 0;
    while (pvalue_sorted_hash_map_entries[k]->second->get_pvalue()
           <= (alpha / ((double)(m + 1 - k))) //alpha is defined when parsing options
          ) {
        mi_sorted_hash_map_entries.push_back(pvalue_sorted_hash_map_entries[k]);
        k++;
    }

    std::sort(  mi_sorted_hash_map_entries.begin(),
                mi_sorted_hash_map_entries.end(),
                [] (const auto entry_one, const auto entry_two) {
                    return entry_one->second->get_mi() > entry_two->second->get_mi();
                }
             );

    auto end_chrono_holm_test = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> holm_test_time = end_chrono_holm_test - start_chrono_holm_test;
    std::clog << "Holm method performed in : " << holm_test_time.count() << " seconds\n";

    std::clog << endl << "======== Filtering redundant motif ========" << endl << endl;

    if (!input.cmdOptionExists(position_options))
        filter_redundant_motif(mi_sorted_hash_map_entries, l);

    std::clog << endl << "======== Results ========" << endl << endl;
    std::cout << "# Motif"
              << "\t" << "MI"
              << "\t" << "PValue"
              << "\t" << "#MotifPos/#TotalPos"
              << "\t" << "#MotifNeg/#TotalNeg"
              << endl;    
    for (auto &entry : mi_sorted_hash_map_entries) {
      std::cout << entry->first
                << "\t" << entry->second->get_mi()
                << "\t" << entry->second->get_pvalue()
                << "\t" << entry->second->get_positive_count() << "/" << global_motif_count_positive
                << "\t" << entry->second->get_negative_count() << "/" << global_motif_count_negative
                << endl;
    }

    if (!input.cmdOptionExists(position_options))
        create_meme_file(mi_sorted_hash_map_entries, *hash_map_holder[0], neg_nuc_count, l, meme_filename);

    // std::clog << endl << "======== Cleaning ========" << endl << endl;
    //
    // for (auto const &hash_map_ptr_reference : hash_map_holder) {
    //     for (auto const &hash_map_value : *hash_map_ptr_reference ) {
    //         delete(hash_map_value.second.second);
    //     }
    //     delete(hash_map_ptr_reference);
    // }
}
