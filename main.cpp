#include "lib/sparsepp.h"
using spp::sparse_hash_map;

#include "optionsParser.hpp"
#include "hash.hpp"
#include "degenerate.hpp"
#include "node.hpp"
#include "mutual_information.hpp"
#include "fisher_test.hpp"
#include "graph_simplification.hpp"

#include <utility>
using std::pair;

#include <chrono>


int main (int argc, char **argv) {

    auto start_chrono_parsing_options = std::chrono::high_resolution_clock::now();

    InputParser input(argc, argv);

    /* Parsing filename */

    const string &positive_filename = input.getCmdOption("-pf");
    if (positive_filename.empty()){
        cerr << "You did not input the positive sequence file, try adding \"-pf path/to/positive\" to the end of your command line" << endl;
        exit(EXIT_FAILURE);
    }

    const string &negative_filename = input.getCmdOption("-nf");
    if (negative_filename.empty()){
        cerr << "You did not input the negative sequence file, try adding \"-nf path/to/negative\" to the end of your command line" << endl;
        exit(EXIT_FAILURE);
    }

    /* Parsing kmere's size */

    const string &length_string = input.getCmdOption("-l");
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

    if(input.cmdOptionExists("-d")) {

        const string &degeneration_level = input.getCmdOption("-d");
        int signed_d;
        if (degeneration_level.empty()) {
            cerr << "You specified the option -d (limit degeneration) but did not input a value." << endl;
            cerr << "Option usage : append \"-d X\" to the end of your command line." << endl;
            cerr << "Description : Limits the degeneration to at most X positions per motif. Please note that X cannot be greater than k." << endl;
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
                cerr << "Please note that at most k positions can be degenerated (k being the motif size you chose)." << endl;
                exit(EXIT_FAILURE);
            }
            if (signed_d < 0 || signed_d > signed_l) {
                cerr << "The degeneration limit must be positive, and <= k ! (k being the motif size you chose)" << endl;
                exit(EXIT_FAILURE);
            }
            d = signed_d;
        }
    }

    unsigned int p = 0;

    if(input.cmdOptionExists("-p")) {

        const string &position = input.getCmdOption("-p");
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
                cerr << "If you believe you should be able to count a motif this far please consider opening a ticket." << endl;
                exit(EXIT_FAILURE);
            }
            if (signed_p < 0) {
                cerr << "The degeneration limit must be positive !" << endl;
                exit(EXIT_FAILURE);
            }
            p = signed_p;
        }
    }

    auto end_chrono_parsing_options = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> parsing_options_time = end_chrono_parsing_options - start_chrono_parsing_options;
    std::clog << "Parsed options in : " << parsing_options_time.count() << " seconds\n";

    //vector of pointers to hash_map
    vector<sparse_hash_map<string, pair<int, Node *>> *> hash_map_holder(d+1);

    //0_degree_motifs_hash_map
    hash_map_holder[0] = new sparse_hash_map<string, pair<int, Node *>>();

    unsigned int global_motif_count_positive;
    unsigned int global_motif_count_negative;

    //counting k-mers

    auto start_chrono_counting = std::chrono::high_resolution_clock::now();
    std::clog << "======== Counting ========" << endl << endl;

    if(input.cmdOptionExists("-p")) {

        std::clog << "\tBeginning to read positive file..." << endl;
        global_motif_count_positive = fill_hash_map_from_pos(*hash_map_holder[0], positive_filename, l, p, true);
        std::clog << "\tPositive file read." << endl;
        std::clog << "\tBeginning to read negative file..." << endl;
        global_motif_count_negative = fill_hash_map_from_pos(*hash_map_holder[0], negative_filename, l, p, false);
        std::clog << "\tNegative file read." << endl;
    }
    else {
        std::clog << "\tBeginning to read positive file..." << endl;
        global_motif_count_positive = fill_hash_map(*hash_map_holder[0], positive_filename, l, true);
        std::clog << "\tPositive file read." << endl;
        std::clog << "\tBeginning to read negative file..." << endl;
        global_motif_count_negative = fill_hash_map(*hash_map_holder[0], negative_filename, l, false);
        std::clog << "\tNegative file read." << endl;
    }

    auto end_chrono_counting = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> counting_time = end_chrono_counting - start_chrono_counting;
    std::clog << "Counting done in : " << counting_time.count() << " seconds" << endl;


    auto start_chrono_degeneration = std::chrono::high_resolution_clock::now();

    //degeneration
    std::clog << endl << "======== Degeneration ========" << endl << endl;

    for (unsigned int i=0; i < d; i++) {

        std::clog << "\tLevel " << i << " : degeneration ongoing" << endl;

        hash_map_holder[i+1] = new sparse_hash_map<string, pair<int, Node *>>();
        degenerate(*(hash_map_holder[i]), *(hash_map_holder[i+1]), l);
    }

    auto end_chrono_degeneration = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> degeneration_time = end_chrono_degeneration - start_chrono_degeneration;
    std::clog << "Degeneration done in : " << degeneration_time.count() << " seconds\n";


    auto start_chrono_simplification = std::chrono::high_resolution_clock::now();
    std::clog << endl << "======== Simplification ========" << endl << endl;

    std::clog << "\tGenerating MIs..." << endl;
    //TODO get the number of nodes to define the vector as an array
    std::vector<pair<const string, pair<int, Node *> > *> mi_sorted_hash_map_entries;
    std::vector<pair<const string, pair<int, Node *> > *> pvalue_sorted_hash_map_entries;

    for (auto const &hash_map_reference : hash_map_holder) {
        for (auto &hash_map_entry_ref : *hash_map_reference ) {
            mi_sorted_hash_map_entries.push_back(&hash_map_entry_ref);
            hash_map_entry_ref.second.second->calculate_mi(global_motif_count_positive, global_motif_count_negative);
        }
    }

    std::clog << "\tDone generating MIs." << endl << endl;

    std::clog << "\tSorting entries by MI..." << endl << endl;
    std::sort(  mi_sorted_hash_map_entries.begin(),
                mi_sorted_hash_map_entries.end(),
                [] (const auto entry_one, const auto entry_two) {
                    return entry_one->second.second->get_mi() > entry_two->second.second->get_mi();
                }
             );


    std::clog << "\tSimplificating graph..." << endl;

    graph_simplification(mi_sorted_hash_map_entries, input.cmdOptionExists("-p"));

    std::clog << "\tSimplification done !" << endl << endl;

    std::clog << "\tGenerating pvalue of the remaining entries..." << endl;
    unsigned int m = 0;
    for (auto &entry : mi_sorted_hash_map_entries) {
        if (entry->second.second->get_state() == validated) {
            entry->second.second->calculate_pvalue(global_motif_count_positive, global_motif_count_negative);
            pvalue_sorted_hash_map_entries.push_back(entry);
            m++;
        }
    }
    std::clog << "\tDone generating pvalues." << endl << endl;

    std::clog << "\tSorting the remaining entries by pvalue..." << endl;
    std::sort(  pvalue_sorted_hash_map_entries.begin(),
                pvalue_sorted_hash_map_entries.end(),
                [] (const auto entry_one, const auto entry_two) {
                    return entry_one->second.second->get_pvalue() < entry_two->second.second->get_pvalue();
                }
             );

    auto end_chrono_simplification = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> simplification_time = end_chrono_simplification - start_chrono_simplification;
    std::clog << "Simplification made in : " << simplification_time.count() << " seconds\n";


    auto start_chrono_holm_test = std::chrono::high_resolution_clock::now();

    std::clog << endl << "======== Results ========" << endl << endl;

    std::clog << "\tApplying Holm-Bonferroni method to determine independance of observed values in both files" << endl;
    std::clog << endl << "The appearance of the following values were found to depend on the file you look in" << endl;
    unsigned int k = 0;
    double alpha = 0.05;
    while (pvalue_sorted_hash_map_entries[k]->second.second->get_pvalue()
           <= (alpha / ((double)(m + 1 - k)))
          ) {
        std::cout << pvalue_sorted_hash_map_entries[k]->first << endl;
        k++;
    }

    auto end_chrono_holm_test = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> holm_test_time = end_chrono_holm_test - start_chrono_holm_test;
    std::clog << "Holm method performed in : " << holm_test_time.count() << " seconds\n";

    std::clog << endl << "======== Cleaning ========" << endl << endl;

    for (auto const &hash_map_ptr_reference : hash_map_holder) {
        for (auto const &hash_map_value : *hash_map_ptr_reference ) {
            delete(hash_map_value.second.second);
        }
        delete(hash_map_ptr_reference);
    }
}
