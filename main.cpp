#include "lib/sparsepp.h"
using spp::sparse_hash_map;

#include <utility>
using std::pair;

#include "optionsParser.hpp"
#include "hash.hpp"
#include "degenerate.hpp"
#include "node.hpp"
#include "mutual_information.hpp"
#include "fisher_test.hpp"
#include "stat_container.hpp"


int main (int argc, char **argv) {

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

    const string &km_string = input.getCmdOption("-k");
    if (km_string.empty()){
        cerr << "Please specify a k-mer size putting \"-k size\" to the end of your command line" << endl;
        exit(EXIT_FAILURE);
    }

    int signed_k;
    try {
      signed_k = stoi(km_string);
    } catch (std::invalid_argument& e) {
        cerr << "The k-mer size you input (" << km_string << ") could not be converted to an int." << endl;
        exit(EXIT_FAILURE);
    } catch (std::out_of_range& e) {
        cerr << "The k-mer size you input (" << km_string << ") is too big !" << endl;
        cerr << "Note that most of the k-mer will be unique when the size is >= 15" << endl;
        exit(EXIT_FAILURE);
    }
    if (signed_k < 2) {
        cerr << "The motif size must be greater than 1 !" << endl;
        exit(EXIT_FAILURE);
    }

    unsigned int k = signed_k;
    unsigned int d = k;

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
            if (signed_d < 0 || signed_d > signed_k) {
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

    //vector of pointers to hash_map
    vector<sparse_hash_map<string, pair<int, Node *>> *> hash_map_holder(d+1);

    //0_degree_motifs_hash_map
    hash_map_holder[0] = new sparse_hash_map<string, pair<int, Node *>>();

    unsigned long long int global_motif_count_positive;
    unsigned long long int global_motif_count_negative;

    //counting k-mers
    std::clog << "======== Counting ========" << endl << endl;

    if(input.cmdOptionExists("-p")) {

        std::clog << "\tBeginning to read positive file..." << endl;
        global_motif_count_positive = fill_hash_map_from_pos_positive(*hash_map_holder[0], positive_filename, k, p);
        std::clog << "\tPositive file read." << endl;
        std::clog << "\tBeginning to read negative file..." << endl;
        global_motif_count_negative = fill_hash_map_from_pos_negative(*hash_map_holder[0], negative_filename, k, p);
        std::clog << "\tNegative file read." << endl;
    }
    else {
        std::clog << "\tBeginning to read positive file..." << endl;
        global_motif_count_positive = fill_hash_map_positive(*hash_map_holder[0], positive_filename, k);
        std::clog << "\tPositive file read." << endl;
        std::clog << "\tBeginning to read negative file..." << endl;
        global_motif_count_negative = fill_hash_map_negative(*hash_map_holder[0], negative_filename, k);
        std::clog << "\tNegative file read." << endl;
    }

    //degeneration
    std::clog << endl << "======== Degeneration ========" << endl << endl;

    for (unsigned int i=0; i < d; i++) {

        std::clog << "\tLevel " << i << " : degeneration ongoing" << endl;

        hash_map_holder[i+1] = new sparse_hash_map<string, pair<int, Node *>>();
        degenerate(*(hash_map_holder[i]), *(hash_map_holder[i+1]), k);
    }

    std::clog << endl << "======== Generating MI and p-value ========" << endl << endl;

    //construction of a vector containing <motif, flag, MI, pvalue>
    vector<StatContainer> stat_values;

    int generated_containers = 0;
    for (auto const &hash_map_reference : hash_map_holder) {
        for (auto const &hash_map_entry_ref : *hash_map_reference ) {
            // std::cout << hash_map_entry_ref.first << "\t" << hash_map_entry_ref.second.second->get_positive_count() << "\t" << hash_map_entry_ref.second.second->get_negative_count() << "\t" << global_motif_count_positive << "\t" << global_motif_count_negative << std::endl;
            stat_values.emplace_back(StatContainer( hash_map_entry_ref.first,
                                                    mutual_information (hash_map_entry_ref.second.second->get_positive_count(),
                                                                        hash_map_entry_ref.second.second->get_negative_count(),
                                                                        global_motif_count_positive,
                                                                        global_motif_count_negative),
                                                    fisher_test_p_value(hash_map_entry_ref.second.second->get_positive_count(),
                                                                        hash_map_entry_ref.second.second->get_negative_count(),
                                                                        global_motif_count_positive,
                                                                        global_motif_count_negative,
                                                                        one_tailed)
                                                  )
                                    );
            generated_containers++;
        }
    }

    std::clog << "\tDone generating " << generated_containers << " containers and their content" << endl;
    std::clog << "\tBytes allocated : " << generated_containers << " x " << sizeof(StatContainer) << endl;

    std::clog << endl << "======== Simplification ========" << endl << endl;

    std::clog << "\tSorting the containers by MI" << endl;
    std::sort(stat_values.begin(), stat_values.end(), StatContainer::compare_by_mi);
    for (auto &stat_container : stat_values)
    {
        std::cout << stat_container.get_motif() << "\t" << stat_container.get_mi() << std::endl;
    }

    std::clog << endl << "======== Cleaning ========" << endl << endl;

    for (auto const &hash_map_ptr_reference : hash_map_holder) {
        for (auto const &hash_map_value : *hash_map_ptr_reference ) {
            delete(hash_map_value.second.second);
        }
        delete(hash_map_ptr_reference);
    }
}
