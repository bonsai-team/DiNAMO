#include "lib/sparsepp.h"
using spp::sparse_hash_map;

#include "optionsParser.hpp"
#include "hash.hpp"

int main (int argc, char **argv) {

    InputParser input(argc, argv);

    /* Parsing filename */

    const string &filename = input.getCmdOption("-f");
    if (filename.empty()){
        cerr << "You did not input the file name, try adding \"-f path/to/your/file\" to the end of your command line" << endl;
        exit(EXIT_FAILURE);
    }

    /* Parsing kmere's size */

    const string &kmString = input.getCmdOption("-k");
    if (kmString.empty()){
        cerr << "Please specify a k-mer size putting \"-k size\" to the end of your command line" << endl;
        exit(EXIT_FAILURE);
    }

    int signedK;
    try {
      signedK = stoi(kmString);
    } catch (std::invalid_argument& e) {
        cerr << "The k-mer size you entered (" << kmString << ") could not be converted to an int." << endl;
        exit(EXIT_FAILURE);
    } catch (std::out_of_range& e) {
        cerr << "The k-mer size you entered (" << kmString << ") is too big !" << endl;
        cerr << "Note that most of the k-mer will be unique when the size is >= 15" << endl;
        exit(EXIT_FAILURE);
    }
    if (signedK < 0) {
        cerr << "K-mer size must be greater than or equal to 2 !" << endl;
        exit(EXIT_FAILURE);
    }

    unsigned int k = signedK;
    sparse_hash_map<string, int> hash_map;

    //counting k-mers
    fill_hash_map(hash_map, filename, k);

    for(auto const &it : hash_map) {
        std::cout << it.first << "\t" << it.second << endl;
    }
}
