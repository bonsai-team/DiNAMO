#if defined(__USE_STD_UNORDERED_MAP__)
    /* Use Standard library */
    #include <unordered_map>
    using std::unordered_map;
    #define __MAP__ unordered_map<string, int>

#else
    /* Use Sparse++ */
    #include "sparsepp.h"
    using spp::sparse_hash_map;
    #define __MAP__ sparse_hash_map<string, int>

#endif

#include <iostream>
using std::endl;
using std::cerr;
using std::cout;

#include <fstream>
using std::ifstream;

#include <string>
using std::string;

#include <vector>
using std::vector;

#include <deque>
using std::deque;

#include <cstdlib>
using std::atoi;

#include "optionsParser.hpp"

int main(int argc, char **argv) {

    InputParser input(argc, argv);

    /* Parsing filename */

    const string &filepath = input.getCmdOption("-f");
    if (filepath.empty()){
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

    ifstream infile;
    infile.open(filepath);

    // checking if the file could be opened
    if (!infile.is_open()) {
        cerr << "Error the file could not be opened, make sure that \"" << filepath << "\" exists and isn't a directory" << endl;
        exit(EXIT_FAILURE);
    }

    string data;
    __MAP__ encounters;  //__MAP__ is defined at compilation time, check Makefile
    deque<char> deque(k);
    bool dequeReady = false;

    while(getline(infile, data)) {
        if(data.begin() == data.end() || *data.begin() == '>') {
            dequeReady = false;
            deque.clear();
        }
        else {

            /* -- Important note --
            The next for-loop is invalid for k=1, because the deque is set to ready at the end of the loop.
            But having k=1 is not supposed to happen since it means user wants to count single nucleotides,
            basically asking the count of A, T, C, G. It is not what the program is meant for.
            */

            for (unsigned int i = 0; i < data.size(); i++) {
                //filtering accepted characters in fasta file
                switch(data[i]) {
                    case 'A':
                    case 'C':
                    case 'G':
                    case 'T':
                        deque.emplace_back(data[i]);
                        break;
                    //lowercases to be added as uppercases
                    case 'a':
                    case 't':
                    case 'g':
                    case 'c':
                        deque.emplace_back(toupper(data[i]));
                        break;
                    //if character is invalid, emptying deque and refilling it
                    default:
                        deque.clear();
                        dequeReady = false;
                        break;
                }
                if (dequeReady) {
                    string kmer(deque.begin(), deque.end());
                    ++encounters[kmer];
                    deque.pop_front();
                }
                else
                    if (deque.size() == (k-1))
                        dequeReady = true;
            }
        }
    }

    infile.close();
    return 0;
}
