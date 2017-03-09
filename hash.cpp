#include <iostream>
#include <ios>
#include <fstream>
#include <string>
#include <vector>
#include <deque>

using std::endl;
using std::cout;
using std::cin;
using std::ifstream;
using std::string;
using std::atoi;
using std::vector;
using std::deque;

/* Standard Library */
#include <map>
using std::map;

/* Sparsepp Library */
#include "lib/sparsepp.h"
using spp::sparse_hash_map;

/* Class to ease parsing 
see http://stackoverflow.com/questions/865668/how-to-parse-command-line-arguments-in-c */
// TODO separate file for parser


class InputParser{
    public:
        InputParser (int &argc, char **argv){
            for (int i=1; i < argc; ++i)
                this->tokens.push_back(string(argv[i]));
        }

        const string& getCmdOption(const string &option) const{
            vector<string>::const_iterator itr;
            itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
            if (itr != this->tokens.end() && ++itr != this->tokens.end()){
                return *itr;
            }
            return empty_string;
        }

        bool cmdOptionExists(const string &option) const{
            return std::find(this->tokens.begin(), this->tokens.end(), option)
                   != this->tokens.end();
        }
    private:
        vector <string> tokens;
        string empty_string;
};

int main(int argc, char **argv) {

    InputParser input(argc, argv);

    /* Parsing filename */

    const string &filename = input.getCmdOption("-f");
    if (filename.empty()){
        exit(EXIT_FAILURE);
    }

    /* Parsing kmere's size */

    const string &kmString = input.getCmdOption("-k");
    if (kmString.empty()){
        exit(EXIT_FAILURE);
    }

    unsigned int k = stoi(kmString);

    ifstream infile;
    infile.open(filename);

    string data; 
    _MAP encounters;  //_MAP is defined at compilation time, check Makefile
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
            basically asking the count of A, T, C, G. It is not the goal of the programm.
            */

            for (unsigned int i = 0; i < data.size(); i++) {
                deque.emplace_back(data[i]); //TODO check whether data[i] is valid
                if (dequeReady) {
                    string kmere(deque.begin(), deque.end());
                    ++encounters[kmere];
                    deque.pop_front();
                }
                else 
                    if (deque.size() == (k-1))
                        dequeReady = true;
            }
        }
    }

    infile.close();

    for (auto const& iterator : encounters) {
        cout << iterator.first << '\t' << iterator.second << endl;
    }
}

