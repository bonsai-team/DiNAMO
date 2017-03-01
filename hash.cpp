#include <iostream>
#include <ios>
#include <fstream>
#include <string>
#include <vector>

using std::endl;
using std::cout;
using std::cin;
using std::ifstream;
using std::string;
using std::atoi;
using std::vector;

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

    int k = stoi(kmString);

    ifstream infile;
    infile.open(filename);

    string data; 
    _MAP encounters;  //_MAP is defined at compilation time, check Makefile

    while(getline(infile, data)) {
        if(data.begin() == data.end() || *data.begin() == '>')
            continue;
        //data.pop_back(); //TODO second type of km parsing

        for (unsigned long int i=0; i <= data.size() - k; i++) {
            string km = data.substr(i, k);
            for (auto & c: km) c = toupper(c);
            ++encounters[km]; //pre-incrementing avoids unnecessary copies
        }// maybe switch to readsome 50 for large files
    }

    infile.close();

    for (auto const& iterator : encounters) {
        cout << iterator.first << '\t' << iterator.second << endl;
    }
}

/* valgrind */
