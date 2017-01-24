#include <iostream>
#include <ios>
#include <fstream>
#include <string>
#include <map>
using std::endl;
using std::cout;
using std::cin;
using std::ifstream;
using std::string;
using std::map;

#define K 5

/* valgrind, sparsehash */

int main(int argc, char **argv) {

    ifstream infile;
    infile.open("./controle.fa");

    string data;
    map<string, int> encounters; //TODO optimize, make it char[5]
    while(getline(infile, data)) {
        if(data.begin() == data.end() || *data.begin() == '>')
        /* filters out descriptive lines and blank lines*/
            continue;
        //data.pop_back(); //TODO

        //TODO optimize? with a sliding array
        for (unsigned int i=0; i <= data.size() - K; i++) {
            string km = data.substr(i, K);
            /* map[x] initialises the value corresponding to x to 0 if non existent */
            ++encounters[km]; //pre-incrementing avoids unnecessary copies
        }
    }
 
    infile.close();

    for (auto const& iterator : encounters) {
        cout << iterator.first << '\t' << iterator.second << endl;
    }

}
