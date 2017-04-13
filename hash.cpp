#include "hash.hpp"


void fill_hash_map(sparse_hash_map<string, pair<int, int>> &encounters, const string &filepath, unsigned int k) {

    ifstream infile;
    infile.open(filepath);

    // checking if the file could be opened
    if (!infile.is_open()) {
        cerr << "Error the file could not be opened, make sure that \"" << filepath << "\" exists and isn't a directory" << endl;
        exit(EXIT_FAILURE);
    }

    string data;
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

            for (char nuc : data) {
                //filtering accepted characters in fasta file
                switch(nuc) {
                    case 'A':
                    case 'C':
                    case 'G':
                    case 'T':
                        deque.emplace_back(nuc);
                        break;
                    //lowercases to be added as uppercases
                    case 'a':
                    case 't':
                    case 'g':
                    case 'c':
                        deque.emplace_back(toupper(nuc));
                        break;
                    //if character is invalid, emptying deque and refilling it
                    default:
                        deque.clear();
                        dequeReady = false;
                        break;
                }
                if (dequeReady) {
                    string motif(deque.begin(), deque.end());
                    auto motif_it = encounters.find(motif);
                    if (motif_it != encounters.end())
                        ++((motif_it->second).first);
                    else
                        encounters.emplace(make_pair(motif, make_pair(1, -1)));
                    deque.pop_front();
                }
                else
                    if (deque.size() == (k-1))
                        dequeReady = true;
            }
        }
    }
    infile.close();
    return;
}

void fill_hash_map_from_pos(sparse_hash_map<string, pair<int, int>> &encounters, const string &filepath, unsigned int k, unsigned int p) {
    ifstream infile;
    infile.open(filepath);

    // checking if the file could be opened
    if (!infile.is_open()) {
        cerr << "Error the file could not be opened, make sure that \"" << filepath << "\" exists and isn't a directory" << endl;
        exit(EXIT_FAILURE);
    }

    string data;
    deque<char> deque(k + p);

    while(getline(infile, data)) {
        if(data.begin() == data.end() || *data.begin() == '>') {
            on_sequence_end(deque, encounters, k, p);
            deque.clear();
        }
        else {
            for (char nuc : data) {
                if(deque.size() == k + p)
                    deque.pop_front();
                switch(nuc) {
                    case 'a':
                    case 'c':
                    case 'g':
                    case 't':
                        deque.emplace_back(toupper(nuc));
                        break;
                    default :
                        deque.emplace_back(nuc);
                        break;
                }
            }
        }
    }
    on_sequence_end(deque, encounters, k, p);
    return;
}

void on_sequence_end(deque<char> &deque, sparse_hash_map<string, pair<int, int>> &encounters, unsigned int k, unsigned int p) {
    if(deque.size() == k + p) {
        string motif(deque.begin(), next(deque.begin(), k));
        if(motif.find_first_not_of("ACGT") != string::npos)
            return;
        auto motif_it = encounters.find(motif);
        if (motif_it != encounters.end())
            ++((motif_it->second).first);
        else
            encounters.emplace(make_pair(motif, make_pair(1, -1)));
    }
    return;
}
