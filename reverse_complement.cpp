#include "reverse_complement.hpp"

sparse_hash_map<char, char> iupac_to_complement {
    {'A','T'},
    {'C','G'},
    {'G','C'},
    {'T','A'},
    {'M','K'},
    {'W','W'},
    {'R','Y'},
    {'K','M'},
    {'S','S'},
    {'Y','R'},
    {'V','B'},
    {'H','D'},
    {'D','H'},
    {'B','V'},
    {'N','N'}
};

string reverse_complement(string &motif) {
    string rv(motif);
    transform(rv.begin(), rv.end(), rv.begin(),
    [](char nuc) {return iupac_to_complement[nuc]; });
    reverse(rv.begin(), rv.end());
    return rv;
}
