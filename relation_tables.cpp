#include "relation_tables.hpp"

const vector<char> nucleotides {'A','C','G','T'};

sparse_hash_map<string, vector<unordered_set<char>>> nucs_to_iupacs {
    {"ACGT", {{'Y','S','K','R','W','M'}, {'B','D','H','V'}, {'N'}}},
    {"ACG" , {{'M','R','S'}, {'V'}}},
    {"ACT" , {{'M','W','Y'}, {'H'}}},
    {"AGT" , {{'R','W','K'}, {'D'}}},
    {"CGT" , {{'S','Y','K'}, {'B'}}},
    {"AC"  , {{'M'}}},
    {"AG"  , {{'R'}}},
    {"AT"  , {{'W'}}},
    {"CG"  , {{'S'}}},
    {"CT"  , {{'Y'}}},
    {"GT"  , {{'K'}}}
};

sparse_hash_map<char, unordered_set<char>> iupacs_dependencies {
    {'N', {'B','D','H','V'}},
    {'B', {'Y','S','K'}},
    {'D', {'R','W','K'}},
    {'H', {'Y','W','M'}},
    {'V', {'R','S','M'}},
    {'R', {'A','G'}},
    {'Y', {'C','T'}},
    {'S', {'G','C'}},
    {'W', {'A','T'}},
    {'K', {'G','T'}},
    {'M', {'A','C'}}
};

sparse_hash_map<char, vector<char>> iupac_to_nucs {
    {'N', {'A', 'C', 'G', 'T'}},
    {'B', {'C', 'G', 'T'}},
    {'D', {'A', 'G', 'T'}},
    {'H', {'A', 'C', 'T'}},
    {'V', {'A', 'C', 'G'}},
    {'Y', {'C', 'T'}},
    {'S', {'C', 'G'}},
    {'K', {'G', 'T'}},
    {'R', {'A', 'G'}},
    {'W', {'A', 'T'}},
    {'M', {'A', 'C'}},
    {'A', {'A'}},
    {'C', {'C'}},
    {'G', {'G'}},
    {'T', {'T'}}
};
