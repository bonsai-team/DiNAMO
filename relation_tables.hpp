#if !defined(__RELATION_TABLES_HPP_INCLUDED__)

    #define  __RELATION_TABLES_HPP_INCLUDED__

    //=====================
    //Included dependencies
    //=====================

    #include <unordered_set>
    using std::unordered_set;

    #include <vector>
    using std::vector;

    #include <string>
    using std::string;

    //=====================
    //  Included libs
    //=====================

    #include "lib/sparsepp.h"
    using spp::sparse_hash_map;

    //=====================

    extern const vector<char> nucleotides;

    extern sparse_hash_map<string, vector<unordered_set<char>>> nucs_to_iupacs;

    extern sparse_hash_map<char, unordered_set<char>> iupacs_dependencies;

    //second member needs to be ordered (ordering used by set_intersection)
    extern sparse_hash_map<char, vector<char>> iupac_to_nucs;

    extern sparse_hash_map<char, vector<char>> iupac_to_successors;

#endif
