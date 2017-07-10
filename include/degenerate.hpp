#if !defined(__DEGENERATE_HPP_INCLUDED__)

    #define  __DEGENERATE_HPP_INCLUDED__

    //=====================
    //     Dependencies
    //=====================

    #include <iostream>

    #include <string>
    using std::string;

    #include <vector>
    using std::vector;

    #include <unordered_set>
    using std::unordered_set;

    #include <utility>
    using std::pair;
    using std::make_pair;

    //=====================
    //     Libs needed
    //=====================

    #include "node.hpp"

    #include "sparsepp.h"
    using spp::sparse_hash_map;

    #include "reverse_complement.hpp"

    #include "relation_tables.hpp"

    //=====================

    const string find_neighbor_motifs(sparse_hash_map<string, Node *> &, sparse_hash_map<char, pair<string, Node *>> &, const string &, unsigned int);

    void degenerate(sparse_hash_map<string, Node *> &, sparse_hash_map<string, Node *> &, unsigned int, bool);

#endif
