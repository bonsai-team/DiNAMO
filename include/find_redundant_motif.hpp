#if !defined(__FIND_REDUNDANT_MOTIF_HPP_INCLUDED__)

    #define  __FIND_REDUNDANT_MOTIF_HPP_INCLUDED__

    //=====================
    //Included dependencies
    //=====================

    #include <vector>
    using std::vector;

    #include <string>
    using std::string;

    #include <utility>
    using std::pair;

    #include <algorithm>
    using std::find_if;
    using std::find;
    using std::remove_if;

    #include <cmath>
    using std::ceil;

    #include <unordered_set>
    using std::unordered_set;

    //=====================
    //  Included libs
    //=====================

    #include "sparsepp.h"
    using spp::sparse_hash_map;

    #include "node.hpp"

    #include "reverse_complement.hpp"

    #include "relation_tables.hpp"

    //=====================

    void filter_redundant_motif(vector<pair<const string, Node *> *> &, unsigned int);

#endif
