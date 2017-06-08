//=====================
//      Guards
//=====================

#if !defined(__HASH_HPP_INCLUDED__)

    #define  __HASH_HPP_INCLUDED__

    //=====================
    //Included dependencies
    //=====================

    #include <iostream>
    using std::endl;
    using std::cerr;

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

    #include <utility>
    using std::pair;
    using std::make_pair;

    #include <iterator>
    using std::next;

    #include <algorithm>
    using std::transform;
    using std::reverse;

    //=====================
    //  Included libs
    //=====================

    #include "node.hpp"

    #include "sparsepp.h"
    using spp::sparse_hash_map;

    #include "reverse_complement.hpp"

    //=====================

    unsigned int fill_hash_map(sparse_hash_map<string, pair<int, Node *>> &, const string &, unsigned int, bool, bool, sparse_hash_map<char, unsigned int> &);

    unsigned int fill_hash_map_from_pos(sparse_hash_map<string, pair<int, Node *>> &, const string &, unsigned int, unsigned int, bool, bool);

#endif
