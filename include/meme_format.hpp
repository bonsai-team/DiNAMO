#if !defined(__MEME_FORMAT_HPP_INCLUDED__)

    #define  __MEME_FORMAT_HPP_INCLUDED__

    //=====================
    //Included dependencies
    //=====================

    #include <vector>
    using std::vector;

    #include <fstream>
    using std::ofstream;

    #include <string>
    using std::string;

    #include <utility>
    using std::pair;

    #include <unordered_set>
    using std::unordered_set;

    #include <algorithm>
    using std::max_element;

    #include <cmath>
    using std::floor;

    //=====================
    //  Included libs
    //=====================

    #include "sparsepp.h"
    using spp::sparse_hash_map;

    #include "node.hpp"

    #include "relation_tables.hpp"

    //=====================

    int create_meme_file(vector<pair<const string, Node *> *> &, sparse_hash_map<string, Node *> &, sparse_hash_map<char, unsigned int> &, unsigned int, unsigned int, unsigned int, string &);

#endif
