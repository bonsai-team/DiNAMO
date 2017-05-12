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

    //=====================
    //  Included libs
    //=====================

    #include "lib/sparsepp.h"
    using spp::sparse_hash_map;

    #include "node.hpp"

    //=====================

    int create_meme_file(vector<pair<const string, pair<int, Node *> > *> &, sparse_hash_map<string, pair<int, Node *>> &, sparse_hash_map<char, unsigned int> &, unsigned int, string &);

#endif
