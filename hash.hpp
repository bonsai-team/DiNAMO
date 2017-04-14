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

    //=====================
    //  Included libs
    //=====================

    #include "node.hpp"

    #include "lib/sparsepp.h"
    using spp::sparse_hash_map;

    //=====================

    void fill_hash_map(sparse_hash_map<string, pair<int, Node *>> &, vector<Node *> &, const string &, unsigned int);

    void fill_hash_map_from_pos(sparse_hash_map<string, pair<int, Node *>> &, vector<Node *> &, const string &, unsigned int, unsigned int);

    void on_sequence_end(deque<char> &, sparse_hash_map<string, pair<int, Node *>> &, vector<Node *> &, unsigned int, unsigned int);

#endif
