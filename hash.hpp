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

    unsigned long long int fill_hash_map_positive(sparse_hash_map<string, pair<int, Node *>> &, vector<Node *> &, const string &, unsigned int);
    unsigned long long int fill_hash_map_negative(sparse_hash_map<string, pair<int, Node *>> &, vector<Node *> &, const string &, unsigned int);

    unsigned long long int fill_hash_map_from_pos_positive(sparse_hash_map<string, pair<int, Node *>> &, vector<Node *> &, const string &, unsigned int, unsigned int);
    unsigned long long int fill_hash_map_from_pos_negative(sparse_hash_map<string, pair<int, Node *>> &, vector<Node *> &, const string &, unsigned int, unsigned int);

    bool on_sequence_end_positive(deque<char> &, sparse_hash_map<string, pair<int, Node *>> &, vector<Node *> &, unsigned int, unsigned int);
    bool on_sequence_end_negative(deque<char> &, sparse_hash_map<string, pair<int, Node *>> &, vector<Node *> &, unsigned int, unsigned int);

#endif
