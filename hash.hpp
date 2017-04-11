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

    //=====================
    //  Included libs
    //=====================

    #include "lib/sparsepp.h"
    using spp::sparse_hash_map;

    //=====================

    void fill_hash_map(sparse_hash_map<string, pair<int, int>> &, const string &, unsigned int);

#endif
