//=====================
//      Guards
//=====================

#if !defined(__REVERSE_COMPLEMENT_HPP_INCLUDED__)

    #define  __REVERSE_COMPLEMENT_HPP_INCLUDED__

    //=====================
    //Included dependencies
    //=====================

    #include <algorithm>
    using std::transform;
    using std::reverse;

    #include <string>
    using std::string;

    //=====================
    //  Included libs
    //=====================

    #include "lib/sparsepp.h"
    using spp::sparse_hash_map;

    //=====================

    string reverse_complement(string &);

#endif
