#if !defined(__GRAPH_SIMPLIFICATION_HPP_INCLUDED__)

    #define  __GRAPH_SIMPLIFICATION_HPP_INCLUDED__

    //=====================
    //     Dependencies
    //=====================

    #include <vector>
    using std::vector;

    #include <utility>
    using std::pair;

    #include <string>
    using std::string;

    //=====================
    //  Lib dependencies
    //=====================

    #include "node.hpp"

    //====================

    void graph_simplification(vector<pair<const string, Node *> *>, bool);

    void suppress_predecessors(Node *);

    void tag_successors(Node *);

    void tag_predecessors(Node *);

    void suppress_successors_and_tag_predecessors(Node *);

#endif
