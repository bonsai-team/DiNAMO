#if !defined(__GRAPH_SIMPLIFICATION_HPP_INCLUDED__)

    #define  __GRAPH_SIMPLIFICATION_HPP_INCLUDED__

    //=====================
    //     Dependencies
    //=====================

    #include <vector>
    using std::vector;

    #include <utility>
    using std::pair;

    //=====================
    //  Lib dependencies
    //=====================

    #include "node.hpp"

    //====================

    void graph_simplification(vector<pair<const string, pair<int, Node *> > *>, bool);

    void suppress_predecessors(Node *);

    void tag_successors(Node *);

    void suppress_successors(Node *);

#endif
