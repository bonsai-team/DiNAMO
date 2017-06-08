#if !defined(__NODE_HPP_INCLUDED__)

    #define __NODE_HPP_INCLUDED__

    //===========================
    //      Dependencies
    //===========================

    // #include <string>
    // using std::string;

    #include <vector>
    using std::vector;

    //===========================
    //      Libs needed
    //===========================

    #include "mutual_information.hpp"
    #include "fisher_test.hpp"

    //===========================

    enum State {validated = 0, unvisited = 1, tagged = 2, deleted = 3};

    class Node {

    private:
        // string motif;
        unsigned int positive_count;
        unsigned int negative_count;
        double mi;
        double pvalue;
        State state;
        vector<Node *> successors;
        vector<Node *> predecessors;

    public:
        Node(unsigned int, unsigned int);

        void add_predecessor(Node *);
        void add_successor(Node *);
        vector<Node *> &get_predecessors();
        vector<Node *> &get_successors();

        void increment_positive_count();
        void increment_negative_count();
        void set_positive_count(unsigned int);
        void set_negative_count(unsigned int);
        unsigned int get_positive_count();
        unsigned int get_negative_count();

        void calculate_mi(unsigned int, unsigned int);
        double get_mi();
        void calculate_pvalue(unsigned int, unsigned int);
        double get_pvalue();

        void reset_state();
        void validate();
        void suppress();
        void tag();
        State get_state();
    };

#endif
