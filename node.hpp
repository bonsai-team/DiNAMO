#if !defined(__NODE_HPP_INCLUDED__)

    #define __NODE_HPP_INCLUDED__

    //===========================
    //      Dependencies
    //===========================

    #include <string>
    using std::string;

    #include <vector>
    using std::vector;

    //===========================
    //      Libs needed
    //===========================


    //===========================

    class Node {

    private:
        string motif;
        unsigned int positive_count;
        unsigned int negative_count;
        vector<Node *> successors;
        vector<Node *> predecessors;

    public:
        Node(string &, unsigned int, unsigned int);
        void add_predecessor(Node *);
        void add_successor(Node *);
        unsigned int get_positive_count();
        unsigned int get_negative_count();
        void increment_positive_count();
        void increment_negative_count();
        void set_positive_count(unsigned int);
        void set_negative_count(unsigned int);
    };

#endif
