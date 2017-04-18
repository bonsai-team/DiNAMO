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
        int count;
        vector<Node *> in;
        vector<Node *> out;

    public:
        Node(string &);
        void add_in_connection(Node *);
        void add_out_connection(Node *);
        unsigned int get_count();
        void increment_count();
        void set_count(unsigned int);
    };

#endif
