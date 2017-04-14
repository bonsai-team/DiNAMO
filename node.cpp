#include "node.hpp"

Node::Node(string &motif) {
    this->motif = string(motif);
    this->count = 1;
    this->in  = vector<Node *>();
    this->out = vector<Node *>();
}

Node::Node(string &motif, unsigned int count) {
    this->motif = string(motif);
    this->count = count;
    this->in  = vector<Node *>();
    this->out = vector<Node *>();
}

void Node::add_in_connection(Node *child) {
    this->in.push_back(child);
}

void Node::add_out_connection(Node *parent) {
    this->out.push_back(parent);
}

unsigned int Node::get_count() {
    return this->count;
}

void Node::increment_count() {
    ++(this->count);
}
