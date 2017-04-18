#include "node.hpp"

Node::Node(string &motif, unsigned int positive_count, unsigned int negative_count) {
    this->motif = string(motif);
    this->positive_count = positive_count;
    this->negative_count = negative_count;
    this->successors  = vector<Node *>();
    this->predecessors = vector<Node *>();
}

void Node::add_predecessor(Node *child) {
    this->predecessors.push_back(child);
}

void Node::add_successor(Node *parent) {
    this->successors.push_back(parent);
}

unsigned int Node::get_positive_count() {
    return this->positive_count;
}

unsigned int Node::get_negative_count() {
    return this->negative_count;
}

void Node::increment_positive_count() {
    ++(this->positive_count);
}

void Node::increment_negative_count() {
    ++(this->negative_count);
}

void Node::set_positive_count(unsigned int positive_count) {
    this->positive_count = positive_count;
}

void Node::set_negative_count(unsigned int negative_count) {
    this->negative_count = negative_count;
}
