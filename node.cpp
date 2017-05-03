#include "node.hpp"

Node::Node(unsigned int positive_count, unsigned int negative_count) {
    this->positive_count = positive_count;
    this->negative_count = negative_count;
    this->successors  = vector<Node *>();
    this->predecessors = vector<Node *>();
    this->state = unvisited;
}

void Node::add_predecessor(Node *child) {
    this->predecessors.push_back(child);
}

void Node::add_successor(Node *parent) {
    this->successors.push_back(parent);
}

vector<Node *> &Node::get_predecessors() {
    return this->predecessors;
}

vector<Node *> &Node::get_successors() {
    return this->successors;
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

void Node::calculate_mi(unsigned int global_motif_count_positive, unsigned int global_motif_count_negative) {
    this->mi = mutual_information (this->get_positive_count(),
                                   this->get_negative_count(),
                                   global_motif_count_positive,
                                   global_motif_count_negative);
}

void Node::calculate_pvalue(unsigned int global_motif_count_positive, unsigned int global_motif_count_negative) {
    this->pvalue = fisher_test_p_value(this->get_positive_count(),
                                       this->get_negative_count(),
                                       global_motif_count_positive,
                                       global_motif_count_negative,
                                       one_sided_greater);
}

double Node::get_mi() {
    return this->mi;
}

double Node::get_pvalue() {
    return this->pvalue;
}

void Node::validate() {
    this->state = validated;
}

void Node::suppress() {
    this->state = deleted;
}

void Node::tag() {
    this->state = tagged;
}

State Node::get_state() {
    return this->state;
}

// void Node::set_motif(string motif) {
//     this->motif = motif;
// }
//
// string &Node::get_motif() {
//     return this->motif;
// }
