#include "graph_simplification.hpp"

void graph_simplification(vector<pair<const string, Node *> *> mi_sorted_entry_ptr, bool fixed_pos) {
    for (auto &entry_ptr : mi_sorted_entry_ptr) {
        Node *node_ptr = entry_ptr->second;
        switch (node_ptr->get_state()) {
        case unvisited:
            // std::cout << entry_ptr->first << " unvisited->validated " << node_ptr << endl;
            node_ptr->validate();
            for (Node *predecessor_ptr: node_ptr->get_predecessors()) {
                suppress_predecessors(predecessor_ptr);
            }
            for (Node *successor_ptr: node_ptr->get_successors()) {
                if(fixed_pos)
                    tag_successors(successor_ptr);
                else {
                    // std::cout << entry_ptr->first << " -X-> " << successor_ptr << endl;
                    suppress_successors_and_tag_predecessors(successor_ptr);
                }
            }
            break;
        case tagged:
            // std::cout << entry_ptr->first << " tagged->suppressed " << node_ptr << endl;
            node_ptr->suppress();
            for (Node *predecessor_ptr: node_ptr->get_predecessors()) {
                suppress_predecessors(predecessor_ptr);
            }
            break;
        default:
            // std::cout << entry_ptr->first << " already validated or supressed " << node_ptr << endl;
            continue;
        }
    }
}

void suppress_predecessors(Node *node_ptr) {
    switch (node_ptr->get_state()) {
    case unvisited:
    case tagged:
        node_ptr->suppress();
        for (Node *predecessor_ptr: node_ptr->get_predecessors()) {
            suppress_predecessors(predecessor_ptr);
        }
        break;
    default:
        return;
    }
}

void tag_successors(Node *node_ptr) {
    switch (node_ptr->get_state()) {
    case unvisited:
        node_ptr->tag();
        for (Node *successor_ptr: node_ptr->get_successors()) {
            // std::cout << node_ptr << " -T-> " << successor_ptr << endl;
            tag_successors(successor_ptr);
        }
        break;
    default:
        return;
    }
}

void tag_predecessors(Node *node_ptr) {
    switch (node_ptr->get_state()) {
    case unvisited:
        node_ptr->tag();
        for (Node *predecessor_ptr: node_ptr->get_predecessors()) {
            // std::cout << node_ptr << " -T-> " << predecessor_ptr << endl;
            tag_predecessors(predecessor_ptr);
        }
        break;
    default:
        return;
    }
}

void suppress_successors_and_tag_predecessors(Node *node_ptr) {
    switch (node_ptr->get_state()) {
    case unvisited:
    case tagged:
        node_ptr->suppress();
        for (Node *predecessor_ptr: node_ptr->get_predecessors()) {
            // std::cout << node_ptr << " -T-> " << predecessor_ptr << endl;
            tag_predecessors(predecessor_ptr);
        }
        for (Node *successor_ptr: node_ptr->get_successors()) {
            // std::cout << node_ptr << " -X-> " << successor_ptr << endl;
            suppress_successors_and_tag_predecessors(successor_ptr);
        }
        break;
    default:
        return;
    }
}
