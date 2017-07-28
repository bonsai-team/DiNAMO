#include "degenerate.hpp"

const string find_neighbours_at_pos(sparse_hash_map<string, Node *> &,
                                    sparse_hash_map<char, Node *> &,
                                    const string,
                                    const unsigned int);

Node *create_successor( Node *,
                        const string &,
                        sparse_hash_map<string, Node *> &,
                        sparse_hash_map<string, Node *> &,
                        bool);

void create_links(Node *, Node *);

void create_graph_from_node(sparse_hash_map<string, Node *> &,
                            sparse_hash_map<string, Node *> &,
                            const string,
                            Node *,
                            bool);

void successor_motifs_at_pos(const string &, unsigned int, vector<string> &);

Node *node_creation(sparse_hash_map<string, Node *> &,
                    sparse_hash_map<char, Node *> &,
                    const string &,
                    bool);

/*      Function prototypes            */

void degenerate(sparse_hash_map<string, Node *> &motifs,
                sparse_hash_map<string, Node *> &degenerated_motifs,
                const unsigned int motif_size,
                const bool rc) {

    for (unsigned int pos = 0; pos < motif_size; pos++) {
        for (auto const &motif_it : motifs){
            if (motif_it.first[pos] != 'A' && motif_it.first[pos] != 'C' && motif_it.first[pos] != 'G' && motif_it.first[pos] != 'T')
                continue;

            const string motif = motif_it.first;
            //find_neighbours at current position
            sparse_hash_map<char, Node *> neighbours;
            neighbours.reserve(4);
            const string neighbours_nuc = find_neighbours_at_pos(motifs, neighbours, motif, pos);


	    
            if (neighbours_nuc.size() <= 1) {
	      continue;
            }


            //generate appropriate degenerated motif
            //this needs the string to be sorted alphabetically
            auto generated_iupac_entry = nucs_to_iupac.find(neighbours_nuc);
            assert(generated_iupac_entry != nucs_to_iupac.end());
            char generated_iupac = generated_iupac_entry->second;

            string degenerated_motif(motif);
            degenerated_motif.replace(pos, 1, 1, generated_iupac);

            if (degenerated_motifs.find(degenerated_motif) != degenerated_motifs.end())
                continue;

            Node *degenerated_motif_ptr = node_creation(degenerated_motifs, neighbours, degenerated_motif, rc);

            create_graph_from_node(motifs, degenerated_motifs, degenerated_motif, degenerated_motif_ptr, rc);
        }
    }
}

const string find_neighbours_at_pos(sparse_hash_map<string, Node *> &motifs,
                                    sparse_hash_map<char, Node *> &neighbours,
                                    const string motif,
                                    const unsigned int pos) {

    string neighbours_nuc; //each nucleotides found on neighbours at the current pos
    string hypothetical_neighbour = motif; //string that we will transform to look for neighbours in the hash table
    for (char nuc : nucleotides) {
        hypothetical_neighbour.replace(pos, 1, 1, nuc);
        auto const neighbour_ptr = motifs.find(hypothetical_neighbour); // checking if the neighbour is present in the hash table
        if (neighbour_ptr != motifs.end()) {
            neighbours_nuc.push_back(nuc);
            neighbours[nuc] = neighbour_ptr->second;
        }
    }
    return neighbours_nuc;
}

Node *create_successor( Node *predecessor_ptr,
                        const string &degenerated_motif,
                        sparse_hash_map<string, Node *> &motifs,
                        sparse_hash_map<string, Node *> &degenerated_motifs,
                        bool rc) {

    auto const entry_ptr = degenerated_motifs.find(degenerated_motif);
    if (entry_ptr != degenerated_motifs.end()) {
        Node *degenerated_motif_ptr = (*entry_ptr).second;
        create_links(predecessor_ptr, degenerated_motif_ptr);
    } else {
        //finding neighbours
        sparse_hash_map<char, Node *> neighbours;
        neighbours.reserve(4);
        unsigned int pos = degenerated_motif.find_first_not_of("ACGT");
        string neighbour(degenerated_motif);
        for (const char nuc : iupac_to_nucs[degenerated_motif[pos]]){
            neighbour.replace(pos, 1, 1, nuc);
            neighbours[nuc] = motifs[neighbour];
        }

        //creating the new node
        Node *degenerated_motif_ptr = node_creation(degenerated_motifs, neighbours, degenerated_motif, rc);

        //creating the link between the newly created node and its predecessor
        create_links(predecessor_ptr, degenerated_motif_ptr);

        create_graph_from_node(motifs, degenerated_motifs, degenerated_motif, degenerated_motif_ptr, rc);
    }
}

void create_links(Node *predecessor_ptr, Node *successor_ptr) {
    successor_ptr->add_predecessor(predecessor_ptr);
    predecessor_ptr->add_successor(successor_ptr);
}

void create_graph_from_node(sparse_hash_map<string, Node *> &motifs,
                            sparse_hash_map<string, Node *> &degenerated_motifs,
                            const string degenerated_motif,
                            Node *degenerated_motif_ptr,
                            bool rc) {

    for (unsigned int pos = 0; pos < degenerated_motif.size(); pos++) {
        char current_iupac = degenerated_motif[pos];
        switch (current_iupac) {

        case 'A':
        case 'C':
        case 'G':
        case 'T':
            continue;
            break;

        //cases where successors may not exist yet
        case 'B':
        case 'D':
        case 'H':
        case 'V':
        case 'N':
            {
                vector<string> successors_at_pos(4);
                successor_motifs_at_pos(degenerated_motif, pos, successors_at_pos);
                for (string successor_motif : successors_at_pos) {
                    create_successor(degenerated_motif_ptr, successor_motif, motifs, degenerated_motifs, rc);
                }
                break;
            }

        //cases where we are sure that sucessors exist
        default:
            {
                vector<string> successors_at_pos(4);
                successor_motifs_at_pos(degenerated_motif, pos, successors_at_pos);
                for (string successor_motif : successors_at_pos) {
                    Node *successor_ptr = motifs[successor_motif];
                    create_links(degenerated_motif_ptr, successor_ptr);
                }
                break;
            }
        }
    }
}

void successor_motifs_at_pos(const string &degenerated_motif, unsigned int pos, vector<string> &successors_at_pos) {
    successors_at_pos.clear();
    string successor(degenerated_motif);
    for (const char iupac : iupacs_dependencies[degenerated_motif[pos]]) {
        successor.replace(pos, 1, 1, iupac);
        successors_at_pos.push_back(successor);
    }
}

Node *node_creation(sparse_hash_map<string, Node *> &degenerated_motifs,
                    sparse_hash_map<char, Node *> &neighbours,
                    const string &degenerated_motif,
                    bool rc) {

    unsigned int degenerated_motif_positive_count = 0;
    unsigned int degenerated_motif_negative_count = 0;

    for (auto const &entry_ptr : neighbours) {

        unsigned int successor_positive_count = entry_ptr.second->get_positive_count();
        unsigned int successor_negative_count = entry_ptr.second->get_negative_count();

        //add positive successor count or raise flag
        if (degenerated_motif_positive_count > ~0 - successor_positive_count) {
            std::cerr << "Error : an overflow occurred while setting the positive count of a degenerated motif. You should consider switching to a bigger unsigned type." << endl;
            exit(EXIT_FAILURE);
        }
        degenerated_motif_positive_count += successor_positive_count;

        //add negative successor count or raise flag
        if (degenerated_motif_negative_count > ~0 - successor_negative_count) {
            std::cerr << "Error : an overflow occurred while setting the negative count of a degenerated motif. You should consider switching to a bigger unsigned type." << endl;
            exit(EXIT_FAILURE);
        }
        degenerated_motif_negative_count += successor_negative_count;
    }
    Node *current_node_ptr = new Node(degenerated_motif_positive_count, degenerated_motif_negative_count);
    degenerated_motifs.emplace(make_pair(degenerated_motif, current_node_ptr));
    if (rc) {
        string degenerated_motif_rc = reverse_complement(degenerated_motif);
        if (degenerated_motif_rc != degenerated_motif)
            degenerated_motifs.emplace(make_pair(degenerated_motif_rc, current_node_ptr));
    }
    return current_node_ptr;
}
