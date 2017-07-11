#include "degenerate.hpp"

//retourne les nucléotides concaténés dans l'ordre (assuré par une itération sur le vecteur nucleotides)
const string find_neighbor_motifs(  sparse_hash_map<string, Node *> &motifs,
                                    sparse_hash_map<char, pair<string, Node *>> &neighbor_motifs,
                                    const string &motif,
                                    unsigned int pos,
                                    sparse_hash_map<string, int> &last_checked_pos) {

    string neighbor_motif(motif);
    string neighbors_nuc;
    neighbors_nuc.reserve(4);
    //pour tous les nucléotides
    for (auto const &nuc : nucleotides) {
        //on remplace la (pos)ième position du motif de base par ce nucléotide
        neighbor_motif.replace(pos, 1, 1, nuc);

        //on passe au suivant si le voisin n'est pas dans la table
        auto const motif_it = motifs.find(neighbor_motif);
        if(motif_it == motifs.end())
            continue;
        last_checked_pos[motif_it->first] = pos;
        //sinon, on insère dans la table des résultats le motif et son compte
        neighbor_motifs[nuc] = *motif_it;//make_pair(neighbor_motif, motifs[neighbor_motif].second);
        neighbors_nuc.push_back(nuc);
    }
    return neighbors_nuc;
}

void degenerate(sparse_hash_map<string, Node *> &motifs,
                sparse_hash_map<string, Node *> &degenerated_motifs,
                const unsigned int kmer_size,
                bool rc) {

    sparse_hash_map<string, int> last_checked_pos;
    last_checked_pos.reserve(motifs.size());
    for (auto const &motif_it : motifs) {
        if (motif_it.second->get_state() == unvisited) {
            last_checked_pos[motif_it.first] = -1;
            motif_it.second->validate();
        }
        else
            last_checked_pos[motif_it.first] = kmer_size;
    }

    //pour chaque position
    for (unsigned int pos = 0; pos < kmer_size; pos++) {

        //pour chaque motif dans la table
        for (auto const &motif_it : motifs) {
            //s'il n'a pas été marqué pour cette position ou si la position a déjà été dégénérée
            if((last_checked_pos[motif_it.first] == kmer_size) || (last_checked_pos[motif_it.first] == pos) || (!string("ACGT").find(motif_it.first[pos]))) {
                continue;
            }

            sparse_hash_map<char, pair<string, Node *>> neighbor_motifs;
            neighbor_motifs.reserve(4);
            //on trouve les voisins du kmer par rapport à la position courante
            const string neighbors_nuc = find_neighbor_motifs(motifs, neighbor_motifs, motif_it.first, pos, last_checked_pos);

            if (neighbors_nuc.size() <= 1) {
                continue;
            }
            //this is where things can explode if the string is not sorted alphabetically
            auto iupacs = nucs_to_iupacs.find(neighbors_nuc);
            assert(iupacs != nucs_to_iupacs.end());

            //string that will be transformed into the degenerated motif later
            string degenerated_motif(motif_it.first);

            //TODO not required since lookup in a hashtable is O(1)
            sparse_hash_map<char, Node *> iupac_to_node;
            iupac_to_node.reserve(11);

            // pour tous les iupacs dérivant des nucléotides
            for (unsigned int degeneration_degree = 0; degeneration_degree < iupacs->second.size(); degeneration_degree++) {
                for (auto const iupac : iupacs->second[degeneration_degree]) {

                    degenerated_motif.replace(pos, 1, 1, iupac);
                    Node *current_node_ptr;

                    //checks if the node has already been created
                    auto entry_ptr = degenerated_motifs.find(degenerated_motif);
                    if(entry_ptr != degenerated_motifs.end()) {
                        current_node_ptr = entry_ptr->second;
                    }
                    else {
                        //node creation
                        unsigned int degenerated_motif_positive_count = 0;
                        unsigned int degenerated_motif_negative_count = 0;

                        for (auto const &nuc : iupac_to_nucs[iupac]) {

                            unsigned int successor_positive_count = neighbor_motifs[nuc].second->get_positive_count();
                            unsigned int successor_negative_count = neighbor_motifs[nuc].second->get_negative_count();

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
                        current_node_ptr = new Node(degenerated_motif_positive_count, degenerated_motif_negative_count);
                        // current_node_ptr->set_motif(degenerated_motif);
                        degenerated_motifs.emplace(make_pair(degenerated_motif, current_node_ptr));
                        if (rc) {
                            degenerated_motifs.emplace(make_pair(reverse_complement(degenerated_motif), current_node_ptr));
                        }
                    }
                    //remembering the node that we created to be able to find them easily
                    //TODO this is not required since lookup in a hash table is O(1)
                    iupac_to_node[iupac] = current_node_ptr;

                    //adding links between nodes

                    //here we create the nodes that depends on the nucleotides first,
                    //then we create the nodes that depends on the previous one, etc.
                    if (degeneration_degree == 0) {
                        for (char nuc : iupacs_dependencies[iupac]) {
                            current_node_ptr->add_successor(neighbor_motifs[nuc].second);
                            neighbor_motifs[nuc].second->add_predecessor(current_node_ptr);
                        }
                    } else {
                        for (char iupac_dependency : iupacs_dependencies[iupac]) {
                            if (iupac_to_node.find(iupac_dependency) == iupac_to_node.end()) {
                                std::cerr << "Tried to link a node to an non-existing one" << std::endl;
                                exit(EXIT_FAILURE);
                            }
                            current_node_ptr->add_successor(iupac_to_node[iupac_dependency]);
                            iupac_to_node[iupac_dependency]->add_predecessor(current_node_ptr);
                        }
                    }
                }
            }
        }
    }
}
