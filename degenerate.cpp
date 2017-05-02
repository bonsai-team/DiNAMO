#include "degenerate.hpp"

const vector<char> nucleotides {'A','C','G','T'};

sparse_hash_map<string, vector<unordered_set<char>>> nucs_to_iupacs {
    {"ACGT", {{'Y','S','K','R','W','M'}, {'B','D','H','V'}, {'N'}}},
    {"ACG" , {{'M','R','S'}, {'V'}}},
    {"ACT" , {{'M','W','Y'}, {'H'}}},
    {"AGT" , {{'R','W','K'}, {'D'}}},
    {"CGT" , {{'S','Y','K'}, {'B'}}},
    {"AC"  , {{'M'}}},
    {"AG"  , {{'R'}}},
    {"AT"  , {{'W'}}},
    {"CG"  , {{'S'}}},
    {"CT"  , {{'Y'}}},
    {"GT"  , {{'K'}}}
};

sparse_hash_map<char, unordered_set<char>> iupacs_dependencies {
    {'N', {'B','D','H','V'}},
    {'B', {'Y','S','K'}},
    {'D', {'R','W','K'}},
    {'H', {'Y','W','M'}},
    {'V', {'R','S','M'}},
    {'R', {'A','G'}},
    {'Y', {'C','T'}},
    {'S', {'G','C'}},
    {'W', {'A','T'}},
    {'K', {'G','T'}},
    {'M', {'A','C'}}
};

sparse_hash_map<char, unordered_set<char>> iupac_to_nucs {
    {'N', {'A', 'C', 'G', 'T'}},
    {'B', {'C', 'G', 'T'}},
    {'D', {'A', 'G', 'T'}},
    {'H', {'A', 'C', 'T'}},
    {'V', {'A', 'C', 'G'}},
    {'Y', {'C', 'T'}},
    {'S', {'C', 'G'}},
    {'K', {'G', 'T'}},
    {'R', {'A', 'G'}},
    {'W', {'A', 'T'}},
    {'M', {'A', 'C'}}
};

//retourne les nucléotides concaténés dans l'ordre (assuré par une itération sur le vecteur nucleotides)
const string find_neighbor_motifs(  sparse_hash_map<string, pair<int, Node *>> &motifs,
                                    sparse_hash_map<char, pair<string, Node *>> &neighbor_motifs,
                                    const string &motif,
                                    unsigned int pos) {

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
        motif_it->second.first = pos;
        //sinon, on insère dans la table des résultats le motif et son compte
        neighbor_motifs[nuc] = make_pair(neighbor_motif, motifs[neighbor_motif].second);
        neighbors_nuc.push_back(nuc);
    }
    return neighbors_nuc;
}

void degenerate(sparse_hash_map<string, pair<int, Node *>> &motifs,
                sparse_hash_map<string, pair<int, Node *>> &degenerated_motifs,
                const unsigned int kmer_size) {

    //pour chaque position
    for (unsigned int pos = 0; pos < kmer_size; pos++) {

        //pour chaque motif dans la table
        for (auto const &motif_it : motifs) {
            //s'il n'a pas été marqué pour cette position ou si la position a déjà été dégénérée
            if(motif_it.second.first == pos || !string("ACGT").find(motif_it.first[pos])) {
                continue;
            }

            sparse_hash_map<char, pair<string, Node *>> neighbor_motifs;
            neighbor_motifs.reserve(4);
            //on trouve les voisins du kmer par rapport à la position courante
            const string neighbors_nuc = find_neighbor_motifs(motifs, neighbor_motifs, motif_it.first, pos);

            if (neighbors_nuc.size() <= 1) {
                continue;
            }
            //this is where things can explode if the string is not sorted alphabetically
            auto iupacs = nucs_to_iupacs.find(neighbors_nuc);
            assert(iupacs != nucs_to_iupacs.end());

            string degenerated_motif(motif_it.first);

            //TODO not required since lookup in a hashtable is O(1)
            sparse_hash_map<char, Node *> iupac_to_node;
            iupac_to_node.reserve(11);

            // pour tous les iupacs dérivant des nucléotides
            for (int degeneration_degree = 0; degeneration_degree < iupacs->second.size(); degeneration_degree++) {
                for (auto const iupac : iupacs->second[degeneration_degree]) {

                    degenerated_motif.replace(pos, 1, 1, iupac);
                    Node *current_node_ptr;

                    //checks if the node has already been created
                    auto entry_ptr = degenerated_motifs.find(degenerated_motif);
                    if(entry_ptr != degenerated_motifs.end()) {
                        current_node_ptr = entry_ptr->second.second;
                    }
                    else {
                        //node creation
                        current_node_ptr = new Node(0, 0);
                        unsigned int degenerated_motif_positive_count = 0;
                        unsigned int degenerated_motif_negative_count = 0;
                        for (auto const &nuc : iupac_to_nucs[iupac]) {
                            //add child count
                            degenerated_motif_positive_count += neighbor_motifs[nuc].second->get_positive_count();
                            degenerated_motif_negative_count += neighbor_motifs[nuc].second->get_negative_count();
                        }
                        current_node_ptr->set_positive_count(degenerated_motif_positive_count);
                        current_node_ptr->set_negative_count(degenerated_motif_negative_count);
                        // current_node_ptr->set_motif(degenerated_motif);
                        degenerated_motifs.emplace(make_pair(degenerated_motif, make_pair(-1, current_node_ptr)));
                    }
                    //remembering the node that we created to be able to find them easily
                    //TODO this is not required since lookup in a hash table is O(1)
                    iupac_to_node[iupac] = current_node_ptr;

                    //adding links between nodes

                    //here we create the nodes that depends on the nucleotides first,
                    //then we create the nodes that depends on the previous one, etc.
                    if (degeneration_degree == 0) {
                        for (char nuc : iupacs_dependencies[iupac]) {
                            // std::cout << degenerated_motif  << "<->" << neighbor_motifs[nuc].second->get_motif() << "\t Position " << pos << endl;
                            current_node_ptr->add_successor(neighbor_motifs[nuc].second);
                            neighbor_motifs[nuc].second->add_predecessor(current_node_ptr);
                        }
                    } else {
                        for (char iupac_dependency : iupacs_dependencies[iupac]) {
                            // std::cout << degenerated_motif  << "<->" << iupac_to_node[iupac_dependency]->get_motif() << "\t Position " << pos << endl;
                            current_node_ptr->add_successor(iupac_to_node[iupac_dependency]);
                            iupac_to_node[iupac_dependency]->add_predecessor(current_node_ptr);
                        }
                    }
                }
            }
        }
    }
}
