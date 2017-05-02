#include "hash.hpp"


unsigned fill_hash_map( sparse_hash_map<string,
                        pair<int, Node *>> &encounters,
                        const string &filepath,
                        unsigned int l,
                        bool is_positive_file) {

    ifstream infile;
    infile.open(filepath);

    // checking if the file could be opened
    if (!infile.is_open()) {
        cerr << "Error the file could not be opened, make sure that \"" << filepath << "\" exists and isn't a directory" << endl;
        exit(EXIT_FAILURE);
    }

    unsigned global_motif_count = 0;
    string data;
    deque<char> deque(l);
    bool dequeReady = false;

    while(getline(infile, data)) {
        if(data.begin() == data.end() || *data.begin() == '>') {
            dequeReady = false;
            deque.clear();
        }
        else {

            /* -- Important note --
            The next for-loop is invalid for k=1, because the deque is set to ready at the end of the loop.
            But having k=1 is not supposed to happen since it means user wants to count single nucleotides,
            basically asking the count of A, T, C, G. It is not what the program is meant for.
            */

            for (char nuc : data) {
                //filtering accepted characters in fasta file
                switch(nuc) {
                    case 'A':
                    case 'C':
                    case 'G':
                    case 'T':
                        deque.emplace_back(nuc);
                        break;
                    //lowercases to be added as uppercases
                    case 'a':
                    case 't':
                    case 'g':
                    case 'c':
                        deque.emplace_back(toupper(nuc));
                        break;
                    //if character is invalid, emptying deque and refilling it
                    default:
                        deque.clear();
                        dequeReady = false;
                        break;
                }
                if (dequeReady) {
                    ++global_motif_count;
                    string motif(deque.begin(), deque.end());
                    auto motif_it = encounters.find(motif);
                    if (motif_it != encounters.end())
                        if (is_positive_file)
                            (motif_it->second).second->increment_positive_count();
                        else (motif_it->second).second->increment_negative_count();
                    else {
                        Node *new_node_ptr;
                        if (is_positive_file)
                            new_node_ptr = new Node(1, 0);
                        else new_node_ptr = new Node(0, 1);
                        // new_node_ptr->set_motif(motif);
                        encounters.emplace(make_pair(motif, make_pair(-1, new_node_ptr)));
                    }
                    deque.pop_front();
                }
                else
                    if (deque.size() == (l-1))
                        dequeReady = true;
            }
        }
    }
    infile.close();
    return global_motif_count;
}

unsigned fill_hash_map_from_pos(sparse_hash_map<string, pair<int, Node *>> &encounters,
                                         const string &filepath,
                                         unsigned int l,
                                         unsigned int p,
                                         bool is_positive_file) {
    ifstream infile;
    infile.open(filepath);

    // checking if the file could be opened
    if (!infile.is_open()) {
        cerr << "Error the file could not be opened, make sure that \"" << filepath << "\" exists and isn't a directory" << endl;
        exit(EXIT_FAILURE);
    }

    unsigned global_motif_count = 0;
    string data;
    deque<char> deque(l + p);

    while(getline(infile, data)) {
        if(data.begin() == data.end() || *data.begin() == '>') {
            if(on_sequence_end(deque, encounters, l, p, is_positive_file))
                ++global_motif_count;
            deque.clear();
        }
        else {
            for (char nuc : data) {
                if(deque.size() == l + p)
                    deque.pop_front();
                switch(nuc) {
                    case 'a':
                    case 'c':
                    case 'g':
                    case 't':
                        deque.emplace_back(toupper(nuc));
                        break;
                    default :
                        deque.emplace_back(nuc);
                        break;
                }
            }
        }
    }
    if (on_sequence_end(deque, encounters, l, p, is_positive_file)) {
        ++global_motif_count;
    }
    return global_motif_count;
}


bool on_sequence_end(deque<char> &deque,
                     sparse_hash_map<string, pair<int, Node *>> &encounters,
                     unsigned int l,
                     unsigned int p,
                     bool is_positive_file) {

    if(deque.size() == l + p) {
        string motif(deque.begin(), next(deque.begin(), l));
        if(motif.find_first_not_of("ACGT") != string::npos)
            return false;

        auto motif_it = encounters.find(motif);
        if (motif_it != encounters.end())
            if (is_positive_file)
                (motif_it->second).second->increment_positive_count();
            else (motif_it->second).second->increment_negative_count();
        else {
            Node *new_node_ptr;
            if (is_positive_file)
                new_node_ptr = new Node(1, 0);
            else new_node_ptr = new Node(0, 1);
            // new_node_ptr->set_motif(motif);
            encounters.emplace(make_pair(motif, make_pair(-1, new_node_ptr)));
        }
    }
    return true;
}
