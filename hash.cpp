#include "hash.hpp"

sparse_hash_map<char, char> nuc_to_complement {
    {'A','T'},
    {'C','G'},
    {'G','C'},
    {'T','A'}
};

string reverse_complement(string &motif) {
    string rv(motif);
    transform(rv.begin(), rv.end(), rv.begin(),
    [](char nuc) {return nuc_to_complement[nuc]; });
    reverse(rv.begin(), rv.end());
    // std::cout << motif << "\t" << rv << std::endl;
    return rv;
}

void open_file(ifstream &infile, const string &filepath) {
    infile.open(filepath);

    // checking if the file could be opened
    if (!infile.is_open()) {
        cerr << "Error the file could not be opened, make sure that \"" << filepath << "\" exists and isn't a directory" << endl;
        exit(EXIT_FAILURE);
    }
}

bool is_new_sequence(string &data) {
    return ((data.begin() == data.end()) || (*data.begin() == '>'));
}

void inc_global_count(unsigned &global_motif_count, bool rc) {
    if (rc)
        global_motif_count += 2;
    else global_motif_count++;
}

void on_valid_motif(string &motif,
                    sparse_hash_map<string, pair<int, Node *>> &encounters,
                    bool is_positive_file,
                    bool rc
                   ) {
    auto motif_it = encounters.find(motif);
    if (motif_it != encounters.end())
        if (is_positive_file) {
            (motif_it->second).second->increment_positive_count();
            if (rc && reverse_complement(motif) == motif)
                (motif_it->second).second->increment_positive_count();
        }
        else {
            (motif_it->second).second->increment_negative_count();
            if (rc && reverse_complement(motif) == motif)
                (motif_it->second).second->increment_negative_count();
        }
    else {
        Node *new_node_ptr;
        if (is_positive_file)
            new_node_ptr = new Node(1, 0);
        else new_node_ptr = new Node(0, 1);
        // new_node_ptr->set_motif(motif);
        encounters.emplace(make_pair(motif, make_pair(-1, new_node_ptr)));
        // case of reverse complements
        if (rc) {
            string rv_cmp(reverse_complement(motif));
            if (rv_cmp != motif)
                encounters.emplace(make_pair(rv_cmp, make_pair(-1, new_node_ptr)));
            else if (is_positive_file)
                new_node_ptr->increment_positive_count();
            else new_node_ptr->increment_negative_count();
        }
    }
}

unsigned fill_hash_map( sparse_hash_map<string, pair<int, Node *>> &encounters,
                        const string &filepath,
                        unsigned int l,
                        bool is_positive_file,
                        bool rc
                        ) {
    ifstream infile;
    open_file(infile, filepath);

    unsigned global_motif_count = 0;
    string data;
    deque<char> deque(l);
    bool dequeReady = (l == 1);

    while(getline(infile, data)) {
        if(is_new_sequence(data)) {
            dequeReady = false;
            deque.clear();
        }
        else {

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
                    //if character is invalid, emptying deque
                    default:
                        deque.clear();
                        dequeReady = false;
                        break;
                }
                if (dequeReady) {
                    inc_global_count(global_motif_count, rc);
                    string motif(deque.begin(), deque.end());
                    on_valid_motif(motif, encounters, is_positive_file, rc);
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

bool on_sequence_end(deque<char> &deque,
                     sparse_hash_map<string, pair<int, Node *>> &encounters,
                     unsigned int l,
                     unsigned int p,
                     bool is_positive_file,
                     bool rc
                     ) {

    if(deque.size() == l + p) {
        string motif(deque.begin(), next(deque.begin(), l));
        if(motif.find_first_not_of("ACGT") != string::npos)
            return false;

        on_valid_motif(motif, encounters, is_positive_file, rc);
    }
    return true;
}


unsigned fill_hash_map_from_pos(sparse_hash_map<string, pair<int, Node *>> &encounters,
                                const string &filepath,
                                unsigned int l,
                                unsigned int p,
                                bool is_positive_file,
                                bool rc
                                ) {
    ifstream infile;
    open_file(infile, filepath);

    unsigned global_motif_count = 0;
    string data;
    deque<char> deque(l + p);

    while(getline(infile, data)) {
        if(is_new_sequence(data)) {
            //returns true only if the motif inside the deque was valid
            if(on_sequence_end(deque, encounters, l, p, is_positive_file, rc)) {
                inc_global_count(global_motif_count, rc);
            }
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
    //returns true only if the motif inside the deque was valid
    if (on_sequence_end(deque, encounters, l, p, is_positive_file, rc)) {
        inc_global_count(global_motif_count, rc);
    }
    return global_motif_count;
}
