#include "meme_format.hpp"

void round_to_hundred_thousandth( double &d ) {
    d = (double) floor((d * 100000.0 ) + 0.5) / 100000;
}

void normalize_freq(vector<double> &freq_vec) {

    double freq_sum = 0.0;
    for (auto &freq : freq_vec) {
        round_to_hundred_thousandth(freq);
        freq_sum += freq;
    }
    double diff = 1.0 - freq_sum;

    *max_element(freq_vec.begin(), freq_vec.end()) += diff;
}


void generate_exact_successor_motifs_rec(vector<string> &res, const string &motif, string &acc, unsigned int pos, unsigned int l) {
    if (pos == l) {
        res.emplace_back(acc);
        return;
    }
    for (auto nuc : iupac_to_nucs[motif[pos]]) {
        acc.push_back(nuc);
        generate_exact_successor_motifs_rec(res, motif, acc, pos+1, l);
        acc.pop_back();
    }
}

void generate_exact_successor_motifs(vector<string> &res, const string &motif, unsigned int l) {
    string acc;
    generate_exact_successor_motifs_rec(res, motif, acc, 0, l);
}

int create_meme_file(vector<pair<const string, Node *> *> &entries,
                     sparse_hash_map<string, Node *> &exact_motifs,
                     sparse_hash_map<char, unsigned int> &neg_nuc_count,
                     unsigned int global_motif_count_positive,
                     unsigned int global_motif_count_negative,
                     unsigned int l,
                     string &filename) {

    ofstream meme_file;

    meme_file.open(filename, std::ios::out | std::ios::trunc);
    // checking if the file could be opened
    if (!meme_file.is_open()) {
        cerr << "Error opening meme file, make sure that the program has write access to \"" << filename << "\"." << endl;
        cerr << "Meme file creation has failed" << endl;
        return 1;
    }

    double total_neg_nuc_count = 0;
    for (auto const &entry : neg_nuc_count) {
        total_neg_nuc_count += entry.second;
    }

    vector<double> background_frequencies{
        neg_nuc_count['A'] / total_neg_nuc_count,
        neg_nuc_count['C'] / total_neg_nuc_count,
        neg_nuc_count['G'] / total_neg_nuc_count,
        neg_nuc_count['T'] / total_neg_nuc_count
    };

    normalize_freq(background_frequencies);

    //en-tête
    meme_file << "MEME version 4" << endl
              << endl
              << "ALPHABET= ACGT" << endl
              << endl
              << "strands: + -" << endl
              << endl
              << "Background letter frequencies" << endl
              << "A " << background_frequencies[0] <<
                " C " << background_frequencies[1] <<
                " G " << background_frequencies[2] <<
                " T " << background_frequencies[3] << endl
              << endl;

    //corps
    for (auto const &entry : entries) {
        meme_file << "MOTIF " << entry->first << endl << endl;
        meme_file << "#\tMI\tP-value\t#MotifPos/#TotalPos\t#MotifPos/#TotalPos" << endl;
        meme_file << "#\t" << entry->second->get_mi()
                  <<  "\t" << entry->second->get_pvalue()
                  <<  "\t" << entry->second->get_positive_count() << "/" << global_motif_count_positive
                  <<  "\t" << entry->second->get_negative_count() << "/" << global_motif_count_negative
                  << endl << endl;
        meme_file << "letter-probability matrix:" << endl;

        vector<string> exact_successor_motifs;
        generate_exact_successor_motifs(exact_successor_motifs, entry->first, l);
        for (size_t i=0; i < l; i++) {
            size_t exact_motifs_count = 0;
            sparse_hash_map<char, double> nuc_counts{
                {'A', 0.0},
                {'C', 0.0},
                {'G', 0.0},
                {'T', 0.0}
            };
            for (auto const &motif : exact_successor_motifs) {
                nuc_counts[motif[i]] += exact_motifs[motif]->get_positive_count();
                exact_motifs_count += exact_motifs[motif]->get_positive_count();
            }

            vector<double> motif_nuc_frequencies{
                nuc_counts['A'] / exact_motifs_count,
                nuc_counts['C'] / exact_motifs_count,
                nuc_counts['G'] / exact_motifs_count,
                nuc_counts['T'] / exact_motifs_count
            };

            normalize_freq(motif_nuc_frequencies);

            meme_file << motif_nuc_frequencies[0] << "\t"
                      << motif_nuc_frequencies[1] << "\t"
                      << motif_nuc_frequencies[2] << "\t"
                      << motif_nuc_frequencies[3] << endl;
        }
        meme_file << endl;
    }
    return 0;
}
