#include "meme_format.hpp"

sparse_hash_map<char, unordered_set<char>> iupac_or_nuc_to_nucs {
    {'A', {'A'}},
    {'C', {'C'}},
    {'G', {'G'}},
    {'T', {'T'}},
    {'R', {'A','G'}},
    {'Y', {'C','T'}},
    {'S', {'G','C'}},
    {'W', {'A','T'}},
    {'K', {'G','T'}},
    {'M', {'A','C'}},
    {'B', {'C','G','T'}},
    {'D', {'A','G','T'}},
    {'H', {'A','C','T'}},
    {'V', {'A','C','G'}},
    {'N', {'A','C','G','T'}}
};

void generate_exact_successor_motifs_rec(vector<string> &res, const string &motif, string &acc, unsigned int pos, unsigned int l) {
    if (pos == l) {
        res.emplace_back(acc);
        return;
    }
    for (auto nuc : iupac_or_nuc_to_nucs[motif[pos]]) {
        acc.push_back(nuc);
        generate_exact_successor_motifs_rec(res, motif, acc, pos+1, l);
        acc.pop_back();
    }
}

void generate_exact_successor_motifs(vector<string> &res, const string &motif, unsigned int l) {
    string acc;
    generate_exact_successor_motifs_rec(res, motif, acc, 0, l);
}

int create_meme_file(vector<pair<const string, pair<int, Node *> > *> &entries,
                     sparse_hash_map<string, pair<int, Node *>> &exact_motifs,
                     sparse_hash_map<char, unsigned int> &neg_nuc_count,
                     unsigned int l) {

    ofstream meme_file;

    string meme_filename("results.meme");

    meme_file.open(meme_filename, std::ios::out | std::ios::trunc);
    // checking if the file could be opened
    if (!meme_file.is_open()) {
        cerr << "Error opening meme file, make sure that the program has write access to \"" << meme_filename << "\"." << endl;
        cerr << "Meme file creation has failed" << endl;
        return 1;
    }

    double total_neg_nuc_count = 0;
    for (auto const &entry : neg_nuc_count) {
        total_neg_nuc_count += entry.second;
    }

    //en-tête
    meme_file << "MEME version 4" << endl
              << endl
              << "ALPHABET= ACGT" << endl
              << endl
              << "strands: + -" << endl
              << endl
              << "Background letter frequencies" << endl
              << "A " << neg_nuc_count['A'] / total_neg_nuc_count << " C " << neg_nuc_count['C'] / total_neg_nuc_count << " G " << neg_nuc_count['G'] / total_neg_nuc_count << " T " << neg_nuc_count['T'] / total_neg_nuc_count << endl
              << endl;

    //corps
    for (auto const &entry : entries) {
        meme_file << "MOTIF " << entry->first << endl << endl;
        meme_file << "#\tMI\tP-value" << endl;
        meme_file << "#\t" << entry->second.second->get_mi() << "\t" << entry->second.second->get_pvalue() << endl << endl;
        meme_file << "letter-probability matrix:" << endl;

        vector<string> exact_successor_motifs;
        generate_exact_successor_motifs(exact_successor_motifs, entry->first, l);
        for (size_t i=0; i < l; i++) {
            size_t exact_motifs_count = 0;
            sparse_hash_map<char, float> nuc_counts{
                {'A', 0.0},
                {'C', 0.0},
                {'G', 0.0},
                {'T', 0.0}
            };
            for (auto const &motif : exact_successor_motifs) {
                nuc_counts[motif[i]] += exact_motifs[motif].second->get_positive_count();
                exact_motifs_count += exact_motifs[motif].second->get_positive_count();
            }
            meme_file << nuc_counts['A'] / exact_motifs_count << "\t"
                      << nuc_counts['C'] / exact_motifs_count << "\t"
                      << nuc_counts['G'] / exact_motifs_count << "\t"
                      << 1.0 - nuc_counts['A'] / exact_motifs_count
                             - nuc_counts['C'] / exact_motifs_count
                             - nuc_counts['G'] / exact_motifs_count
                      << endl;
        }
        meme_file << endl;
    }
    return 0;
}
