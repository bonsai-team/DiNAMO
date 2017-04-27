#include "stat_container.hpp"

StatContainer::StatContainer(string motif, double mi, double pvalue) {
    this->motif = motif;
    this->flag = unvisited;
    this->mi = mi;
    this->pvalue = pvalue;
}

double StatContainer::get_mi() {
    return this->mi;
}

string StatContainer::get_motif() {
    return this->motif;
}

bool StatContainer::compare_by_mi(const StatContainer &st1, const StatContainer &st2) {
    return (st1.mi) > (st2.mi);
}

bool StatContainer::compare_by_pvalue(const StatContainer &st1, const StatContainer &st2) {
    return (st1.pvalue) < (st2.pvalue);
}
