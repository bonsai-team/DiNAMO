#include "mutual_information.hpp"
#include <iostream>

double mutual_information(  unsigned positive_motif_count,        /////////// positif // negatif ////////////////////
                            unsigned negative_motif_count,        // motif //   pm    //   nm    //  pm + nm       //
                            unsigned positive_total_count,        //~motif // pt - pm // nt - nm //  pt-pm + nt-nm //
                            unsigned negative_total_count) {      // total //   pt    //   nt    //  t             //
    double pm = positive_motif_count;
    double nm = negative_motif_count;
    double pt = positive_total_count;
    double nt = negative_total_count;

    double t = positive_total_count + negative_total_count;

    double pvalue = (pm / t) * log2((pm / pt) * (t / (pm + nm))) +
                    (nm / t) * log2((nm / nt) * (t / (pm + nm))) +
                    ((pt - pm) / t) * log2(((pt - pm) / ((pt - pm) + (nt - nm))) * (t / pt)) +
                    ((nt - nm) / t) * log2(((nt - nm) / ((pt - pm) + (nt - nm))) * (t / nt));

    return pvalue;
}
