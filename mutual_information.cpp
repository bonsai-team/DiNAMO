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
    
    return  (pm / t) * log((pm / t) / (((pm + nm) / t) * (pt / t))) +
            (nm / t) * log((nm / t) / (((pm + nm) / t) * (nt / t))) +
            ((pt - pm) / t) * log(((pt - pm) / t) / ((((pt - pm) + (nt - nm)) / t) * (pt / t))) +
            ((nt - nm) / t) * log(((nt - nm) / t) / ((((pt - pm) + (nt - nm)) / t) * (nt / t)));
}
