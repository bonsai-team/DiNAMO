#include "fisher_test.hpp"

double fisher_test_p_value( unsigned positive_motif_count,
                            unsigned negative_motif_count,
                            unsigned positive_total_count,
                            unsigned negative_total_count,
                            Methods method) {

    hypergeometric_distribution<> hgd ( positive_total_count, //r (number of defective inside the total population)
                                        positive_motif_count + negative_motif_count, //n (sample ie. occurrences of the motif inside both files)
                                        positive_total_count + negative_total_count);//N (total population)
    /*
    the cumulative distribution function of the hypergeometric distribution
    allows us to know the probability that k < positive_motif_count.
    (ie. the sum of p(k), k < positive_motif_count )
    This is already implemented in boost, and gives us the traditional pvalue.
    */
    if (method == one_sided_less) {
        return cdf(hgd, positive_motif_count);
    }


    if (method == one_sided_greater) {
        return cdf(complement(hgd, positive_motif_count));
    }

    /*
    To match R fisher.test function, p(k = positive_motif_count) is chosen as the critical limit
    */
    if (method == two_tailed) {
        return cdf(hgd, positive_motif_count) + cdf(complement (hgd, quantile(complement(hgd, pdf(hgd, positive_motif_count)))));
    }
}
