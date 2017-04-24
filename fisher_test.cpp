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
    if (method == one_tailed)
        return cdf(hgd, positive_motif_count);

    /*
    Another approach is to select the values by their odds ratio.
    If it is greater than the actual distribution's odds, add their probability to it.
    Odds ratio = positive_motif_count * (negative_total_count - negative_motif_count)
                 --------------------------------------------------------------------
                 negative_motif_count * (positive_total_count - positive_motif_count)
    */
    if (method == odds_ratio) {
        std::cerr << "Fisher exact test : odds ratio method isn't implemented yet" << std::endl;
        exit(EXIT_FAILURE);
    }

    /*
    To match R fisher.test function, p(k = positive_motif_count) is chosen as the critical limit
    */
    if (method == two_tailed) {
        return cdf(hgd, positive_motif_count) + cdf(complement (hgd, quantile(complement(hgd, pdf(hgd, positive_motif_count)))));
    }
}
