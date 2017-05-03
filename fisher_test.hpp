#if !defined(__FISHER_TEST_HPP_INCLUDED__)

    #define  __FISHER_TEST_HPP_INCLUDED__

    //=====================
    //     Dependencies
    //=====================

    #include <iostream>
    using std::cerr;
    using std::endl;

    //=====================
    //  Lib dependencies
    //=====================

    #include <boost/math/distributions/hypergeometric.hpp>
    using boost::math::hypergeometric_distribution;
    using boost::math::cdf;
    using boost::math::pdf;
    using boost::math::complement;
    using boost::math::quantile;

    //====================

    enum Methods {one_sided_less, one_sided_greater, two_tailed};

    double fisher_test_p_value(unsigned int, unsigned int, unsigned int, unsigned int, Methods);

#endif
