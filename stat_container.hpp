#if !defined(__STAT_CONTAINER_HPP_INCLUDED__)

    #define  __STAT_CONTAINER_HPP_INCLUDED__

    //=====================
    //     Dependencies
    //=====================

    #include <string>
    using std::string;

    //=====================
    //  Lib dependencies
    //=====================



    //====================

    enum Flag {unvisited, successor, predecessor};

    class StatContainer {

    private:
        string motif;
        Flag flag;
        double mi;
        double pvalue;

    public:
        StatContainer(string, double, double);

        static bool compare_by_mi(const StatContainer &, const StatContainer &);
        static bool compare_by_pvalue(const StatContainer &, const StatContainer &);

        double get_mi();
        string get_motif();
    };

#endif
