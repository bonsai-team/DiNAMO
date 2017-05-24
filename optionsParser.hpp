//========================
//         Guards
//========================
#if !defined(__OPTIONSPARSER_HPP_INCLUDED__)

    #define      __OPTIONSPARSER_HPP_INCLUDED__

    //========================
    //  Included dependencies
    //========================

    #include <vector>
    using std::vector;

    #include <string>
    using std::string;

    #include <algorithm>
    using std::find;

    //========================

    class InputParser {

        private:
            vector<string> tokens;
            string empty_string;

        public:
            InputParser (int &, char **);

            const string& getCmdOption(const std::initializer_list<const string> &) const;

            bool cmdOptionExists(const std::initializer_list<const string> &) const;
    };

#endif
