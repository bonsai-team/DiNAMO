//========================
//         Guards
//========================
#if !defined __OPTIONSPARSER_H_INCLUDED__
#define      __OPTIONSPARSER_H_INCLUDED__

//========================
//  Included dependencies
//========================

#include <vector>
#include <string>

using std::string;
using std::vector;

//========================

class InputParser {

    private:
        vector<string> tokens;
        string empty_string;

    public:
        InputParser (int &, char **);

        const string& getCmdOption(const string &) const;

        bool cmdOptionExists(const string &) const;
};

#endif
