#include <algorithm>
#include "optionsParser.hpp"


InputParser::InputParser(int &argc, char **argv){
    for (int i=1; i < argc; ++i)
        this->tokens.emplace_back(string(argv[i]));
}

const string& InputParser::getCmdOption(const string &option) const{
    vector<string>::const_iterator itr;
    itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
    if (itr != this->tokens.end() && ++itr != this->tokens.end()){
        return *itr;
    }
    return empty_string;
}

bool InputParser::cmdOptionExists(const string &option) const{
    return std::find(this->tokens.begin(), this->tokens.end(), option)
           != this->tokens.end();
}
