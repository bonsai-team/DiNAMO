#include "optionsParser.hpp"


InputParser::InputParser(int &argc, char **argv){
    for (int i=1; i < argc; ++i)
        this->tokens.emplace_back(string(argv[i]));
}

const string &InputParser::getCmdOption(const std::initializer_list<const string> &options) const{
    vector<string>::const_iterator itr;
    for (auto const &option : options) {
        itr = find(this->tokens.begin(), this->tokens.end(), option);
        if (itr != this->tokens.end() && ++itr != this->tokens.end()){
            return *itr;
        }
    }
    return empty_string;
}

bool InputParser::cmdOptionExists(const std::initializer_list<const string> &options) const{
    for (auto const &option : options) {
        if (find(this->tokens.begin(), this->tokens.end(), option) != this->tokens.end())
            return true;
    }
    return false;
}
