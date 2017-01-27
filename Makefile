CXXFLAGS = -Ofast -std=c++14 -Wall
#CXXFLAGS = -g -std=c++14 -Wall

all: hash

hash: hash.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

clean:
	rm hash

.PHONY: all hash
