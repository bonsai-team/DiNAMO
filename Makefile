CXXFLAGS = -std=c++14 -Wall

all: hash

hash: hash.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

clean:
	rm hash

.PHONY: all hash
