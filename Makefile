CXXFLAGS = -Ofast -std=c++14 -Wall 
#CXXFLAGS = -g -std=c++14 -Wall

all: sparse

sparse: hash.cpp
	$(CXX) $(CXXFLAGS) -D_MAP="sparse_hash_map<string, int>" $^ -o $@

std: hash.cpp
	$(CXX) $(CXXFLAGS) -D_MAP="map<string, int>" $^ -o $@

clean:
	rm hash

.PHONY: all sparse std clean
