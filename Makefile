CXXFLAGS = -Ofast -std=c++14 -Wall 
#CXXFLAGS = -g -std=c++14 -Wall

all: sparse

sparse: hash.cpp lib/sparsepp.h
	$(CXX) $(CXXFLAGS) -D_MAP="sparse_hash_map<string, int>" hash.cpp -o $@

std: hash.cpp
	$(CXX) $(CXXFLAGS) -D_MAP="map<string, int>" $^ -o $@

exectime: hash.cpp test.sh
	bash ./test.sh exec_time

maxmemory: hash.cpp test.sh
	bash ./test.sh max_memory

clean:
	rm std
	rm sparse
	rm graph.html

.PHONY: all clean exectime maxmemory
