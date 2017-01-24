COMP=g++
OPT=-std=c++14 -Wall

all: hash

hash: hash.cpp
	$(COMP) $(OPT) $^ -o $@

clean:
	rm hash

.PHONY: all hash
