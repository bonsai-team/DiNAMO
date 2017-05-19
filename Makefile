#CXX = g++
CXXFLAGS += -std=c++14 -Wall -Ofast

SRCS=main.cpp hash.cpp optionsParser.cpp degenerate.cpp node.cpp mutual_information.cpp fisher_test.cpp graph_simplification.cpp reverse_complement.cpp meme_format.cpp find_redundant_motif.cpp
OBJS=$(subst .cpp,.o,$(SRCS))
DEPS=$(subst .cpp,.d,$(SRCS))

all: dinamo

dinamo: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

clean:
	rm -f $(OBJS) $(DEPS)

realclean: clean
	rm -f ./dinamo

.PHONY: all clean realclean artifact
