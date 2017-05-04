CXX = g++
CXXFLAGS = -std=c++14 -Wall -Ofast
CXXADDFLAGS = -MP -MD

SRCS=main.cpp hash.cpp optionsParser.cpp degenerate.cpp node.cpp mutual_information.cpp fisher_test.cpp graph_simplification.cpp
OBJS=$(subst .cpp,.o,$(SRCS))
DEPS=$(subst .cpp,.d,$(SRCS))

all: artifact

dinamo: $(OBJS)
	$(CXX) $(CXXFLAGS) $(CXXADDFLAGS) -o $@ $(OBJS)

clean:
	rm -f $(OBJS) $(DEPS)

realclean: clean
	rm -f ./dinamo

.PHONY: all clean realclean artifact
