CXX = g++
CXXFLAGS = -std=c++14 -Wall -Ofast
CXXADDFLAGS = -MP -MD

SRCS=main.cpp hash.cpp optionsParser.cpp degenerate.cpp node.cpp mutual_information.cpp fisher_test.cpp stat_container.cpp
OBJS=$(subst .cpp,.o,$(SRCS))
DEPS=$(subst .cpp,.d,$(SRCS))

all: artifact

artifact: $(OBJS)
	$(CXX) $(CXXFLAGS) $(CXXADDFLAGS) -o $@ $(OBJS)

clean:
	rm -f $(OBJS) $(DEPS)

realclean: clean
	rm -f ./artifact

.PHONY: all clean realclean artifact
