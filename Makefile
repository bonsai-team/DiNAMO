CXX = g++
CXXFLAGS = -Ofast -std=c++14 -Wall
CXXADDFLAGS = -MP -MD

SRCS=hash.cpp optionsParser.cpp
OBJS=$(subst .cpp,.o,$(SRCS))
DEPS=$(subst .cpp,.d,$(SRCS))

all: sparse

sparse: $(OBJS)
	$(CXX) $(CXXFLAGS) $(CXXADDFLAGS) -o $@ $(OBJS)

std: hash.cpp optionsParser.cpp
	$(CXX) $(CXXFLAGS) -D __USE_STD_UNORDERED_MAP__ -c -o hash.o hash.cpp
	$(CXX) $(CXXFLAGS) -c -o optionsParser.o optionsParser.cpp
	$(CXX) -o $@ optionsParser.o hash.o
	rm -f $(OBJS)

exectime: test.sh
	./test.sh exec_time

maxmemory: test.sh
	./test.sh max_memory

clean:
	rm -f $(OBJS) $(DEPS)

realclean: clean
	rm -f ./std
	rm -f ./sparse
	rm -f ./graph.html

.PHONY: all clean realclean exectime maxmemory sparse std
