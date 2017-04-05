CXX = g++
CXXFLAGS = -Ofast -std=c++14 -Wall
CXXADDFLAGS = -MP -MD

SRCS=libtest.cpp optionsParser.cpp
OBJS=$(subst .cpp,.o,$(SRCS))
DEPS=$(subst .cpp,.d,$(SRCS))

all: sparse std

sparse: $(OBJS)
	$(CXX) $(CXXFLAGS) $(CXXADDFLAGS) -o $@ $(OBJS)

std: libtest.cpp optionsParser.cpp
	$(CXX) $(CXXFLAGS) -D __USE_STD_UNORDERED_MAP__ -c -o libtest.o libtest.cpp
	$(CXX) $(CXXFLAGS) -c -o optionsParser.o optionsParser.cpp
	$(CXX) -o $@ optionsParser.o libtest.o
	rm -f $(OBJS)

exectime: test.sh
	./test.sh exec_time

maxmemory: test.sh
	./test.sh max_memory

clean:
	rm -f $(OBJS) $(DEPS)
	rm -f ./std
	rm -f ./sparse

realclean: clean
	rm -f ./graph.html


.PHONY: all clean realclean exectime maxmemory sparse std
