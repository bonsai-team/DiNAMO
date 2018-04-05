SRCDIR := src
BUILDDIR := build
EXECDIR := bin
TARGET := bin/dinamo
INC := -I include 'C:\msys64\mingw64\include' -L'C:\msys64\mingw64\lib'

SRCEXT := cpp
SOURCES := $(wildcard $(SRCDIR)/*.$(SRCEXT))
OBJECTS := $(patsubst %.$(SRCEXT),$(BUILDDIR)/%.o,$(notdir $(SOURCES)))

CXXFLAGS := --std=c++14 -Wall -Ofast

all: $(TARGET)

$(TARGET): $(OBJECTS)
	@echo " Linking..."
	@mkdir -p $(EXECDIR)
	@echo " $(CXX) $^ -o $(TARGET) "; $(CXX) $^ -o $(TARGET)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	@echo " $(CXX) $(CXXFLAGS) $(INC) -c -o $@ $<"; $(CXX) $(CXXFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning...";
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)

# Tests
# tester:
#	$(CCX) $(CFLAGS) test/tester.cpp $(INC) $(LIB) -o bin/tester

# Spikes
# ticket:
#	$(CCX) $(CFLAGS) spikes/ticket.cpp $(INC) $(LIB) -o bin/ticket

.PHONY: clean
