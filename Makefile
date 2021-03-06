# compiler (needs to support -std=c++20)
CXX = g++-10

# location of third party headers (Eigen, CGAL, boost)
THIRD_PARTY_INC = -I${mkEigenInc} -I${mkBoostInc} -I${mkCgalInc}

# compiler optimization flags
OPTCXXFLAGS = -Ofast -fopenmp

# preprocessor optional flags (set to -DVERBOSE for timings)
OPTCPPFLAGS = -DVERBOSE

# compiler flags
CXXFLAGS = -std=c++20 -fno-pic -no-pie $(OPTCXXFLAGS)

# preprocessor flags (header paths and dependencies flags)
CPPFLAGS = -I./utilities -I./mesh -I./core -I./functional -I./external \
	$(THIRD_PARTY_INC) \
	$(OPTCPPFLAGS) \
	-MP -MD \

# linker flags
LDFLAGS = $(CXXFLAGS)

# linker default rule
LINK.o = $(CXX) $(LDFLAGS) $(LDLIBS) $(TARGET_ARCH)

# names of all executables
EXE_ALL = hconv pconv test

# names of executables with path
EXE_ALL_WITH_PATH = $(EXE_ALL:%=./tests/%)

# source files with path for the executables
SRC_EXE = $(EXE_ALL_WITH_PATH:=.cpp)

# objects files with path for the executables
OBJ_EXE = $(SRC_EXE:.cpp=.p)

# directories where to search for library source and header files
DIRS = ./mesh ./core ./functional ./utilities
vpath %.cpp $(DIRS)
vpath %.hpp $(DIRS)

# library source files with path
SRC = $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.cpp))

# library source files without path
SRC_WITHOUT_PATH = $(notdir $(SRC))

# library object files with path (all inside the build directories along with .d files)
OBJ = $(SRC_WITHOUT_PATH:%.cpp=./build/%.o)

# names of the library header files with paths
INC = $(wildcard *.hpp) $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.hpp))

# rules that are always run no matter what
.PHONY: all clean distclean doc shared static library $(EXE_ALL)

# this rule is to print the help
help:
	@echo "type in the current directory one of these commands:\n"
	@echo "make all\n    build all executables directly linked with object files (inside the ./tests folder) and the documentation\n"
	@echo "make <name-of-executable>\n    build <name-of-executable> (inside the ./tests folder, where <name-of-executable> can be hconv | pconv | test)\n"
	@echo "make shared\n    build a shared library in the current directory (this first runs make clean to rebuild object files with -fPIC)\n"
	@echo "make static\n    build a static library in the current directory (this first runs make clean to rebuild object files without -fPIC)\n"
	@echo "headers will be in ./wave/include, libraries in ./wave/lib, the library global header is wave.hpp, the soname is libwave.so\n"
	@echo "to run an executable go in ./tests and type\n    ./hconv -f hconv.pot\n    ./pconv -f pconv.pot\n    ./test ./data/<meshfile>\nwhere <meshfile> is a file in ./tests/data and can be voro27.voro | voro125.voro | voro1000.voro | voro8000.voro\n"
	@echo "make doc\n    builds the documentation in html format with doxygen\n"
	@echo "make clean\n    removes all .o files and .d files\n"
	@echo "make distclean\n    restores the current folder as it was removing all constructed files and folders (.o, .d, exectuables, documentation, libraries)\n"
	@echo "make run\n    after building the test executable with make test, it can be run for quick testing"

# this rule will make all executables linked directly to object files
all: $(EXE_ALL) doc

# this rule is to run a quick test
run:
	./tests/test ./tests/data/voro1000.voro

# this rule is to be able to type make <executable> from the current directory while the true executable is inside the ./tests folder
$(EXE_ALL): %: ./tests/%

# this rule will build the executable from the library object files and its own .o file (the $ symbol needs to be expanded twice)
.SECONDEXPANSION:
$(EXE_ALL_WITH_PATH): $(OBJ) $$@.o

# this rule is for building object files (except object files for executables)
./build/%.o: %.cpp
	mkdir -p ./build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $^ -o $@

# this rule prepares the library directories and include
library: clean
	mkdir -p ./wave/lib
	mkdir -p ./wave/include
	cp $(INC) ./wave/include

# this rule is to build the shared library in the current directory
shared: library ./wave/lib/libwave.so

# this rule will add the position independent code flag to compiler flags
shared: CXXFLAGS=-std=c++20 -fPIC $(OPTCXXFLAGS)

./wave/lib/libwave.so: $(OBJ)
	$(CXX) -shared -Wl,-soname,libwave.so -o $@ $^ 

# this rule is to build the static library in the current directory
static: library ./wave/lib/libwave.a

./wave/lib/libwave.a: $(OBJ)
	$(AR) -rs $@ $^

# rule to build the documentation in html format
doc:
	mkdir -p ./doc
	doxygen Doxyfile

# this rule cleans all .o and .d files
clean:
	$(RM) ./build/*.o
	$(RM) ./build/*.d
	$(RM) ./tests/*.o
	$(RM) ./tests/*.d

# this rule restores the folder as it was
distclean: clean
	$(RM) -r -f ./doc
	$(RM) -r -f ./wave
	$(RM) -r -f ./build
	$(RM) $(EXE_ALL_WITH_PATH)
	$(RM) -r -f ./tests/*.vtk

# include the autogenerated rules for dependencies of sources from headers
-include $(SRC_ALL:.cpp=.d)