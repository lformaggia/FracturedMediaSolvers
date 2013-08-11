CXX = g++
WARFLAGS = -Wall
OPTFLAGS = -O3 -DNDEBUG
#OPTFLAGS = -g
#OPTFLAGS = -g 
CPPFLAGS=
# Use verbose for extra info
#CPPFLAGS=-DVERBOSE
OMPFLAGS = -fopenmp
#OMPFLAGS = 
CXXFLAGS = $(OPTFLAGS) $(WARFLAGS)  $(OMPFLAGS) $(CPPFLAGS) -I./../src -I../../external/include
LDFLAGS= $(OPTFLAGS) -Wl,-rpath=$(PWD)/../lib -Wl,-rpath=$(PWD)/../../external/lib 
LDLIBS=-L../external/lib -lqhull $(OMPFLAGS)
LOG_LINKER = link.log
LOG_COMPILER = compile.log
SRCS = $(wildcard ./../src/*.cpp)
OBJ = $(SRCS:.cpp=.o)
HEADERS = $(wildcard ./../src/*.hpp)

LIBDIR=$(PWD)/../lib
INSTALLDIR = $(PWD)/../lib

.PHONY = test run clean

all:
	@echo
	@echo To compile the test
	@echo make main_test
	@echo To run the test
	@echo main run
	@echo To clean up
	@echo make clean 
	@echo  "(to remove: >make clean & >make clean_lib)"
	@echo

main_test: $(INSTALLDIR)/libCPgeom.so.1.0 main_test.o  
	$(CXX) $(LDFLAGS) -o main_test main_test.o  $(LDLIBS) -L$(LIBDIR) -lCPgeom

main_test.o: main_test.cpp
	$(CXX) $(CXXFLAGS) -c -o main_test.o main_test.cpp


$(OBJ): %.o: %.cpp $(HEADERS)
	$(CXX) -c $(CXXFLAGS) -fPIC $< -o $@

chrono.o: chrono.cpp chrono.hpp
	$(CXX) -c $(CXXFLAGS) -o chrono.o chrono.cpp

clean_doc:
	@rm -fr ./doc/

clean:
	@rm -f *.o
	@rm -f main_test

run:
	./runtest.sh
