CXX = g++
WARFLAGS = -Wall
#OPTFLAGS = -O3
OPTFLAGS = -g
OMPFLAGS = -fopenmp
CXXFLAGS = $(OPTFLAGS) $(WARFLAGS)  $(OMPFLAGS) -I./src -I../external/include
LDFLAGS= $(OPTFLAGS) -Wl,-rpath=$(PWD)/lib -Wl,-rpath=$(PWD)/../external/lib 
LDLIBS=-L../external/lib -lqhull $(OMPFLAGS)
LOG_LINKER = link.log
LOG_COMPILER = compile.log
SRCS = $(wildcard ./src/*.cpp)
OBJ = $(SRCS:.cpp=.o)
HEADERS = $(wildcard ./src/*.hpp)

LIBDIR=$(PWD)/lib
INSTALLDIR = $(PWD)/lib

.PHONY = install libdev tests doc uninstall clean clean_lib clean_doc

all:
	@echo
	@echo To install the library use:
	@echo >make install
	@echo  "(to remove: >make uninstall)"
	@echo To install the dev version in the current directory use:
	@echo >make libdev
	@echo  "(to remove: >make clean & >make clean_lib)"
	@echo

doc:
	doxygen Doxyfile

install: $(INSTALLDIR)/libCPgeom.so.1.0

$(INSTALLDIR)/libCPgeom.so.1.0: $(OBJ)
	$(CXX) $(LDFLAGS) -shared -Wl,-soname,libCPgeom.so.1 -o libCPgeom.so.1.0 $(OBJ) $(LDLIBS)
	$(AR) rs libCPgeom.a $(OBJ)
	@mv libCPgeom.so.1.0 $(INSTALLDIR)/libCPgeom.so.1.0
	@mv libCPgeom.a $(INSTALLDIR)/
	ln -s -f $(INSTALLDIR)/libCPgeom.so.1.0 $(INSTALLDIR)/libCPgeom.so.1	# soname
	ln -s -f $(INSTALLDIR)/libCPgeom.so.1.0 $(INSTALLDIR)/libCPgeom.so	# link name

libdev: install

#libdev: $(OBJ) 
#	$(CXX) $(CXXFLAGS) $(LDFLAGS) -shared -Wl,-soname,libCPgeom.so.1 -o libC#Pgeom.so.1.0 $(OBJ)  $(LDLIBS)
#	ln -s libCPgeom.so.1.0 libCPgeom.so.1	# soname
#	ln -s libCPgeom.so.1.0 libCPgeom.so	# link name

tests: main_anna main_tutorial 
#main_Bashkirian main_testRefined main_tutorial 

main_anna: main_anna.o
	$(CXX) $(LDFLAGS) -o main_anna main_anna.o  $(LDLIBS) -L -lCPgeom 

main_anna.o: main_anna.cpp
	$(CXX) $(CXXFLAGS) -c -o main_anna.o main_anna.cpp

main_tutorial: $(INSTALLDIR)/libCPgeom.so.1.0 main_tutorial.o  
	$(CXX) $(LDFLAGS) -o main_tutorial main_tutorial.o  $(LDLIBS) -L$(LIBDIR) -lCPgeom 

main_tutorial.o: main_tutorial.cpp
	$(CXX) $(CXXFLAGS) -c -o main_tutorial.o main_tutorial.cpp

main_Bashkirian: main_Bashkirian.o chrono.o
	$(CXX) $(LDFLAGS) -o main_Bashkirian main_Bashkirian.o $(LDLIBS) -lCPgeom

main_Bashkirian.o: main_Bashkirian.cpp
	$(CXX) $(CXXFLAGS) -c -o main_Bashkirian.o main_Bashkirian.cpp

main_testRefined: main_testRefined.o
	$(CXX) $(LDFLAGS) -o main_testRefined main_testRefined.o  $(LDLIBS) -lCPgeom

main_testRefined.o: main_testRefined.cpp
	$(CXX) $(CXXFLAGS) -c -o main_testRefined.o main_testRefined.cpp

$(OBJ): %.o: %.cpp $(HEADERS)
	$(CXX) -c $(CXXFLAGS) -fPIC $< -o $@

chrono.o: chrono.cpp chrono.hpp
	$(CXX) -c $(CXXFLAGS) -o chrono.o chrono.cpp

clean_doc:
	@rm -fr ./doc/

clean:
	@rm -f *.o
	@rm -f ./src/*.o
	@rm -f main_Bashkirian
	@rm -f main_tutorial
	@rm -f main_testRefined
	@rm -f main_anna

clean_lib:
	@rm -f libCPgeom.so*

uninstall:
	rm -f $(INSTALLDIR)libCPgeom.so
	rm -f $(INSTALLDIR)libCPgeom.so.1
	rm -f $(INSTALLDIR)libCPgeom.so.1.0
