CXX = g++
WARFLAGS = -Wall
OPTFLAGS = -O3
OMPFLAGS = -fopenmp
FLAGS = $(OPTFLAGS) $(WARFLAGS)  $(OMPFLAGS)
LOG_LINKER = link.log
LOG_COMPILER = compile.log
SRCS = $(wildcard ./src/*.cpp)
OBJ = $(SRCS:.cpp=.o)
HEADERS = $(wildcard ./src/*.hpp)

INSTALLDIR = /usr/lib/

.PHONY = install libdev tests doc uninstall clean clean_lib clean_doc

all:
	@echo
	@echo To install the library use:
	@echo >sudo
	@echo >make install
	@echo  (to remove: >make uninstall)
	@echo To install the dev version in the current directory use:
	@echo >make libdev
	@echo  (to remove: >make clean & >make clean_lib)
	@echo

doc:
	doxygen Doxyfile

install: $(OBJ)
	$(CXX) $(FLAGS) -shared -Wl,-soname,libCPgeom.so.1 -o libCPgeom.so.1.0 $(OBJ)
	@mv libCPgeom.so.1.0 $(INSTALLDIR)libCPgeom.so.1.0
	ln -s $(INSTALLDIR)libCPgeom.so.1.0 $(INSTALLDIR)libCPgeom.so.1	# soname
	ln -s $(INSTALLDIR)libCPgeom.so.1.0 $(INSTALLDIR)libCPgeom.so	# link name
	$(MAKE) clean
	$(MAKE) tests

libdev: $(OBJ) 
	$(CXX) $(FLAGS) -shared -Wl,-soname,libCPgeom.so.1 -o libCPgeom.so.1.0 $(OBJ)  -rdynamic  /usr/lib/libqhull.so.5
	ln -s libCPgeom.so.1.0 libCPgeom.so.1	# soname
	ln -s libCPgeom.so.1.0 libCPgeom.so	# link name
	$(MAKE) clean

tests: main_anna main_tutorial 
#main_Bashkirian main_testRefined main_tutorial 

main_anna: main_anna.o chrono.o
	$(CXX) $(FLAGS) -o main_anna main_anna.o chrono.o   -rdynamic  /usr/lib/libqhull.so.5 -L. -lCPgeom 

main_anna.o: main_anna.cpp
	$(CXX) $(FLAGS) -c -o main_anna.o main_anna.cpp

main_tutorial: main_tutorial.o chrono.o
	$(CXX) $(FLAGS) -o main_tutorial main_tutorial.o chrono.o   -rdynamic  /usr/lib/libqhull.so.5 -L. -lCPgeom 

main_tutorial.o: main_tutorial.cpp
	$(CXX) $(FLAGS) -c -o main_tutorial.o main_tutorial.cpp

main_Bashkirian: main_Bashkirian.o chrono.o
	$(CXX) $(FLAGS) -o main_Bashkirian main_Bashkirian.o chrono.o -L. -lCPgeom

main_Bashkirian.o: main_Bashkirian.cpp
	$(CXX) $(FLAGS) -c -o main_Bashkirian.o main_Bashkirian.cpp

main_testRefined: main_testRefined.o chrono.o
	$(CXX) $(FLAGS) -o main_testRefined main_testRefined.o chrono.o -L. -lCPgeom

main_testRefined.o: main_testRefined.cpp
	$(CXX) $(FLAGS) -c -o main_testRefined.o main_testRefined.cpp

$(OBJ): %.o: %.cpp $(HEADERS)
	$(CXX) -c $(FLAGS) -fpic $< -o $@

chrono.o: chrono.cpp chrono.hpp
	$(CXX) -c $(FLAGS) -o chrono.o chrono.cpp

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
