CXX = g++
WARFLAGS = -Wall
OPTFLAGS = -O3 -DNDEBUG
#OPTFLAGS = -g -fno-inline
#OPTFLAGS = -g 
CPPFLAGS=
# Use verbose for extra info
#CPPFLAGS=-DVERBOSE
OMPFLAGS = -fopenmp
#OMPFLAGS = 
CXXFLAGS = $(OPTFLAGS) $(WARFLAGS)  $(OMPFLAGS) $(CPPFLAGS) -I./src -I../external/include
LDFLAGS= $(OPTFLAGS) -Wl,-rpath=$(PWD)/lib -Wl,-rpath=$(PWD)/../external/lib 
LDLIBS=-L../external/lib -lqhull $(OMPFLAGS)
LOG_LINKER = link.log
LOG_COMPILER = compile.log
SRCS = $(wildcard ./src/*.cpp)
OBJ = $(SRCS:.cpp=.o)
HEADERS = $(wildcard ./src/*.hpp)

LIBDIR=$(PWD)/lib
INSTALLDIR = $(PWD)/lib

.PHONY = install libdev tests doc uninstall clean clean_lib clean_doc all old_properties

all:
	@echo
	@echo To install the library use:
	@echo make install
	@echo  "(to remove: >make uninstall)"
	@echo  "(to remove all object files: >make clean & >make clean_lib)"
	@echo
	@echo  "to compile the program:"
	@echo  "make main_edfm"
	@echo  "to compile the program with the old procedule:"
	@echo  "make old_properties"
	@echo  " to create the documentation"
	@echo  " make doc"
doc:
	doxygen Doxyfile

install: $(INSTALLDIR)/libCPgeom.so.1.0

old_properties:
	$(MAKE) main_edfm -DOLD_PROPERTIES

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


main_edfm: $(INSTALLDIR)/libCPgeom.so.1.0 main_edfm.o  
	$(CXX) $(LDFLAGS) -o main_edfm main_edfm.o  $(LDLIBS) -L$(LIBDIR) -lCPgeom 

main_edfm.o: main_edfm.cpp
	$(CXX) $(CXXFLAGS) -c -o main_edfm.o main_edfm.cpp

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
