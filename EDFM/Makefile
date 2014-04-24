CXX = g++
WARFLAGS = -Wall -pedantic
# Optimised version
OPTFLAGS = -O3 -DNDEBUG
# For debugging
#OPTFLAGS = -g -fno-inline
CPPFLAGS=
# Use verbose for extra info
#CPPFLAGS=-DVERBOSE
# Activate OpenMP
OMPFLAGS = -fopenmp
# When debugging take away openmp
#OMPFLAGS =
CXXFLAGS = $(OPTFLAGS) $(WARFLAGS)  $(OMPFLAGS) $(CPPFLAGS) -I./src -I../external/include -I${mkQhullHome}/include
# Release version
#LDFLAGS= $(OPTFLAGS)
# For development version we hard link the loading of shared libs
LDFLAGS= $(OPTFLAGS) -Wl,-rpath=$(PWD)/lib -Wl,-rpath=${mkQhullHome}/lib
LDLIBS=-L${mkQhullHome}/lib -lqhull $(OMPFLAGS)
LOG_LINKER = link.log
LOG_COMPILER = compile.log
SRCS = $(wildcard ./src/*.cpp)
OBJ = $(SRCS:.cpp=.o)
HEADERS = $(wildcard ./src/*.hpp)

LIBDIR=$(PWD)/lib
INSTALLDIR = $(PWD)/lib

.PHONY : install libdev tests doc uninstall clean clean_lib clean_doc all old_properties doc

all:
	@echo
	@echo To install the library use:
	@echo make install
	@echo  "to remove: >make uninstall)"
	@echo  "to remove all object files: make clean"
	@echo  "to remove library flies files: make clean_lib"
	@echo  "to remove doc files: make clean_doc"
	@echo  "to remove all object, libs and doc: make dist_clean"
	@echo
	@echo  "to compile the program:"
	@echo  "make main_edfm"
	@echo  "to create the documentation"
	@echo  "make doc"
	@echo  "to compile the program with the old procedure:"
	@echo  "make old_properties"
doc:
	doxygen Doxyfile

install: $(INSTALLDIR)/libCPgeom.so.1.0

old_properties:
	$(MAKE) main_edfm_old

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

main_edfm_old: $(INSTALLDIR)/libCPgeom.so.1.0 main_edfm_old.o
	$(CXX) $(LDFLAGS) -o main_edfm_old main_edfm_old.o  $(LDLIBS) -L$(LIBDIR) -lCPgeom

main_edfm_old.o: main_edfm.cpp
	$(CXX) $(CXXFLAGS) -DOLD_PROPERTIES -c -o main_edfm_old.o main_edfm.cpp

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
	@rm -f main_edfm

clean_lib:
	@rm -f lib/libCPgeom.*

dist_clean: clean clean_lib clean_doc

uninstall:
	rm -f $(INSTALLDIR)libCPgeom.so
	rm -f $(INSTALLDIR)libCPgeom.so.1
	rm -f $(INSTALLDIR)libCPgeom.so.1.0
