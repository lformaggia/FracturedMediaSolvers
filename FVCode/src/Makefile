##Francesco Della Porta
#makefile for pacs project

CXX = g++
WARFLAGS = -Wall
CGALFLAGS = -frounding-math
EIGENLIB = -I /home/LOCAL/include/eigen3
CGALLIB = -lCGAL -lgmp

SRCS = $(wildcard ./src_Turconi/*.cpp)
OBJ = $(SRCS:.cpp=.o)
DP_HEADERS = $(wildcard ./DP_src/*.hpp)
T_HEADERS = $(wildcard ./src_Turconi/*.hpp)
HEADERS = $(DP_HEADERS) $(T_HEADERS)

ROOT_FOLDER=$(PWD)/..
LIB_ROOT=$(ROOT_FOLDER)/lib/
INCLUDE_ROOT=$(ROOT_FOLDER)/headers/

MYHEADERS = -I $(PWD)

export DEBUG
ifeq ($(DEBUG),no)
  OPTFLAGS = -O3 -std=c++0x -funroll-loops -Wall
  MESH2DLIB = -L. -L$(LIB_ROOT) -lMesh2D
else
  OPTFLAGS = -g -O0 -std=c++0x -Wall
  MESH2DLIB = -L.  -L$(LIB_ROOT) -lMesh2D
endif

FLAGS = $(OPTFLAGS) $(WARFLAGS) $(CGALFLAGS)

.PHONY = doc clean clean_doc main

.DEFAULT_GOAL = all

all:
	@echo "make library makes static library for Mesh2D"
	@echo "make install moves Mesh2D library and headers"
	@echo "make Test compiles executable with a test."
	@echo "make Oil compiles executable with a second test."
	@echo "make clean cleans all not library not Doxyfiles"
	@echo "make clean_doc cleans Doxyfiles"
	@echo "make distclean cleans all"
	@echo "macro DEBUG=no deactivates debugging and uses libraries in the"
	@echo "root lib when compiling executable"

Test: mainTest

mainTest: mainTest.o
	$(CXX) $(FLAGS) -o mainTest mainTest.o $(CGALLIB) $(MESH2DLIB)

mainTest.o: maintest.cpp
	$(CXX) $(FLAGS) $(EIGENLIB) -c -o mainTest.o maintest.cpp

Oil: mainOil

mainOil: mainOil.o
	$(CXX) $(FLAGS) -o mainOil mainOil.o $(CGALLIB) $(MESH2DLIB)

mainOil.o: mainoil.cpp
	$(CXX) $(FLAGS) $(EIGENLIB) -c -o mainOil.o mainoil.cpp

$(OBJ): %.o: %.cpp $(HEADERS)
	$(CXX) -c $(FLAGS) $(MYHEADERS) $< -o $@

doc:
	doxygen Doxyfile

library: $(OBJ)
	ar -rv libMesh2D.a $(OBJ)
	ranlib libMesh2D.a

install:
	cp $(T_HEADERS) $(INCLUDE_ROOT)/src_Turconi/
	cp -r DP_Darcy $(INCLUDE_ROOT)
	cp -r DP_src $(INCLUDE_ROOT)	
	mv -f libMesh2D.a $(LIB_ROOT)

clean_doc:
	@rm -fr ./doc/

clean:
	@rm -f *.o
	@rm -f ./src_Turconi/*.o
	@rm -f mainTest mainOil

distclean:
	$(MAKE) clean
	-\rm -r -f ./doc *.a

