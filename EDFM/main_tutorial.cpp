#include "./src/CPgeom.hpp"
#include <GetPot.hpp>

  extern "C" {
#ifdef _MSC_VER
# include <libqhull/qhull_a.h>
#else
# include <libqhull/libqhull.h>
//# include <qhull/mem.h>
# include <libqhull/qset.h>
# include <libqhull/geom.h>
# include <libqhull/merge.h>
# include <libqhull/poly.h>
# include <libqhull/io.h>
# include <libqhull/stat.h>
#endif
}
#include "gmm/gmm.h"
#include <stdlib.h>
#include <sstream>

void printHelp(){
  using std::cout;
  using std::endl;
  cout<<"*** Line Options ***"<<endl;
  cout<<"[-h  --help] This help"<<endl;
  cout<<"[InputFile=string] Input file name (data.pot)"<<endl<<endl;
}

int main(int argc, char** argv){ 
	 
	// (0) load namespace
	using namespace Geometry;
	using namespace Intersect;
	// Process help
	GetPot key_input(argc,argv);
	if (key_input.search(2, "--help", "-h")){
	  printHelp();
	  exit(0);
	}
	std::string parameterFile=key_input("InputFile","data.pot");
	
	GetPot   cl(parameterFile);
	std::string gridFile=cl("GridFile","");
	if(gridFile==std::string(""))
	  {
	    std::cerr<<"Wrong grid file "<<gridFile<<std::endl;
	    std::exit(1);
	  }
	std::string fractureFile=cl("FractureFile","");
	if(fractureFile==std::string(""))
	  {
	    std::cerr<<"Wrong fracture file "<<fractureFile<<std::endl;
	    std::exit(1);
	  }
	// (1) Import grid
//	CPgrid grid("./data/EP_3D_grid.GRDECL",0);
	CPgrid grid(gridFile,0);
		
//	CPgrid grid("./data/gridTest",0);
	
	// (2) Points definition


	Fractures lista(fractureFile);
//	Fractures lista("lista_mia3.fab");

	std::vector<GridIntersections> intMegaStore;	
	std::vector<CProp> propStore;	

	// ciclo sulle fratture

	for (gmm::size_type  i=0; i<lista.M_nfractures;++i){
		
		// (3) Fault definition
		Fracture ff(lista.M_fractures[i]);
		Fault f(ff);

		// (4) To store computed intersection
		GridIntersections intersectionStore(grid.Nx(),grid.Ny(),grid.Nz());
	
		// (5) Solver definition
		IntersectionSolver newton_FOR3(NEWTON,FOR3);
	
		// (6) Set solver properties
		newton_FOR3.setMaxIteration(100);
	
		// (7) Solve intersection problem
		newton_FOR3(f,grid,intersectionStore);
		intMegaStore.push_back(intersectionStore);

		std::stringstream ss4 (std::stringstream::in | std::stringstream::out);
	        ss4<< "./data/Tutorial/intersections"<<i<<".vtk";
		intersectionStore.exportVtk(ss4.str());
	}

	lista.computeIntersections(0);
	

	for (gmm::size_type i=0; i<lista.M_nfractures;++i){

		Fracture ff(lista.M_fractures[i]);
		Fault f(ff);
		
		CProp propfaglia(intMegaStore[i], &grid, &ff);
		propfaglia.setProperties();
		propStore.push_back(propfaglia);

		// (8) Export solutions
		
		std::stringstream ss (std::stringstream::in | std::stringstream::out);
		ss<< "./data/Tutorial/fault"<<i<<".vtk";
		f.exportVtk(ss.str());
	
		std::stringstream ss1 (std::stringstream::in | std::stringstream::out);
	        ss1<< "./data/Tutorial/grid_Aree"<<i<<".vtk";
		grid.exportVtk(ss1.str(), propfaglia.getAreas(),"area",0);
	     
		std::stringstream ss2 (std::stringstream::in | std::stringstream::out);
	        ss2<< "./data/Tutorial/grid_vol"<<i<<".vtk";
		grid.exportVtk(ss2.str(), propfaglia.getVolumes(),"volumi",0);
		
		std::stringstream ss3 (std::stringstream::in | std::stringstream::out);
	        ss3<< "./data/Tutorial/grid_d"<<i<<".vtk";
		grid.exportVtk(ss3.str(), propfaglia.getDmedio(),"d",0);

	}

	grid.exportVtk("./data/Tutorial/grid.vtk");
	std::string nomefile("./data/Tutorial/fratture");
	std::ofstream myfile;
	myfile.open(nomefile.c_str());
	//ricerca delle intersezioni  tra fratture

	for (gmm::size_type i=0;i<lista.M_nfractures;++i){

	(lista.M_fractures[i]).setGeoProp(propStore[i]);
	std::stringstream ss4 (std::stringstream::in | std::stringstream::out);
	        ss4<< "./data/Tutorial/meshF"<<i<<".vtk";
		(lista.M_fractures[i]).exportVtk(ss4.str());
		std::cout << "frattura   "<<i<<"  esportata"<<std::endl;
	}

	lista.computeIntersections(1);

	for (gmm::size_type i=0;i<lista.M_nfractures;++i){
	(lista.M_fractures[i]).exportFracture(myfile,i);
	}
	std::cout << "Esportate le fratture"<<std::endl;
	return 0;
	
}
