#include "./src/CPgeom.hpp"
#include "./src/chrono.hpp"
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

	std::string outpath=cl("outpath","");
	bool noOutput=(outpath == std::string(""));
	if(noOutput) std::cout<<" No output directory specified. No output will be produced"<<std::endl;
	
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
	Real conv_z=cl("conv_z",1);
	
	bool rotate_z(false);

	if (cl("rotate_z","no")=="yes") {rotate_z=true;} 

	std::string direzione=cl("direzione", "x");
	Real angle(0.);
	std::string direzioneF=cl("direzioneF", "x");
	Real angleF(0.);
	if (direzione=="x") {angle=cl("theta_x",0);}
	if (direzione=="y") {angle=cl("theta_y",0);}
	if (direzioneF=="x") {angleF=cl("theta_xF",0);}
	if (direzioneF=="y") {angleF=cl("theta_yF",0);}

	
	CPgrid grid(gridFile,0,direzione, angle, rotate_z);
        
	std::string permFile=cl("PermFile","");
	PermeabilityField permfield(grid);
        
	if(permFile==std::string(""))
	  {
	    std::vector<Real>  uni(grid.Nx()*grid.Ny()*grid.Nz(),1.0);
	    std::cerr<<"No perm file file "<<fractureFile<<std::endl;
	    permfield.setPermx(uni);
	    permfield.setPermy(uni);
	    permfield.setPermz(uni);
	  }
	else  {
		permfield.importFromFile (permFile);
	}
	
	if (cl("exportGrid","no")=="yes"){	
	std::stringstream ss (std::stringstream::in | std::stringstream::out);
	ss<<outpath<<"/grid.vtk";
	grid.exportVtk(ss.str());
	std::stringstream ss_perm  (std::stringstream::in | std::stringstream::out);
	ss_perm<< outpath<<"/permX"<<".vtk";
	grid.exportVtk(ss_perm.str(), permfield.getPermxVec(),"PermX",0);
        ss_perm.str("");
        ss_perm<< outpath<<"/permY"<<".vtk";
	grid.exportVtk(ss_perm.str(), permfield.getPermyVec(),"PermY",0);
	ss_perm.str("");
        ss_perm<< outpath<<"/permZ"<<".vtk";
	grid.exportVtk(ss_perm.str(), permfield.getPermzVec(),"PermZ",0);
	}

	Fractures lista(fractureFile,conv_z, direzioneF, angleF);
	std::cout << "List of fractures has been created"<<std::endl;
		
	std::vector<GridIntersections> intMegaStore;	
	std::vector<CProp> propStore;	

	GridStrategy strategy(FOR3);
	if (cl("isOPT","yes")=="yes"){	
		strategy=FOR3OPT;
	}
	// ciclo sulle fratture
	Timings::Chrono calcoloint, calcologeom, esportazione;
	calcoloint.start();

	for (gmm::size_type  i=0; i<lista.M_nfractures;++i){
		
		if ((i+1)%10==0) {std::cout << "Processing fracture  "<<i+1<<"  of  "<<lista.M_nfractures<<"\n";}		
	
		// (3) Fault definition
		Fracture ff(lista.M_fractures[i]);
		Fault f(ff);
		// (4) To store computed intersection
		GridIntersections intersectionStore(grid.Nx(),grid.Ny(),grid.Nz());

		if (ff.areaFault()>0){
	
		// (5) Solver definition
		IntersectionSolver newton_FOR3(NEWTON,strategy);
	
		// (6) Set solver properties
		newton_FOR3.setMaxIteration(100);
	
		// (7) Solve intersection problem
		newton_FOR3(f,grid,intersectionStore);
		}
		intMegaStore.push_back(intersectionStore);

		
	}
	calcoloint.stop();
	
	std::cout << "All intersections computed in "<<calcoloint <<std::endl;

	lista.computeIntersections(0);

	calcologeom.start();
	for (gmm::size_type i=0; i<lista.M_nfractures;++i){

		Fracture ff(lista.M_fractures[i]);
		Fault f(ff);
		CProp propfaglia(intMegaStore[i], &grid, &ff, &permfield);
#ifdef OLD_PROPERTIES
		propfaglia.setProperties_old();
#else
		propfaglia.setProperties();
#endif
		propStore.push_back(propfaglia);
	}

//ricerca delle intersezioni  tra fratture

	for (gmm::size_type i=0;i<lista.M_nfractures;++i){
	(lista.M_fractures[i]).setGeoProp(propStore[i]);
	
	}

	lista.computeIntersections(1);

	calcologeom.stop();
	std::cout << "All fracture properties computed in "<<calcologeom <<std::endl;

	if(!noOutput)
	  {
	    esportazione.start();
	    for (gmm::size_type i=0; i<lista.M_nfractures;++i){
	      
	      Fracture ff(lista.M_fractures[i]);
	      Fault f(ff);
	      
	      if (cl("exportIntersections","no")=="yes"){	
		std::stringstream ss4 (std::stringstream::in | std::stringstream::out);
		ss4<< outpath<<"/intersections"<<i+1<<".vtk";
		intMegaStore[i].exportVtk(ss4.str());
	      }
	      
	      // (8) Export solutions
	      
	      if (cl("exportFractures","no")=="yes"){		
		std::stringstream ss (std::stringstream::in | std::stringstream::out);
		ss<< outpath<<"/fault"<<i+1<<".vtk";
		f.exportVtk(ss.str());
	      }
	      if (cl("exportGeoprop","no")=="yes"){	

		std::stringstream ss1 (std::stringstream::in | std::stringstream::out);
		ss1<< outpath<< "/grid_Aree"<<i+1<<".vtk";
		grid.exportVtk(ss1.str(), propStore[i].getAreas(),"area",0);
		
		std::stringstream ss2 (std::stringstream::in | std::stringstream::out);
		ss2<< outpath<<"/grid_vol"<<i+1<<".vtk";
		grid.exportVtk(ss2.str(), propStore[i].getVolumes(),"volumi",0);
		
		std::stringstream ss3 (std::stringstream::in | std::stringstream::out);
		ss3<< outpath<<"/grid_d"<<i+1<<".vtk";
		grid.exportVtk(ss3.str(), propStore[i].getDmedio(),"d",0);
	      }
	      if (cl("exportFrGrids","no")=="yes"){		
		std::stringstream ss4 (std::stringstream::in | std::stringstream::out);
	        ss4<<outpath<< "/meshF"<<i+1<<".vtk";
		(lista.M_fractures[i]).exportVtk(ss4.str());
	      }
	      
	    }
	    
	    for (gmm::size_type i=0;i<lista.M_nfractures;++i){
	      std::stringstream ss5 (std::stringstream::in | std::stringstream::out);
	      ss5<< outpath <<"/frattura"<<i+1<<".txt";
	      std::ofstream myfile;
	      myfile.open(ss5.str().c_str());
	      if(myfile.fail())
		{
		  std::cerr<<"Error opening file "<<ss5.str()<<std::endl;
		}
	      else
		{
		  (lista.M_fractures[i]).exportFracture2(myfile,i);
		}
	      myfile.close();
	    }
	    
	    std::stringstream ss6 (std::stringstream::in | std::stringstream::out);
	    ss6<< outpath <<"/intMatrix.txt";
	    std::ofstream myfile;
	    myfile.open(ss6.str().c_str());
	      if(myfile.fail())
		{
		  std::cerr<<"Error opening file "<<ss6.str()<<std::endl;
		}
	      else
		{
		  lista.exportIntersectionMatrix(myfile);
		}
	      myfile.close();
	    esportazione.stop();
	    std::cout << "Export in  "<<esportazione <<std::endl;
	  }
	
	return 0;
	
}
