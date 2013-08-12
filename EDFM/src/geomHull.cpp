 /*!
 *	@file geomTetra.cpp
 *	@brief Class for Tetra in 3D space (definition). 
 *
 *	
 */

#include<cmath>
#include<limits>
#include<fstream>
#include<iomanip>
#include<vector>
#include<algorithm>
 
#include "geomTetra.hpp"
#include "geomHull.hpp"

namespace Geometry
{
	
	// --------------------   Class Hull  --------------------

// ==================================================
// Constructors & Destructor
// ==================================================
	
		
	//! Empty constructor
	Hull::Hull() :  M_simplexes(), M_points(), M_tetra(), M_Ntetra(0)  {}
//------------------costruttore 1: dati i punti e la matrice di connettività costruisce i tetraedri----------	

	//presuppone una chiamata a qhull

	Hull::Hull(const gmm::dense_matrix<gmm::size_type> simplexes, const std::vector<coordT> coordinates) 
	{  	gmm::copy(simplexes, M_simplexes);
		gmm::size_type npunti(coordinates.size()/3);
		for (gmm::size_type i=0; i<npunti;++i)
		{
			M_points[i].x=coordinates[i*3];
			M_points[i].y=coordinates[i*3+1];
			M_points[i].z=coordinates[i*3+2];
		}
		for (gmm::size_type i=0; i<simplexes.ncols();++i){
			Tetra elemento(M_points[simplexes(0,i)],M_points[simplexes(1,i)],M_points[simplexes(2,i)],M_points[simplexes(3,i)]);
			M_tetra[i]=elemento;
		}
		M_Ntetra=M_tetra.size();
	} //fine primo costruttore
	
	
//------------------costruttore 2: dati i punti costruisce i tetraedri   ----------	

	//chiama qhull al suo interno

	Hull::Hull(std::vector<coordT> coordinates) 
	{  
		gmm::size_type npunti(coordinates.size()/3);
		M_points.resize(npunti);
			
		for (gmm::size_type i=0; i<npunti;++i)
		{	M_points[i].x=coordinates[i*3];
			M_points[i].y=coordinates[i*3+1];
			M_points[i].z=coordinates[i*3+2];
		}

		bool exitcode=true;
 	       if (npunti <= 3) { gmm::resize(M_simplexes, 4, 0);  }
    	       if (npunti == 4) {
      			gmm::resize(M_simplexes, 4, 1);
      			for (gmm::size_type i=0; i <= 3; ++i) M_simplexes(i, 0) = i;
			
			exitcode=false;
    		}
	       
	       if (npunti>4){
	 	       exitcode=call_qhull(npunti, M_simplexes, coordinates);
		}
    	       if (!exitcode){
	       for (gmm::size_type i=0; i<M_simplexes.ncols();++i){
			Tetra elemento(M_points[M_simplexes(0,i)],M_points[M_simplexes(1,i)],M_points[M_simplexes(2,i)],M_points[M_simplexes(3,i)]);
			M_tetra.push_back(elemento);
		}
		M_Ntetra=M_tetra.size();
		}
	}  //fine secondo costruttore


//------------------costruttore 3: prende un vettore di point3d e poi lo trasforma, poi è = al 2 ----------	

	Hull::Hull(std::vector<Point3D> punti){
		gmm::size_type npunti(punti.size());
		std::vector<coordT> Pts(3*npunti);	
		M_points.resize(npunti);

		for (gmm::size_type i=0; i<npunti;++i){

			Pts[i*3]=punti[i].x;
			Pts[i*3+1]=punti[i].y;
			Pts[i*3+2]=punti[i].z;	
			M_points[i]=punti[i];
	
		}

		bool exitcode(true);
 	       if (npunti <= 3) { gmm::resize(M_simplexes, 4, 0);  }
    	       if (npunti == 4) {
      			gmm::resize(M_simplexes, 4, 1);
      			for (gmm::size_type i=0; i <= 3; ++i) M_simplexes(i, 0) = i;
			exitcode=false;
    			}
		if (npunti>4){
	       exitcode=call_qhull(npunti, M_simplexes, Pts);
}

    		if (!exitcode){
	       for (gmm::size_type i=0; i<M_simplexes.ncols();++i){
			Tetra elemento(M_points[M_simplexes(0,i)],M_points[M_simplexes(1,i)],M_points[M_simplexes(2,i)],M_points[M_simplexes(3,i)]);
			M_tetra.push_back(elemento);
		}
		
		
M_Ntetra=M_tetra.size();
		}
else{
M_Ntetra=0;

}

	}


	
//------------------distruttore-------------------------	

	Hull::~Hull() {}

//-----------------------wrapper per la chiamata a qhull---------------------------------------

	bool Hull::call_qhull(gmm::size_type & npunti, gmm::dense_matrix<gmm::size_type> & simplexes, std::vector<coordT> & coordinates){

	FILE *outfile= 0;    /* output from qh_produce_output()
                          *  use NULL to skip qh_produce_output() */ 
        FILE *errfile= stderr;    /* error messages from qhull code */ 
    	int exitcode;             /* 0 if no error from qhull */
 
	facetT *facet;	          /* set by FORALLfacets */
  	vertexT *vertex, **vertexp;

	char flags[]= "qhull QJ Q10 d Qbb Pp T0";
        
	exitcode = qh_new_qhull (3, npunti, &coordinates[0], 0, flags, outfile, errfile);
		
	gmm::size_type nbf;
 	if (!exitcode) { /* if no error */ 
  	     	nbf=0;
 	     	FORALLfacets { if (!facet->upperdelaunay) nbf++; }
      		gmm::resize(simplexes, 4, nbf);
       		nbf=0;
      		FORALLfacets {

        		if (!facet->upperdelaunay) {
	 			gmm:: size_type s=0;
			        FOREACHvertex_(facet->vertices) {
					    assert(s < (unsigned)(4));
            				    simplexes(s++,nbf) = qh_pointid(vertex->point);
				}
				nbf++;
			}
		}
	}
	return exitcode;
	}  //fine chiamata a qhull


	Real Hull::getVolume() const{
		Real volume(0);
		for (unsigned int ntr=0; ntr<M_simplexes.ncols();++ntr){
			Tetra tetraedrino(this->getTetra()[ntr]);
			volume=volume + tetraedrino.volume();
		}
	return volume;
	}

	std::vector<gmm::size_type> Hull::getPointsSimplex(gmm::size_type quale)
	const{
		std::vector<gmm::size_type> punti_tetra;
		for (gmm::size_type i=0; i<4;++i){
			punti_tetra.push_back(this->M_simplexes(i,quale));
		}
		return punti_tetra;
	}


} // namespace Geometry
