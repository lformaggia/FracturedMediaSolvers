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
#include "geomHull2D.hpp"

namespace Geometry
{
	
	// --------------------   Class Hull  --------------------

// ==================================================
// Constructors & Destructor
// ==================================================
	
		
	//! Empty constructor
	Hull2D::Hull2D() :  M_simplexes(), M_points(), M_triangle(), M_Ntriangle(0)  {}

	
//------------------costruttore 2: dati i punti costruisce i triangoli   ----------	

	//chiama qhull al suo interno

	Hull2D::Hull2D(std::vector<Point3D> punti, Fault* faglia) 
	{  

		gmm::size_type npunti(punti.size());
		M_points.resize(npunti);

		Point3D p1=punti[1]-punti[0];
		Point3D p2=punti[2]-punti[0];

		Point3D no(p1.cross(p2));

		std::vector<coordT> coordinates(npunti*2,0.0);

		for (gmm::size_type i=0; i<npunti;++i)
		{	M_points[i].x=punti[i].x;
			M_points[i].y=punti[i].y;
			M_points[i].z=punti[i].z;


			std::vector<Real> uv(faglia->inv_param(punti[i], no));
		
			coordinates[2*i]=uv[0];
			coordinates[2*i+1]=uv[1];
		}

	       bool exitcode=true;
 	       if (npunti <= 2) { gmm::resize(M_simplexes, 3, 0);  }
    	       if (npunti == 3) {
      			gmm::resize(M_simplexes, 3, 1);
      			for (gmm::size_type i=0; i <= 2; ++i) M_simplexes(i, 0) = i;
			exitcode=false;
    		}
	       
	       if (npunti>3){
	 	       exitcode=call_qhull2D(npunti, M_simplexes, coordinates);
		}
    	       if (!exitcode){
	       for (gmm::size_type i=0; i<M_simplexes.ncols();++i){
			Triangle elemento(M_points[M_simplexes(0,i)],M_points[M_simplexes(1,i)],M_points[M_simplexes(2,i)]);
			M_triangle.push_back(elemento);
		}
		
		
M_Ntriangle=M_triangle.size();
	}
else{
M_Ntriangle=0;

}
	}  //fine secondo costruttore


//------------------distruttore-------------------------	

	Hull2D::~Hull2D() {}

//-----------------------wrapper per la chiamata a qhull---------------------------------------

	bool Hull2D::call_qhull2D(gmm::size_type & npunti, gmm::dense_matrix<gmm::size_type> & simplexes, std::vector<coordT> & coordinates){

	FILE *outfile= 0;    /* output from qh_produce_output()
                          *  use NULL to skip qh_produce_output() */ 
        FILE *errfile= stderr;    /* error messages from qhull code */ 
    	int exitcode;             /* 0 if no error from qhull */
 
	facetT *facet;	          /* set by FORALLfacets */
  	vertexT *vertex, **vertexp;

	char flags[]= "qhull Qt d Qbb Pp T0 ";
        
	exitcode = qh_new_qhull (2, npunti, &coordinates[0], 0, flags, outfile, errfile);
		
	gmm::size_type nbf;
 	if (!exitcode) { /* if no error */ 
  	     	nbf=0;
 	     	FORALLfacets { if (!facet->upperdelaunay) nbf++; }
      		gmm::resize(simplexes, 3, nbf);
       		nbf=0;
      		FORALLfacets {

        		if (!facet->upperdelaunay) {
	 			gmm:: size_type s=0;
			        FOREACHvertex_(facet->vertices) {
					    assert(s < (unsigned)(3));
            				    simplexes(s++,nbf) = qh_pointid(vertex->point);
				}
				nbf++;
			}
		}
	}
	return exitcode;
	}  //fine chiamata a qhull


	Real Hull2D::getArea(){
		Real area(0);
		for (unsigned int ntr=0; ntr<M_simplexes.ncols();++ntr){
			Triangle triangolino(this->getTriangle()[ntr]);
			area=area + triangolino.area();
		}
	return area;
	}

	std::vector<gmm::size_type> Hull2D::getPointsSimplex(gmm::size_type quale){
		std::vector<gmm::size_type> punti_tetra;
		for (gmm::size_type i=0; i<3;++i){
			punti_tetra.push_back(this->M_simplexes(i,quale));
		}
		return punti_tetra;
	}


} // namespace Geometry
