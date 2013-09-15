 /*!
 *	@file DarcyMass.hpp
 *	@brief Class for building a Mass-matrix for finite volume discretization.
 *
 *	@author Francesco Della Porta 
 *
 */ 

#ifndef __DARCYMASS_HPP__
#define __DARCYMASS_HPP__

#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>
#include "DarcyTypeDefinitions.hpp"
#include "MatrixHandler.hpp"
	
namespace Darcy{

/*!
	@class MassMatrix
	This class constructs the mass-matrix. The adopted technique is the one of the finite volume method: it hence represents the volume of the cell. Also fractures represented by Facets are considered as cells.   
*/
template <class T> 
class MassMatrix: public MatrixHandler<T>{

	typedef typename T::Generic_Point Generic_Point;
	typedef typename T::Generic_Vector Generic_Vector;
	typedef typename T::Fracture_Juncture Fracture_Juncture;

public:
	//! @name Constructor & Destructor
	//@{

	//! Constructor for a Mass-Matrix, given a Geometry::Rigid_Mesh.
	/*!
		@param rigid_mesh A Geometry::Rigid_Mesh on which the matrix is constructed
	*/
	MassMatrix(const Geometry::Rigid_Mesh<T> &rigid_mesh);
	//! Copy-Constructor deleted
	MassMatrix(const MassMatrix&) = delete;
	//! Empty-Constructor deleted
	MassMatrix() = delete;
	//! Default destructor
	~MassMatrix() = default;
	//@}

	//! @name Methods
	//@{
		
	//! assemble
	/*!
	 * @return Constructs the Mass matrix
	 */
	void assemble();
	//@}
};


// --------------------   Class MassMatrix   --------------------
 
// ==================================================
// Constructors & Destructor
// ==================================================

template<class T>
MassMatrix<T>::MassMatrix (const Geometry::Rigid_Mesh<T> &rigid_mesh): MatrixHandler<T>(rigid_mesh,D_Cell)
{}

 
// ==================================================
// Methods
// ==================================================

template<class T>
void MassMatrix<T>::assemble()
{
	std::vector<Triplet> Matrix_elements;
	for (auto cell_it : this->M_mesh.getCellsVector())
	{		 
		Matrix_elements.emplace_back(Triplet (cell_it.getId(), cell_it.getId(), cell_it.getVolume()));
	}	

	for (auto facet_it : this->M_mesh.getFractureFacetsIdsVector())
	{
		D_Real _volume = facet_it.Aperture()*facet_it.getFacet().size();
		Matrix_elements.emplace_back(Triplet (facet_it.getIdasCell(), facet_it.getIdasCell(), _volume));
		Matrix_elements.emplace_back(Triplet (facet_it.getSeparated()[0], facet_it.getSeparated()[0], -_volume/2));
		Matrix_elements.emplace_back(Triplet (facet_it.getSeparated()[1], facet_it.getSeparated()[1], -_volume/2));
	}
	this->M_Matrix->setFromTriplets(Matrix_elements.begin(), Matrix_elements.end());
}



}

#endif
