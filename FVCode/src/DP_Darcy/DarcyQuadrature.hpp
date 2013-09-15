/*!
 *	@file DarcyQuadrature.hpp
 *	@brief Class for integrating functions on a Rigid_Mesh.
 *
 *	@author Francesco Della Porta 
 *
 */ 

#ifndef __DARCYQUADRATURE_HPP__
#define __DARCYQUADRATURE_HPP__

#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>
#include "DarcyTypeDefinitions.hpp"
#include "DarcyQuadratureRules.hpp"
	

namespace Darcy{

/*!
	@class Quadrature
	This class integrates solutions to PDE solved with the finite volume method, computes the L2Norm of such a solution, and integrates continous functions. By passing a QuadratureRule to the constructor it is possible to choose the quadrature rule you want to use. 
*/
template <class T> 
class Quadrature {

	typedef typename T::Generic_Point Generic_Point;
	typedef typename T::Generic_Vector Generic_Vector;
	typedef typename T::Fracture_Juncture Fracture_Juncture;
	typedef typename QuadratureRule<T>::QuadratureRuleHandler QR_Handler;

public:
	//! @name Constructor & Destructor
	//@{
	//! Constructor for Quadrature, given a Geometry::Rigid_Mesh.
	/*!
		@param rigid_mesh A Geometry::Rigid_Mesh on which we want to integrate
	*/
	Quadrature (const Geometry::Rigid_Mesh<T> &rigid_mesh);
	//! Constructor for Quadrature, given a Geometry::Rigid_Mesh and a QuadratureRule
	/*!
		@param rigid_mesh A Geometry::Rigid_Mesh on which we want to integrate
		@param quadrature A Darcy::QuadratureRule that we want to use in order to integrate std::function
	*/
	Quadrature (const Geometry::Rigid_Mesh<T> &rigid_mesh, const QuadratureRule<T>& quadrature);
	//! Constructor for Quadrature, given a Geometry::Rigid_Mesh, a QuadratureRule, and a QuadratureRule to integrate fractures
	/*!
		@param rigid_mesh A Geometry::Rigid_Mesh on which we want to integrate
		@param quadrature A Darcy::QuadratureRule that we want to use in order to integrate std::function
		@param fracturequadrature A Darcy::QuadratureRule that we want to use in order to integrate std::function on fractures
	*/
	Quadrature (const Geometry::Rigid_Mesh<T> &rigid_mesh, const QuadratureRule<T>& quadrature, const QuadratureRule<T>& fracturequadrature);
	//@}


	//! @name Methods
	//@{
		
	//! Integrate std::function
	/*!
		@param Integrand A std::function which returns a scalar
	 	@return The integral of the considered Integrand function
	 */
	D_Real Integrate (const std::function<D_Real(Generic_Point)>& Integrand); 
	//! Integrate std::function and return the integral cell by cell
	/*!
		@param Integrand A std::function which returns a scalar
	 	@return A vector with in the i-th component the integral of the considered Integrand function on the i-th cell 
	 */
	Vector CellIntegrate (const std::function<D_Real(Generic_Point)>& func); 
	//! Integrate discrete function
	/*!
		@param Integrand A vector with in the i-th component the value of the function on the i-th cell 
	 	@return The integral of the considered Integrand function
	 */
	D_Real Integrate (const Vector& Integrand);
	//! L2 Norm of a discrete function
	/*!
		@param Integrand A vector with in the i-th component the value of the function on the i-th cell 
	 	@return The L2 norm of the considered Integrand function
	 */
	D_Real L2Norm (const Vector& Integrand);
	//@}

	//! @name Get Methods
	//@{
		
	//! Get Mesh size (const)
	/*!
	 * @return A the number of the cells in the Rigid_Mesh considered
	 */
	UInt getMeshSize() const
		{return M_size;}
	//@}

protected:
	//! A reference to a Rigid_Mesh
	const Geometry::Rigid_Mesh<T> & M_mesh;
	//! The number of cells in M_mesh
	D_UInt M_size;
	//! A pointer to a QuadratureRule
	QR_Handler M_quadrature;
	//! A pointer to a QuadratureRule which is used to integrate on fractures
	QR_Handler M_fractureQuadrature;
};

// --------------------   Class Quadrature   --------------------
 
// ==================================================
// Constructors & Destructor
// ==================================================

template<class T>
Quadrature<T>::Quadrature (const Geometry::Rigid_Mesh<T> &rigid_mesh):
	 M_mesh (rigid_mesh), M_size (rigid_mesh.getCellsVector().size()+ rigid_mesh.getFractureFacetsIdsVector().size()),
	 M_quadrature(std::move(std::unique_ptr<QuadratureRule<T> >(new CentroidQuadrature<T>))), M_fractureQuadrature(std::move(std::unique_ptr<QuadratureRule<T> >(new CentroidQuadrature<T>)))
{}

template<class T>
Quadrature<T>::Quadrature (const Geometry::Rigid_Mesh<T> &rigid_mesh, const QuadratureRule<T>& quadrature):M_mesh (rigid_mesh), M_size (rigid_mesh.getCellsVector().size()+ rigid_mesh.getFractureFacetsIdsVector().size()), M_quadrature(std::move(quadrature.clone())), M_fractureQuadrature(std::move(std::unique_ptr<QuadratureRule<T> >(new CentroidQuadrature<T>)))
{}

template<class T>
Quadrature<T>::Quadrature (const Geometry::Rigid_Mesh<T> &rigid_mesh, const QuadratureRule<T>& quadrature, const QuadratureRule<T>& fracturequadrature):M_mesh (rigid_mesh), M_size (rigid_mesh.getCellsVector().size()+ rigid_mesh.getFractureFacetsIdsVector().size()), M_quadrature(std::move(quadrature.clone())), M_fractureQuadrature(std::move(fracturequadrature.clone())){}





// ==================================================
// Methods
// ==================================================


template<class T>
D_Real Quadrature<T>::Integrate(const Vector& Integrand)
{
	UInt IntSize = Integrand.size();
	if(IntSize != M_size)
		std::cerr << "ERROR: DIMENSION OF INTEGRAND FUNCTION DIFFERS FROM DIMENSION OF MESH" << std::endl;
	D_Real integral = 0;
	for (auto cell_it : M_mesh.getCellsVector())
	{		 
		integral += cell_it.getVolume()*Integrand(cell_it.getId());
	}	

	for (auto facet_it : M_mesh.getFractureFacetsIdsVector())
	{
		D_Real _volume = facet_it.Aperture()*facet_it.getFacet().size();
		integral += _volume*Integrand(facet_it.getIdasCell());
		integral -= _volume/2.*Integrand(facet_it.getSeparated()[0]);
		integral -= _volume/2.*Integrand(facet_it.getSeparated()[1]);
	}

	return integral;
}


template<class T>
D_Real Quadrature<T>::L2Norm(const Vector& Integrand)
{
	D_Real integral = 0;
	for (auto cell_it : M_mesh.getCellsVector())
	{		 
		integral += cell_it.getVolume()*Integrand(cell_it.getId())*Integrand(cell_it.getId());
	}	

	for (auto facet_it : M_mesh.getFractureFacetsIdsVector())
	{
		D_Real _volume = facet_it.Aperture()*facet_it.getFacet().size();
		integral += _volume*Integrand(facet_it.getIdasCell())*Integrand(facet_it.getIdasCell());
		integral -= _volume/2.*Integrand(facet_it.getSeparated()[0])*Integrand(facet_it.getSeparated()[0]);
		integral -= _volume/2.*Integrand(facet_it.getSeparated()[1])*Integrand(facet_it.getSeparated()[1]);
	}

	return sqrt(integral);
}


template<class T>
Vector Quadrature<T>::CellIntegrate (const std::function<D_Real(Generic_Point)>& func)
{
	UInt N = M_mesh.getCellsVector().size() + M_mesh.getFractureFacetsIdsVector().size();
	UInt counter = 0;
	Vector result(N);

	UInt NeighboursId;
	std::vector<Generic_Point> v_Nodes;

	for (auto cell_it : M_mesh.getCellsVector())
	{		 
		for (auto vertex_it : cell_it.getVertexesIds())
			v_Nodes.push_back (M_mesh.getNodesVector()[vertex_it]);

		result (counter) = M_quadrature->apply(v_Nodes, cell_it.getVolume(), func);

		v_Nodes.clear();
		++counter;
	}	

	for (auto facet_it : M_mesh.getFractureFacetsIdsVector())
	{
		D_Real _volume = facet_it.Aperture()*facet_it.getFacet().size();

		for (auto vertex_it : M_mesh.getFacetsVector()[facet_it.getId()].getVertexesIds())
			v_Nodes.push_back (M_mesh.getNodesVector()[vertex_it]);

		result (counter) = M_quadrature->apply(v_Nodes, _volume/2., func);
		v_Nodes.clear();

		NeighboursId = M_mesh.getFacetsVector()[facet_it.getId()].getSeparatedCellsIds()[0];
		result (NeighboursId) -= _volume/2.* result(NeighboursId)*M_mesh.getCellsVector()[NeighboursId].getVolume();

		NeighboursId = M_mesh.getFacetsVector()[facet_it.getId()].getSeparatedCellsIds()[1];
		result (NeighboursId) -= _volume/2.* result(NeighboursId)*M_mesh.getCellsVector()[NeighboursId].getVolume();
		
		++counter;
	}

	return result;
}	  

template<class T>
D_Real Quadrature<T>::Integrate(const std::function<D_Real(Generic_Point)>& Integrand)
{
	UInt N = M_mesh.getCellsVector().size() + M_mesh.getFractureFacetsIdsVector().size();
	Vector result(N);
	D_Real integral = 0;
	result = CellIntegrate(Integrand);

	for (UInt it = 0; it < N; ++it)
		integral += result (it);
	return integral;
}

//questa scelta non e' la migliore ma mi permette generalita'...se avessi solo triangoli potrei ottimizzare la scelta


}

#endif
