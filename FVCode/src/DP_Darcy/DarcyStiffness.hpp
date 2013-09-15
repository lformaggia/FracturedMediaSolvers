 /*!
 *	@file DarcyStiffness.hpp
 *	@brief Class for building a Stiffness-matrix for finite volume discretization of the Darcy problem.
 *
 *	@author Francesco Della Porta 
 *
 */ 

#ifndef __DARCYSTIFFNESS_HPP__
#define __DARCYSTIFFNESS_HPP__

#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>
#include "DarcyTypeDefinitions.hpp"
#include "MatrixHandler.hpp"
#include "DarcyBC.hpp"

namespace Darcy{

/*!
	@class StiffMatrix
	This class constructs the stiffness-matrix for a Darcy problem. The adopted technique is a two point finite volume method. Also fractures represented by Facets are considered as cells and take part to discretization.   
*/
template <class T> 
class StiffMatrix: public MatrixHandler<T> {

	typedef typename T::Generic_Point Generic_Point;
	typedef typename T::Generic_Vector Generic_Vector;
	typedef typename T::Fracture_Juncture Fracture_Juncture;
	typedef typename Geometry::Rigid_Mesh<T>::Facet_ID Facet_ID;

public:
	//! @name Constructor & Destructor
	//@{

	//! Constructor for a stiffness-Matrix, given a Geometry::Rigid_Mesh, the mobility of the fluid, boundary conditions and the permeability of the material
	/*!
		@param rigid_mesh A Geometry::Rigid_Mesh on which the matrix is constructed
		@param BC Boundary condition given in the container Darcy::BoundaryConditions
		@param permeability The permeabibility of the ground is given through a function F: Generic_Point->D_Real
		@param mobility is the mobility of the considered fluid under the ground
	*/
	StiffMatrix(const Geometry::Rigid_Mesh<T> &rigid_mesh, BoundaryConditions<T>& Bc, std::function<D_Real(Generic_Point)> permeability, D_Real mobility = 1.);
	//! Copy-Constructor deleted
	StiffMatrix(const StiffMatrix&) = delete;
	//! Empty-Constructor deleted
	StiffMatrix() = delete;
	//! Default destructor
	~StiffMatrix() = default;
	//@}


	//! @name Get Methods
	//@{
	//! Get BC vector (const)
	/*!
	 * @return A reference to const vector wich has to be summed to f in the linear-system Ab=f. This vector is important in order to impose the BC.
	 */
	const Vector& getBCVector() const
		{return *_b;}
	//@}

	//! @name Methods
	//@{
	//! assemble
	/*!
	 * @return Constructs the Mass matrix
	 */
	void assemble();
	//@}

protected:

	//! @name Protected Methods
	//@{
		
	//! vector product in 3D
	/*!
	 * @param tangent_1 A 3D vector
	 * @param tangent_2 A 3D vector
	 * @return The vector product in 3D tangent_1 x tangent_2
	 */
	Geometry::Vector3D vector_product (Geometry::Vector3D tangent_1, Geometry::Vector3D tangent_2) const
		{return (tangent_1.y()*tangent_2.z() - tangent_1.z()*tangent_2.y(),
				tangent_1.z()*tangent_2.x() - tangent_1.x()*tangent_2.z(),
				tangent_1.x()*tangent_2.y() - tangent_1.y()*tangent_2.x());};

	//! scalar product
	/*!
	 * @param vect1 A vector
	 * @param vect2 A vector
	 * @return The scalar product between vect1 and vect2
	 */
	Real scalar_product (Generic_Vector vect1, Generic_Vector vect2) const 
		{return (vect1*vect2);};

	//! vector lenght
	/*!
	 * @param vect A vector
	 * @return The lenght of vect
	 */
	Real vector_lenght (Generic_Vector vect) const
		 {return sqrt(vect.square_lenght());};

	//! border center in 2D
	/*!
	 * @param fj the Id of the point which is a juncture of two fractures in 2D
	 * @return The center of the juncture between two Fracure_Facet (in 2D it is the point itself)
	 */
	Generic_Point border_center (UInt fj) const
		 {return this->M_mesh.getNodesVector()[fj];};

	//! border center in 3D
	/*!
	 * @param fj the Id of the juncture of two Fracture_Facet in 3D
	 * @return The center of the juncture between two Fracure_Facet
	 */
	Generic_Point border_center (std::pair<UInt,UInt> fj) const
		 {return (this->M_mesh.getNodesVector()[fj.first] + (this->M_mesh.getNodesVector()[fj.second] - this->M_mesh.getNodesVector()[fj.first])/2.);};

	//! It is called by the method assemble() and it computes the coefficient alpha
	/*!
	 * @param cellId the Id of a Cell
	 * @param facet A pointer to a Geometry::Rigid_mesh::Facet_ID
	 * @return The computed coefficient alpha
	 */
	D_Real Findalpha (const UInt& cellId, Facet_ID *const facet) const;

	//! It is called by the method assemble() and it computes the coefficient alpha in the case of Dirichlet Bc
	/*!
	 * @param cellId the Id of a Cell
	 * @param facet A pointer to a Geometry::Rigid_mesh::Facet_ID
	 * @return The computed coefficient alpha
	 */
	D_Real FindDirichletalpha (const UInt& cellId, Facet_ID *const facet) const;

	//! It is called by the method assemble() and it computes the coefficient alpha in the case of a fracture in 2D
	/*!
	 * @param fj is a Fracture_Juncture
	 * @param n_Id The Id of the Fracture_Facet
	 * @return The computed coefficient alpha
	 */
	D_Real Findfracturesalpha (const UInt fj, const UInt n_Id) const;

	//! It is called by the method assemble() and it computes the coefficient alpha in the case of a fracture in 2D
	/*!
	 * @param fj is a Fracture_Juncture
	 * @param n_Id The Id of the Fracture_Facet
	 * @return The computed coefficient alpha
	 */
	D_Real Findfracturesalpha (const std::pair<UInt,UInt> fj, const UInt n_Id) const;
	//@}
	
protected:
	//! The vector with the BC
	std::unique_ptr<Vector> _b;
	//! A function for the permeability of the ground
	std::function<D_Real(Generic_Point)> M_Permeability;
	//! The mobility of the fluid
	D_Real M_Mobility;
	//! The container with the BC
	BoundaryConditions<T>& m_Bc;
};


// --------------------   Class StiffMatrix   --------------------
 
// ==================================================
// Constructors & Destructor
// ==================================================


template<class T>
StiffMatrix<T>::StiffMatrix (const Geometry::Rigid_Mesh<T> &rigid_mesh, BoundaryConditions<T>& Bc,
	 std::function<D_Real(Generic_Point)> permeability, D_Real mobility): MatrixHandler<T>(rigid_mesh),
	   _b (new Vector(this->M_size)), M_Permeability(permeability), M_Mobility(mobility), m_Bc(Bc)
{}


// ==================================================
// Methods
// ==================================================

template<class T>
void StiffMatrix<T>::assemble()
{
	std::vector<Triplet> Matrix_elements;
	D_Real alpha1, alpha2, alphaF, alphaf;
	UInt neighbor1id, neighbor2id;
	std::vector<D_Real> alphas;
	Generic_Vector f1, f2;
	Real Df;
	D_Real T12, Q12, T1f, Q1f, T2f, Q2f, QFf, Q1o;

	for (auto facet_it : this->M_mesh.getInternalFacetsIdsVector())
	{		 
		neighbor1id = facet_it.getSeparated()[0];
		neighbor2id = facet_it.getSeparated()[1];

		alpha1 = Findalpha(neighbor1id, &facet_it);	
		alpha2 = Findalpha(neighbor2id, &facet_it);	
		
		T12 = alpha1*alpha2/(alpha1 + alpha2);
		Q12 = T12*M_Mobility;

		Matrix_elements.emplace_back(Triplet (neighbor1id, neighbor1id, Q12));
		Matrix_elements.emplace_back(Triplet (neighbor1id, neighbor2id, -Q12));
		Matrix_elements.emplace_back(Triplet (neighbor2id, neighbor1id, -Q12));
		Matrix_elements.emplace_back(Triplet (neighbor2id, neighbor2id, Q12));
	}	

	for (auto facet_it : this->M_mesh.getFractureFacetsIdsVector())
	{
		D_Real F_permeability = facet_it.Permeability();
		D_Real F_aperture = facet_it.Aperture();
		Df = F_aperture/2.;
		alphaf = facet_it.getSize()*F_permeability/Df;
		
		neighbor1id = facet_it.getSeparated()[0];
		neighbor2id = facet_it.getSeparated()[1];

		alpha1 = Findalpha(neighbor1id, &facet_it);	
		alpha2 = Findalpha(neighbor2id, &facet_it);	

		T1f = alpha1*alphaf/(alpha1 + alphaf);
		Q1f = T1f*M_Mobility;

		T2f = alphaf*alpha2/(alphaf + alpha2);
		Q2f = T2f*M_Mobility;

		Matrix_elements.emplace_back(Triplet (neighbor1id, neighbor1id, Q1f));
		Matrix_elements.emplace_back(Triplet (neighbor1id, facet_it.getIdasCell(), -Q1f));
		Matrix_elements.emplace_back(Triplet (facet_it.getIdasCell(), neighbor1id, -Q1f));
		Matrix_elements.emplace_back(Triplet (facet_it.getIdasCell(), facet_it.getIdasCell(), Q1f));

		Matrix_elements.emplace_back(Triplet (neighbor2id, neighbor2id, Q2f));
		Matrix_elements.emplace_back(Triplet (neighbor2id, facet_it.getIdasCell(), -Q2f));
		Matrix_elements.emplace_back(Triplet (facet_it.getIdasCell(), neighbor2id, -Q2f));
		Matrix_elements.emplace_back(Triplet (facet_it.getIdasCell(), facet_it.getIdasCell(), Q2f));

		for (auto facetborders_it : facet_it.getFractureNeighbors ())
		{
			alphaF = Findfracturesalpha (facetborders_it.first, facet_it.getId());

			for (auto neighbors_it : facetborders_it.second)
				alphas.emplace_back(Findfracturesalpha (facetborders_it.first, neighbors_it));
			
			D_Real a_sum = alphaF;
			auto sum_maker = [&a_sum](D_Real item){a_sum += item;};
			std::for_each(alphas.begin(), alphas.end(), sum_maker);

			#ifdef __MY__DEBUG__
			assert(alphas.size() == facetborders_it.second.size());
			D_Real sum2 = alphaF;
			for (auto it : alphas)
				sum2 += it;
			assert(a_sum == sum2);
			#endif

			for (UInt counter = 0; counter < alphas.size(); ++counter)
			{
				QFf = alphaF*alphas[counter]*M_Mobility/a_sum;
				Matrix_elements.emplace_back(Triplet (facet_it.getIdasCell(), facet_it.getIdasCell(), QFf));
				Matrix_elements.emplace_back(Triplet (facet_it.getIdasCell(), this->M_mesh.getFractureFacetsIdsVector()[facetborders_it.second[counter]].getIdasCell(), -QFf));
			}
			alphas.clear();
		}
	}

	for(UInt i=0; i<this->M_size; ++i)
		_b->operator()(i)=0;

	for (auto facet_it : this->M_mesh.getBorderFacetsIdsVector())
	{
		neighbor1id = facet_it.getSeparated()[0];
		UInt borderId = facet_it.getBorderId();

		if(m_Bc.getBordersBCVector()[borderId].getBCType() == Neumann )
		{
			Q1o = m_Bc.getBordersBCVector()[borderId].getBC()(facet_it.getCenter())*facet_it.getSize();
			_b->operator()(neighbor1id) += Q1o;
		}

		if(m_Bc.getBordersBCVector()[borderId].getBCType() == Dirichlet)
		{
			alpha1 = Findalpha (neighbor1id, &facet_it);	
			alpha2 = FindDirichletalpha (neighbor1id, &facet_it);

			T12 = alpha1*alpha2/(alpha1 + alpha2);
			Q12 = T12*M_Mobility;
			Q1o = Q12*m_Bc.getBordersBCVector()[borderId].getBC()(facet_it.getCenter());

			Matrix_elements.emplace_back(Triplet (neighbor1id, neighbor1id, Q12));
			_b->operator()(neighbor1id) = _b->operator()(neighbor1id) + Q1o;
		}				
	}

	this->M_Matrix->setFromTriplets(Matrix_elements.begin(), Matrix_elements.end());
}


// ==================================================
// Protected Methods
// ==================================================

template<class T>
D_Real StiffMatrix<T>::Findfracturesalpha (const UInt fj, const UInt n_Id) const
{
	Generic_Point borderCenter = border_center(fj);
	Generic_Point cellCenter = this->M_mesh.getFractureFacetsIdsVector()[n_Id].getCenter();
	D_Real A = this->M_mesh.getFractureFacetsIdsVector()[n_Id].Aperture();
	D_Real k = this->M_mesh.getFractureFacetsIdsVector()[n_Id].Permeability();
	Generic_Vector f;
	D_Real alpha;
	Real D;
	f = borderCenter - cellCenter;
	D = sqrt(scalar_product(f, f));
	alpha = A*k/D;
	return alpha;
}	



template<class T>
D_Real StiffMatrix<T>::Findfracturesalpha (const std::pair<UInt,UInt> fj, const UInt n_Id) const
{
	Generic_Point borderCenter = border_center(fj);
	Generic_Point cellCenter = this->M_mesh.getFractureFacetsIdsVector()[n_Id].getCenter();
	D_Real A = this->M_mesh.getFractureFacetsIdsVector()[n_Id].Aperture()*vector_lenght( this->M_mesh.getNodesVector()[fj.second]-this->M_mesh.getNodesVector()[fj.first]);
	D_Real k = this->M_mesh.getFractureFacetsIdsVector()[n_Id].Permeability();
	Generic_Vector f;
	D_Real alpha;
	Real D;
	f = fj - cellCenter;
	D = sqrt(scalar_product(f, f));
	f /= D;	
	Generic_Vector normal = vector_product(f, this->M_mesh.getNodesVector()[fj.second]-this->M_mesh.getNodesVector()[fj.first]);
	normal = vector_product (this->M_mesh.getNodesVector()[fj.second]-this->M_mesh.getNodesVector()[fj.first], normal);
	normal /= vector_lenght(normal);
	D_Real scalprod = fabs(scalar_product(normal, f));
	alpha = A*k*scalprod/D;
	return alpha;
}	


template<class T>
D_Real StiffMatrix<T>::Findalpha (const UInt& cellId, Facet_ID *const facet) const
{
	Generic_Point facetCenter = facet->getCenter();
	Generic_Point cellCenter = this->M_mesh.getCellsVector()[cellId].getCentroid();
 	Generic_Vector normal = facet->getUNormal();
	D_Real alpha;
	Real D;
	Generic_Vector f;
	Real scalprod;
	f = cellCenter - facetCenter;
	D = sqrt(scalar_product(f, f));
	f = f/D;

	#ifdef __MY__DEBUG__
	assert(scalar_product(f, f) < 1. + 1.e-14);
	assert(scalar_product(f, f) > 1. - 1.e-14);
	#endif

	scalprod = fabs(scalar_product(f, normal));
	alpha = facet->getSize()*M_Permeability(cellCenter)*scalprod/D;
	return alpha;	
}

template<class T>
D_Real StiffMatrix<T>::FindDirichletalpha (const UInt& cellId, Facet_ID *const facet) const
{
	Generic_Point facetCenter = facet->getCenter();
	Generic_Point cellCenter = this->M_mesh.getCellsVector()[cellId].getCentroid();
	D_Real alpha;
	Real D;
	Generic_Vector f;
	f = cellCenter - facetCenter;
	D = sqrt(scalar_product(f, f));

	alpha = facet->getSize()*M_Permeability(cellCenter)/D;
	return alpha;	
}


}

#endif
