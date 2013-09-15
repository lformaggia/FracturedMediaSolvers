 /*!
 *	@file DarcyBC.hpp
 *	@brief Class to handle Boundary Conditions of Differential Problems.
 *
 *	@author Francesco Della Porta 
 *
 */ 

#ifndef __DARCYBC_HPP__
#define __DARCYBC_HPP__

#include <vector>
#include <cmath>
#include <assert.h>
#include <functional>
#include <algorithm>
#include <iostream>
#include "DarcyTypeDefinitions.hpp"
#include "../DP_src/PolygonalDomain.hpp"
#include "../DP_src/Dimension.hpp"

	

namespace Darcy{

/*!
	@enum BCType
	It is possible to choose if the Boundary conditions are of "Dirichlet" type or of "Neumann" type.
*/
enum BCType{Dirichlet, Neumann};

/*!
	@class BoundaryConditions
	This class implements the boundary consitions on a class of type Geometry::Domain. 
	Thanks to the template parameter T ths class works on bi- and tri-dimensional domains.
*/
template <class T> 
class BoundaryConditions
{
	typedef typename T::Generic_Point Generic_Point;
	typedef typename T::Generic_Segment Generic_Segment;
	typedef typename T::Generic_Border Generic_Border;

public:

	/*!
		@class BorderBC
		This class implements the boundary consitions on a given border of thedomain. 
	*/
	class BorderBC {
	protected:
		//! Id of the domain-border 
		UInt m_Id;
		//! Type of Boundary Condition of BCType: Dirichlet or Neumann are implemented
		BCType bcType;
		//! Boundary condition, hence a function
		std::function<D_Real(Generic_Point)> BC;
		//! Pointer to the container of the BorderBC: BoundaryConditions
		BoundaryConditions<T>* m_bcContainer;

	public:
		//! @name Constructor & Destructor
		//@{

		//! Constructor for a BorderBC, given a border-Id, a BCType and a function.
		/*!
			@param Id is the Id of the border on which we want to impose the boundary condition.
			@param bctype is of type BCType and can be Dirichlet or Neumann
			@param bc is the function that is effectively our BC
		*/
		BorderBC (UInt Id, BCType bctype, std::function<D_Real(Generic_Point)> bc);
		//! Default Copy Constructor
		BorderBC (const BorderBC&) = default;
		//! Default Destructor
		~BorderBC () = default;
		//@}

		//! @name Get Methods
		//@{
		//! Get Id (const)
		/*!
		 * @return A the Id of the Border
		 */
		const UInt getId () const
			{return m_Id;}
		//! Get container (const)
		/*!
		 * @return A const pointer to the containing BoundaryConditions object
		 */
		BoundaryConditions<T>* const getContainer () const
			{return m_bcContainer;}
		//! Get BCType (const)
		/*!
		 * @return Dirichlet or Neumann
		 */
		const BCType getBCType () const
			{return bcType;}
		//! Get border (const)
		/*!
		 * @return the Generic_Border related to this BC
		 */
		const Generic_Border getBorder() const
			{return m_bcContainer->getBordersVector()[m_Id];}
		//! Get Bc (const)
		/*!
		 * @return A function which is effectively the Boundary condition
		 */
		const std::function<D_Real(Generic_Point)> getBC() const
			{return BC;}
		//@}

	friend class BoundaryConditions;
	};

	//! @name Get Methods
	//@{
	//! Get Borders Number (const)
	/*!
	 * @return The number of borders of the domain
	 */
	const UInt getBordersnumber() const
		{return BordersBCVector.size();}
	//! Get BorderBC vector (const)
	/*!
	 * @return A reference to the vector of BorderBC
	 */
	const std::vector<BorderBC>& getBordersBCVector() const
		{return BordersBCVector;}
	//@}

	//! @name Constructor & Destructor
	//@{

	//! Constructor for a BoundaryConditions, given a vector of BorderBC and a Geometry::Domain.
	/*!
		@param borderbc is a vector of BorderBC containing the boundary conditions.
		@param domain is the domain contained in the class Geometry::Domain
	*/
	BoundaryConditions (std::vector<BorderBC>& borderbc, Geometry::Domain<T>& domain);	
	//! DEfault copy constructor
	BoundaryConditions (const BoundaryConditions&) = default;
	//! Default destructor
	~BoundaryConditions () = default;
	//@}

protected:
	//! vector of BorderBC
	std::vector<BorderBC> BordersBCVector;
	//! Domain
	Geometry::Domain<T>& M_domain; 
};


// --------------------   Class BoundaryConditions   --------------------
 
// ==================================================
// Constructors & Destructor
// ==================================================

template <class T> 
BoundaryConditions<T>::BoundaryConditions(std::vector<BorderBC>& borderbc, Geometry::Domain<T>& domain):
	BordersBCVector (borderbc), M_domain(domain)
{
	assert(BordersBCVector.size() == M_domain.getBordersVector().size());

	for (auto it : BordersBCVector)
		it.m_bcContainer = this;
	
	auto cmp = [](BorderBC bc1, BorderBC bc2){return (bc1.getId()<bc2.getId());};
	std::sort (BordersBCVector.begin(),BordersBCVector.end(), cmp);

	for (UInt i = 0; i < BordersBCVector.size(); ++i)
		if (BordersBCVector[i].getId() != i)
			std::cerr << "ERRORE: Imposta due volte una BC sullo stesso lato..." << std::endl;
}


// --------------------   Class BorderBC   --------------------
 
// ==================================================
// Constructors & Destructor
// ==================================================

template <class T> 
BoundaryConditions<T> ::BorderBC::BorderBC (UInt Id, BCType bctype, std::function<D_Real(Generic_Point)> bc): m_Id(Id), bcType(bctype), BC(bc)
{}

}

#endif
