/*!
 *	@file BC.hpp
 *	@brief This class handles the boundary conditions of the Darcy problem.
 */ 

#ifndef __DARCYBC_HPP__
#define __DARCYBC_HPP__

#include <vector>
#include <cmath>
#include <assert.h>
#include <functional>
#include <algorithm>
#include <iostream>

class Mesh3D;

namespace Darcy{

//! Select the type of the BC
/*!
	@enum BCType
	It is possible to choose the type of the boundary conditions: "Dirichlet" type or "Neumann" type.
*/
enum BCType{Dirichlet, Neumann};

//! Class for implementing the BCs
/*!
	@class BoundaryConditions
	This class implements the boundary conditions.
*/
class BoundaryConditions
{
	typedef Geometry::Point3D Generic_Point;
	typedef Geometry::Mesh3D::Facet3D Generic_Border;

public:

	//! Class that implements the BC on a portion of the border
	/*!
		@class BorderBC
		This class implements the boundary conditions on a portion of the domain.
	*/
	class BorderBC {
	protected:
		//! Id of the BC
		UInt m_Id;
		//! Type of Boundary Condition as BCType: Dirichlet or Neumann are implemented
		BCType bcType;
		//! Boundary condition as a function
		std::function<Real(Generic_Point)> BC;
		//! Pointer to the container of the BorderBC: BoundaryConditions
		BoundaryConditions * m_bcContainer;

	public:
		//! @name Constructor & Destructor
		//@{

		//! Constructor for a BorderBC, given a border-Id, a BCType and a function.
		/*!
			@param Id the Id of the borders on which we want to impose the boundary condition.
			@param bctype is of type BCType and can be Dirichlet or Neumann
			@param bc is the function that is actually our BC
		*/
		BorderBC(UInt Id, BCType bctype, std::function<Real(Generic_Point)> bc):
			m_Id(Id), bcType(bctype), BC(bc), m_bcContainer(0) {}
		//! Default Copy Constructor
		BorderBC(const BorderBC &) = default;
		//! Default Destructor
		~BorderBC() = default;
		//@}

		//! @name Get Methods
		//@{
		//! Get Id (const)
		/*!
		 * @return The Id of the Border
		 */
		UInt getId () const
			{ return m_Id; }
		//! Get container (const)
		/*!
		 * @return A const pointer to the containing BoundaryConditions object
		 */
		BoundaryConditions * getContainer () const
			{ return m_bcContainer; }
		//! Get BCType (const)
		/*!
		 * @return Dirichlet or Neumann
		 */
		BCType getBCType () const
			{ return bcType; }
		//! Get Bc (const)
		/*!
		 * @return A function which is actually the Boundary condition
		 */
		const std::function<Real(Generic_Point)> getBC() const
			{ return BC; }
		//@}

		friend class BoundaryConditions;
	};

	//! @name Get Methods
	//@{
	//! Get Borders Number (const)
	/*!
	 * @return The number of borders of the domain
	 */
	UInt getBorderConditionsNumber() const
		{ return BordersBCMap.size(); }
	//! Get BorderBC vector (const)
	/*!
	 * @return A reference to the vector of BorderBC
	 */
	const std::map<UInt,BorderBC> & getBordersBCMap() const
		{ return BordersBCMap; }
	//@}

	//! @name Constructor & Destructor
	//@{

	//! Constructor for a BoundaryConditions, given a vector of BorderBC.
	/*!
		@param borderbc is a vector of BorderBC containing the boundary conditions.
	*/
	BoundaryConditions(std::vector<BorderBC> & borderbc);
	//! DEfault copy constructor
	BoundaryConditions(const BoundaryConditions &) = default;
	//! Default destructor
	~BoundaryConditions() = default;
	//@}

protected:
	//! Map of BorderBC. First -> BC id, second -> BorderBC
	std::map<UInt,BorderBC> BordersBCMap;
};

}// namespace Darcy

#endif
