 /*!
 *	@file PolygonalDomain.hpp
 *	@brief Class for polygonal domains.
 *
 *	@author Francesco Della Porta 
 *
 */ 

#ifndef __DOMAIN_HPP__
#define __DOMAIN_HPP__

#include <vector>
#include <iostream>
#include "Dimension.hpp"
	

namespace Geometry{

/*!
	@class Domain
	This class contains the borders of the domain. It is necessary in order to impose Boundary Conditions to PDE.
	Thanks to the template parameter T this class works for bi- and tri-dimensional domains
*/
template <class T> 
class Domain
{
	typedef typename T::Generic_Point Generic_Point;
	typedef typename T::Generic_Segment Generic_Segment;

	/*!
		@typedef Generic_Facet
   		This typedefinition permits to treat a vector of Point3D in 3D and a Segment2D in 2D as a Generic_Border.
   	*/
	typedef typename T::Generic_Border Generic_Border;

public:

	/*!
		@class DomainBorder
		This class describes one borders of the domain.
	*/
	class DomainBorder {
	protected:
		//! The Id of the Border
		UInt m_Id;
		//! The Border: a segment in 2D, a vector of points in 3D. In 3D the point in position "i" is connected by an edge to the point in position "i+1" and "i-1"
		Generic_Border borderFacet;
	public:

		//! @name Constructor & Destructor
		//@{

		//! Constructor for a DomainBorder given a Generic_Border.
		/*!
			@param gBorder a segment in 2D, a vector of points in 3D. In 3D the point in position "i" is connected by an edge to the point in position "i+1" and "i-1"
		*/
		DomainBorder (const Generic_Border gBorder);
		//! Copy Constructor for a DomainBorder.
		DomainBorder (const DomainBorder& domainborder);
		//! Default destructor.
		~DomainBorder () = default;
		//@}

		//! @name Get Methods
		//@{		
		//! Get Id (const)
		/*!
		 * @return The Id of the Border
		 */
		const UInt getId () const
			{return m_Id;}
		//! Get Border (const)
		/*!
		 * @return The Generic_Border
		 */
		const Generic_Border getBorder() const
			{return borderFacet;}
		//@}

		friend class Domain;
	private:
		//! @name Private Methods
		//@{
		
		//! Is a function called by the copy-constructor in the 3D case
		void M_copyConstructor (const std::vector<Geometry::Point3D> border);
		//! Is a function called by the copy-constructor in the 2D case
		void M_copyConstructor (const Geometry::Segment2D border);
		//@}
	};

	//! @name Get Methods
	//@{
	//! Get Borders numbers (const)
	/*!
	 * @return The number of borders of the domain
	 */
	const UInt getBordersnumber() const
		{return BordersVector.size();}
	//! Get Borders vector (const)
	/*!
	 * @return The vector of borders of the domain
	 */
	const std::vector<DomainBorder>& getBordersVector() const
		{return BordersVector;}

	//! Get Border Id (const)
	/*!
	 * @return The Id of a border given a vector of Points in the domain-border in 2D case
	 */
	const UInt getBorderId(const std::vector<Point2D> facetVertexes) const;
	//! Get Border Id (const)
	/*!
	 * @return The Id of a border given a vector of Points in the domain-border in 3D case
	 */
	const UInt getBorderId(const std::vector<Point3D> facetVertexes) const;
	//@}


	//! @name Constructor & Destructor
	//@{
	//! Constructor for a Domain given a vector of Generic_Border.
	/*!
		@param border is a vector containing the DomainBorder of the Domain 
	*/
	Domain(const std::vector<DomainBorder>& border);
	//! Default Copy Constructor	
	Domain(const Domain<T> & domain) = default;	
	//! Default destructor
	~Domain() = default;	
	//@}

protected:
	//! Vector containing DomainBorder
	std::vector<DomainBorder> BordersVector;
};


// --------------------   Class Domain   --------------------
 
// ==================================================
// Constructors & Destructor
// ==================================================


template <class T> 
Domain<T>::Domain(const std::vector<DomainBorder>& border):
	BordersVector (border)
{
	for(UInt N = 0; N < BordersVector.size(); ++N)
		BordersVector[N].m_Id = N;
}

 
// ==================================================
// Get Methods
// ==================================================

//achtung in 3d bisogna stare attenti al segmento che gli si passa!
template <class T> 
const UInt Domain<T> :: getBorderId(const std::vector<Point2D> facetVertexes) const
{
	UInt IdCounter = 0;
	for (auto iter = BordersVector.begin(); iter != BordersVector.end(); ++iter)
	{
		if(iter->getBorder().has_on(facetVertexes[0]))	
			if(iter->getBorder().has_on(facetVertexes[0]))	
				return 	IdCounter;
		++IdCounter;
	}
	std::cerr << "Border not Found" << std::endl;
	return -1;
}


template <class T> 
const UInt Domain<T> ::getBorderId(const std::vector<Point3D> facetVertexes) const
{
	UInt IdCounter = 0;
	for (auto iter = BordersVector.begin(); iter != BordersVector.end(); ++iter)
	{
		if(coplanar(iter->getBorder()[0], iter->getBorder()[1], iter->getBorder()[2], facetVertexes[0]))	
			if(coplanar(iter->getBorder()[0], iter->getBorder()[1], iter->getBorder()[2], facetVertexes[1]))
				if(coplanar(iter->getBorder()[0], iter->getBorder()[1], iter->getBorder()[2], facetVertexes[2]))
					return 	IdCounter;
		++IdCounter;
	}
	std::cerr << "Border not Found" << std::endl;
	return -1;
}


// --------------------   Class DomainBorder   --------------------
 
// ==================================================
// Constructors & Destructor
// ==================================================

template <class T> 
Domain<T> ::DomainBorder::DomainBorder (const Generic_Border gBorder): borderFacet(gBorder)
{}

template <class T> 
Domain<T> ::DomainBorder::DomainBorder (const DomainBorder& domainborder): m_Id(domainborder.getId())
{
	M_copyConstructor(domainborder.getBorder());
}

 
// ==================================================
// Private Methods
// ==================================================

template <class T> 
void Domain<T> ::DomainBorder::M_copyConstructor (const Geometry::Segment2D border)
{
	Geometry::Segment2D segmentcopy(border.source(), border.target());
	borderFacet = segmentcopy;
}

template <class T> 
void Domain<T> ::DomainBorder::M_copyConstructor (const std::vector<Geometry::Point3D> border)
{
	for (auto it : border)
		borderFacet.emplace_back(it);
}


}

#endif
