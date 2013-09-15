 /*!
 *	@file Fracture2D.hpp
 *	@brief Class for fracture in 2D space.
 *
 *	@author Luca Turconi <luca.turconi@moxoff.com>
 *
 */ 

#ifndef FRACTURE2D_HPP_
#define FRACTURE2D_HPP_

#include "TypeDefinition.hpp"

#include<CGAL/squared_distance_2.h>

#include<vector>
#include<string>
#include<cmath>

namespace Geometry{
		
	/*!
		@class Fracture2D
    	This class implements the concept of fracture in 2D space.
    */
class Fracture2D{
public:
	//! @name Constructor & Destructor
	//@{
	
	//! Empty constructor
	Fracture2D();
	
	//! Copy constructor
	Fracture2D(const Geometry::Fracture2D & f);
	
	//! Constructor, getting two points
	/*!
	 * @param source The source point
	 * @param target The target point
	 */
	Fracture2D(const Geometry::Point2D source, const Geometry::Point2D target);
	
	//! Destructor
	~Fracture2D();
	
	//@}
	
	//! @name Get Methods
	//@{
	
	//! Get source (const)
	/*!
	 * @return A constant reference to the source point
	 */
	const Geometry::Point2D & source() const
		{ return M_segment.source(); }
	
	//! Get target (const)
	/*!
	 * @return A constant reference to the target point
	 */
	const Geometry::Point2D & target() const
		{ return M_segment.target(); }
	
	//! Get segment (const)
	/*!
	 * @return A constant reference to the segment from source to target
	 */
	const Geometry::Segment2D & segment() const
		{ return M_segment; }
		
	//! Get segment
	/*!
	 * @return A reference to the segment from source to target
	 */
	Geometry::Segment2D & segment()
		{ return M_segment; }
	
	//! Get discretization step (const)
	/*!
	 * @return A constant reference to the discretization step
	 */
	const Real & h() const
		{ return M_h; }
	
	//! Get discretization step
	/*!
	 * @return A reference to the discretization step
	 */
	Real & h()
		{ return M_h; }
	
	//! Get the number of elements used for the fracture discretization (const)
	/*!
	 * @return A constant reference to the number of elements used for the fracture discretization
	 */
	const UInt & Nelements() const
		{ return M_Nelements; }
	
	//! Get the number of elements used for the fracture discretization
	/*!
	 * @return A reference to the number of elements used for the fracture discretization
	 */
	UInt & Nelements()
		{ return M_Nelements; }
	
	//! Get the vector of discretization nodes (const)
	/*!
	 * @return A constant reference to the vector M_nodes
	 */
	const std::vector<Geometry::Point2D> & nodes() const
		{ return M_nodes; }
		
	//! Get the vector of discretization nodes
	/*!
	 * @return A reference to the vector M_nodes
	 */
	std::vector<Geometry::Point2D> & nodes()
		{ return M_nodes; }

	//! Get the fracture permeability (const)
	/*!
	 * @return A constant reference to the fracture permeability
	 */
	const Real & permeability() const
		{ return M_permeability; }

	//! Get the fracture compressibility (const)
	/*!
	 * @return A constant reference to the fracture compressibility
	 */
	const Real & compressibility() const
		{ return M_compressibility; }

	//! Get the fracture aperture (const)
	/*!
	 * @return A constant reference to the fracture aperture
	 */
	const Real & aperture() const
		{ return M_aperture; }
	
	//! Get the fracture length
	/*!
	 * @return The fracture length
	 */
	Real length() const
		{ return std::sqrt(M_segment.squared_length()); }
	
	//@}
	
	//! @name Set Methods
	//@{
	
	//! Set the fracture permeability
	/*!
	 * @param perm The fracture permeability
	 */
	void setPermeability( const Real & perm )
		{ M_permeability = perm; }

	//! Set the fracture compressibility
	/*!
	 * @param compr The fracture compressibility
	 */
	void setCompressibility( const Real & compr )
		{ M_compressibility = compr; }
	
	//! Set the fracture aperture
	/*!
	 * @param aperture The fracture aperture
	 */
	void setAperture( const Real & aperture )
		{ M_aperture = aperture; }
	
	//@}
	
	//! @name Methods
	//@{
	
	//! Compute the distance between the fracture and the point p
	/*!
	 * @param p The point
	 * @return The distance between this and p
	 */
	Real distance(const Geometry::Point2D & p) const
		{ return CGAL::squared_distance(p,M_segment); }
	
	//! Build the fracture discretization using elements of dimension as similar as possible to h
	/*!
	 * @param h The desidered elements dimension
	 * @return The real element dimension used for the discretization
	 */
	Real buildFractureDiscretization(const Real & h);
	
	//! Build the fracture discretization using a fixed number of elements
	/*!
	 * @param Nelements The numeber of elements to be used for the discretization
	 * @return The dimension of the elements used for the discretization
	 */
	Real buildFractureDiscretization(const UInt & Nelements);
	
	//! Export in vtk format the segment representing the fracture
	/*!
	 * (Use paraview to open the vtk file)
	 * @param filename The name of the vtk file created by this method
	 * @return TRUE -> operation ended correctly
				FALSE -> an error occurred
	 */
	bool exportVtk(const std::string & filename) const;
	
	//! Export in vtk format the set of the nodes used for the fracture discretization
	/*!
	 * (Use paraview to open the vtk file)
	 * @param filename The name of the vtk file created by this method
	 * @return TRUE -> operation ended correctly
				FALSE -> an error occurred
	 */
	bool exportFractureDiscretizationVtk(const std::string & filename) const;
	
	//! Display general information about the content of the class
	/*!
	 * @param out Specify the output format (std::cout by default)
	 */	
	void showMe(std::ostream & out=std::cout) const;
	
	//@}
	
private:

	//! It stores the segment representing the fracture
	Geometry::Segment2D M_segment;
	//! It stores the discretization step
	Real M_h;
	//! It stores the number of nodes used for the discretization
	UInt M_Nelements;
	//! It stores the nodes used for the discretization
	std::vector<Geometry::Point2D> M_nodes;
	
	//! @name Physical properties
	//@{
	//! The fracture permeability
	Real M_permeability;
	//! The fracture compressibility
	Real M_compressibility;
	//! The fracture aperture
	Real M_aperture;
	//@}
};

} // namespace Geometry


#endif /* FRACTURE2D_HPP_ */
