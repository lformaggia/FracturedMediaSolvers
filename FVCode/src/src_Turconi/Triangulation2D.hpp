 /*!
 *	@file Triangulation2D.hpp
 *	@brief Class for unstructured triangular mesh in 2D space.
 *
 *	@author Luca Turconi <luca.turconi@moxoff.com>
 *
 */
 
#ifndef TRIANGULATION2D_HPP_
#define TRIANGULATION2D_HPP_

#include "TypeDefinition.hpp"
#include "FractureNetwork2D.hpp"
#include "Delaunay_mesh_adaptiveSize_criteria_2.hpp"

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

#include <iostream>
#include <utility>
 
namespace Geometry{

	/*!
		@class Triangulation2D
    	This class implements the concept of triangulation in 2D space.
    */
class Triangulation2D{
public:
// -----------------------------------
	typedef Kernel K;
// -----------------------------------
	typedef CGAL::Triangulation_vertex_base_2<K> Vb;
	typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
	typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
// -----------------------------------
//	typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds, CGAL::Exact_intersections_tag> CDT;
	typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds, CGAL::Exact_predicates_tag> CDT; 
// -----------------------------------
	typedef CGAL::Constrained_triangulation_plus_2<CDT> CDTplus;
	typedef CGAL::Constrained_triangulation_plus_2<CDT>::Vertices_in_constraint_iterator
			Vertices_in_constraint_iterator;

	typedef CGAL::Delaunay_mesh_size_criteria_2<CDTplus> Criteria;
	typedef CGAL::Delaunay_mesh_adaptiveSize_criteria_2<CDTplus> AdaptiveCriteria;
	typedef CGAL::Delaunay_mesher_2<CDTplus, Criteria> Mesher;
	typedef CGAL::Delaunay_mesher_2<CDTplus, AdaptiveCriteria> AdaptiveMesher;

	typedef CDT::Vertex_handle Vertex_handle;
	typedef CDT::Face_handle Face_handle;
	typedef CDT::Point Point;
	
	typedef std::vector<Vertex_handle> ApproxFracture;
	//! A pair of Vertex_handle pointing to the vertexes of a constraint segment
	typedef std::pair<Vertex_handle,Vertex_handle> FnConstraint;

public:
	//! @name Constructor & Destructor
	//@{
	
	//! Empty constructor
	Triangulation2D();
	
	//! Copy constructor
	Triangulation2D(const Triangulation2D & t);
	
	//! Constructor, getting the border vertexes
	/*!
	 * @param borderNodes A vector containing the vertexes of the border
	 */
	Triangulation2D(const std::vector<Geometry::Point2D> & borderNodes);
	
	//! Destructor
	~Triangulation2D();
	
	//@}
	
	//! @name Get Methods
	//@{
	
	//! Get the object storing the CGAL triangulation (const)
	/*!
	 * @return A constant reference to the triangulation
	 */
	const CDTplus & getTriangulation() const
		{ return M_cdt; }

	//! Get the object storing the CGAL triangulation
	/*!
	 * @return A reference to the triangulation
	 */
	CDTplus & getTriangulation()
		{ return M_cdt; }
	
	//! Get the associated fracture network (const)
	/*!
	 * @return A constant reference to the associated fracture network
	 */
	const Geometry::FractureNetwork2D & getFractureNetwork() const
		{ return M_fn; }	
		
	//! Get the vector of the constraints representing the associated fracture network
	/*!
	 * @return A constant reference to the vector of the constraint
	 */
	const std::vector<FnConstraint> & getFractureNetworkConstraints() const
		{ return M_fnConstraints; }
		
	//! Get the vector of the approximation of the associated fracture network (const)
	/*!
	 * @return A constant reference to the approximation of the associated fracture network
	 */
	const std::vector<ApproxFracture> & getApproxFractureNetwork() const
		{ return M_approxFractureNetwork; }
		
	//@}
	
	//! @name Methods
	//@{
	
	//! The number of elements of the triangulation
	/*!
	 * @return The number of elements of the triangulation
	 */
	UInt Nelements() const
		{ return M_cdt.number_of_faces(); }
	
	//! Compute the maximum error between the length of a fracture and its approximation
	/*!
	 * @return The pair containing the maximum error and the id of the related fracture
	 */
	std::pair<UInt,Real> maxErrorOnFactureLength() const;
	
	//! Compute the maximum relative error between the length of a fracture and its approximation
	/*!
	 * @return The pair containing the maximum relative error and the id of the related fracture
	 */
	std::pair<UInt,Real> maxRelativeErrorOnFactureLength() const;
	
	//! Add a constrained point to the triangulation
	/*!
	 * @return The Vertex_handle to the constrained point
	 */
	Vertex_handle addConstraint( const Geometry::Point2D & p );
	
	//! Add a constrained segment to the triangulation by its vertexes
	/*!
	 * @return The pair of Vertex_handles to the vertexes of the constrained segment
	 */
	FnConstraint addConstraint
		( const Geometry::Point2D & p1, const Geometry::Point2D & p2 );

	//! Add a constrained segment to the triangulation
	/*!
	 * @return The pair of Vertex_handles to the vertexes of the constrained segment
	 */
	FnConstraint addConstraint( const Geometry::Segment2D & s );
	
	//! Add a fracture network to the triangulation
	/*!
	 * @param fn The fracture network to be added to the triangulation
	 */
	void addFractureNetwork( const Geometry::FractureNetwork2D & fn );
	
	//! Build the triangulation
	/*!
	 * @param alphaMin The minimum amplitude (in degrees) for the elements angles
	 * @param maxElementEdge The maximum length of the elements edges
	 */
	void buildTriangulation( const Real alphaMin=20.705, const Real maxElementEdge=0 );
	
	//! Build the adaptive triangulation
	/*!
	 * @param alphaMin The minimum amplitude (in degrees) for the elements angles
	 * @param hmin The maximum length of the edges of the elements near the fractures
	 * @param hmax The maximum length of the edges of the elements away from the fractures
	 * @param transitionRegion The amplitude of the transition region
	 * @param powerLawExponent The exponent of the power law governing the maximum edges length in the transition region
	 */
	void buildAdaptiveTriangulation( const Real alphaMin = 20.705,
					 const Real hmin = 0,
					 const Real hmax = 0,
					 const Real transitionRegion = 0,
					 const Real powerLawExponent = 1 );
	
	//! Find the approximation of the fracture network in term of the triangulation nodes
	void findApproxFractureNetwork();
	
	//! Export in a single vtk file the nodes representing the fracture network
	/*!
	 * (Use paraview to open the vtk file)
	 * @param filename The name of the vtk file created by this method
	 * @return TRUE -> operation ended correctly
				FALSE -> an error occurred
	 */
	bool exportApproxFractureNetworkUniqueVtk(const std::string & fileName) const;
	
	//! Export in different vtk files the nodes representing each fracture of the network
	/*!
	 * (Use paraview to open the vtk file)
	 * @param filename The prefix to the vtk files created by this method
	 * @return TRUE -> operation ended correctly
				FALSE -> an error occurred
	 */
	bool exportApproxFractureNetworkVtk(const std::string & fileName) const;
	
	//! Export in a vtk file the nodes representing a fracture of the network
	/*!
	 * (Use paraview to open the vtk file)
	 * @param filename The prefix to the vtk files created by this method
	 * @param f The id of the fracture to be exported
	 * @return TRUE -> operation ended correctly
				FALSE -> an error occurred
	 */
	bool exportApproxFractureVtk(const std::string & fileName, const UInt & f) const;

protected:
	//! Compute the length of an approximated fracture
	/*!
	 * @param appF The approximated fracture
	 * @return The length of appF
	 */
	Real approxFractureLength( const ApproxFracture & appF ) const;

protected:

	//! It stores the CGAL constrained Delaunay triangulation object
	CDTplus M_cdt;
	//! It stores the fracture network associated with the triangulation
	Geometry::FractureNetwork2D M_fn;
	//! It stores the handles to the constraint inserted in the triangulation
	std::vector<FnConstraint> M_fnConstraints;
	//! It stores the approximations of the fracture in the network
	std::vector<ApproxFracture> M_approxFractureNetwork;
	
};
 
 
} // namespace Geometry
  
#endif /* TRIANGULATION2D_HPP_ */
