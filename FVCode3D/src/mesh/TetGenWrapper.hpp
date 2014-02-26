/*!
 *	@file TetGenWrapper.hpp
 *	@brief Class for tetrahedralizing a polyhedron.
 */

#ifndef TETGENWRAPPER_HPP_
#define TETGENWRAPPER_HPP_

#include "core/TypeDefinition.hpp"

namespace Geometry{

//! Class that tetrahedralizes a polyhedron and computes its volume
/*!
 * @class TetGenWrapper
 * This class generates a tetrahedral mesh moving from a set of nodes and faces that describe a polyhedron.
 * It permits to compute the volume of the polyhedron.
 */
class TetGenWrapper{
public:

	//! Constructor
	/*!
	 * @param nodes vector of nodes that define the polyhedron
	 * @param faces vector of vector of ids of the nodes. Each vector define a face(polygon) of the polyhedron.
	 */
	TetGenWrapper(std::vector<Point3D> nodes, std::vector< std::vector<UInt> > faces);

	//! Generate a tetrahedralization of the polyhedron
	void generateMesh();

	//! Compute the volume of the polyhedron and return it
	/*!
	 * @return the volume of the polyhedron
	 * @pre call generateMesh()
	 */
	Real computeVolume();

	//! Get the volume of the mesh
	/*!
	 * @return the volume of the polyhedron
	 * @pre call computeVolume()
	 */
	Real getVolume() const
		{ return M_volume; };

	//! Destructor
	~TetGenWrapper() {};

private:

	//! No default constructor
	TetGenWrapper();

	//! Vector of points (input)
	std::vector<Point3D> M_inNodes;
	//! Vector of faces (input)
	std::vector< std::vector<UInt> > M_faces;
	//! Vector of points (output)
	std::vector<Point3D> M_outNodes;
	//! Vector of elements (output)
	std::vector< std::vector<UInt> > M_elements;
	//! Volume of the mesh
	Real M_volume;

};

}// namespace Geometry

#endif /* TETGENWRAPPER_HPP_ */
