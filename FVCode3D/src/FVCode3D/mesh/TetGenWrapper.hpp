/*!
 *  @file TetGenWrapper.hpp
 *  @brief Class for tetrahedralizing a polyhedron.
 */

#ifndef TETGENWRAPPER_HPP_
#define TETGENWRAPPER_HPP_

#include <FVCode3D/core/TypeDefinition.hpp>

namespace FVCode3D
{

//! Class that tetrahedralizes a polyhedron and computes its volume
/*!
 * @class TetGenWrapper
 * This class generates a tetrahedral mesh moving from a set of nodes and faces that describe a polyhedron.
 * It permits to compute the volume of the polyhedron.
 */
class TetGenWrapper
{
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

    //! Compute the center of mass of the polyhedron by considering a constant density
    /*!
     * C = 1/|Omega| int_\Omega x dx
     * @return the center of mass
     * @pre call generateMesh()
     */
    Point3D computeCenterOfMass();

    //! Get the volume of the mesh
    /*!
     * @return the volume of the polyhedron
     * @pre call computeVolume()
     */
    Real getVolume() const
        { return M_volume; }

    //! Get the center of mass of the polyhedron
    /*!
     * @return the center of mass of the polyhedron
     * @pre call computeCenterOfMass()
     */
    Real getCenterOfMass() const
        { return M_volume; }

    //! Get i-th element
    /*!
     * @param i i-th element
     * @return vector that contains the 4 points that define the tetrahedron
     */
    const std::vector<Point3D> getElement(const UInt i) const;

    //! Destructor
    ~TetGenWrapper() = default;

private:

    //! No default constructor
    TetGenWrapper() = delete;

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
    //! Center of mass
    Point3D M_centerOfMass;
};

} // namespace FVCode3D

#endif /* TETGENWRAPPER_HPP_ */
