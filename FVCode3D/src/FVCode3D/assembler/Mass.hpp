/*!
 *  @file mass.hpp
 *  @brief This class build a Mass-matrix for a finite volume method.
 */

#ifndef __DARCYMASS_HPP__
#define __DARCYMASS_HPP__

#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>
#include <FVCode3D/assembler/MatrixHandler.hpp>

namespace FVCode3D
{

class PropertiesMap;

//! Class for assembling a mass matrix
/*!
 * @class MassMatrix
 * This class constructs the mass-matrix.
 * The adopted technique is the one of the finite volume method:
 * it hence represents the volume of the cell. The fractures are considered as cells.
 */
class MassMatrix: public MatrixHandler
{
public:
    //! @name Constructor & Destructor
    //@{

    //! Construct a Mass-Matrix, given a Rigid_Mesh.
    /*!
        @param rigid_mesh A Rigid_Mesh used to build the matrix
    */
    MassMatrix(const Rigid_Mesh & rigid_mesh):
        MatrixHandler(rigid_mesh) {}
    //! No Copy-Constructor
    MassMatrix(const MassMatrix&) = delete;
    //! No Empty-Constructor
    MassMatrix() = delete;
    //! Destructor
    ~MassMatrix() = default;
    //@}

    //! @name Methods
    //@{

    //! Assemble method
    /*!
     * Assemble the Mass matrix
     */
    void assemble();
    //@}

}; // class MassMatrix

} //namespace FVCode3D

#endif
