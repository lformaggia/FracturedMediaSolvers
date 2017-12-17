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
 * @class MassMatHandlerFV
 * This class constructs the mass-matrix.
 * The adopted technique is the one of the finite volume method:
 * it hence represents the volume of the cell. The fractures are considered as cells.
 */
class MassMatHandlerFV: public MatrixHandlerFV
{
public:
    //! @name Constructor & Destructor
    //@{

    //! Construct a Mass-Matrix, given a Rigid_Mesh.
    /*!
        @param rigid_mesh A Rigid_Mesh used to build the matrix
        @param Mat The matrix
    */
    MassMatHandlerFV(const Rigid_Mesh & rigid_mesh, SpMat & Mat):
        MatrixHandlerFV(rigid_mesh, Mat) {}
    //! No Copy-Constructor
    MassMatHandlerFV(const MassMatHandlerFV &) = delete;
    //! No Empty-Constructor
    MassMatHandlerFV() = delete;
    //! Destructor
    ~MassMatHandlerFV() = default;
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
