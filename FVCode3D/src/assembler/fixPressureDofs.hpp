/*!
 * @file fixPressureDofs.hpp
 * @brief Class that fix the pressure in the fractures.
 */

#ifndef FIXPRESSUREDOFS_HPP_
#define FIXPRESSUREDOFS_HPP_

#include "core/TypeDefinition.hpp"

template <typename ProblemType>
class FixPressureDofs
{
public:

    //! No empty constructor
    FixPressureDofs() = delete;

    //! Constructor
    /*!
     * @param problem problem to modify
     */
    FixPressureDofs(ProblemType * problem):
        M_problem(problem) {}

    //! Get problem
    /*!
     * @return the problem
     */
    ProblemType & getProblem()
        { return *M_problem; }

    //! Set the pressure dofs
    /*!
     * @param pressure pressure value to set inside the fractures
     */
    void apply(const Real pressure);

    //! Default destructor
    ~FixPressureDofs() {};

private:

    //! Problem
    ProblemType * M_problem;

}; // class MSR

template <typename ProblemType>
void FixPressureDofs<ProblemType>::apply(const Real pressure)
{
    SpMat & A = M_problem->getA();
    Vector & b = M_problem->getb();

    const UInt dofTotal = A.outerSize();
    const UInt dofFractures = M_problem->getMesh().getFractureFacetsIdsVector().size();
    const UInt dofPorousMatrix = dofTotal - dofFractures;

    // Clear the fractures rows
    A.bottomRows( dofFractures ) *= 0.;

    const Real weight = A.diagonal().sum() / dofPorousMatrix;

    b.tail( dofFractures ).setConstant( weight * pressure );

    for(UInt index = dofPorousMatrix; index < dofTotal; ++index )
        A.coeffRef( index, index ) = weight;

    A.prune(0.);
}

#endif /* FIXPRESSUREDOFS_HPP_ */
