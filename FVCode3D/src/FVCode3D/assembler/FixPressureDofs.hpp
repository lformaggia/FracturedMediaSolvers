/*!
 * @file fixPressureDofs.hpp
 * @brief Class that fix the pressure in the fractures.
 */

#ifndef FIXPRESSUREDOFS_HPP_
#define FIXPRESSUREDOFS_HPP_

#include <FVCode3D/core/TypeDefinition.hpp>

namespace FVCode3D
{

//! Class for fixing a pressure value inside the fractures
/*!
 * @class FixPressureDofs
 * This class replaces each row of the algebraic matrix associated with
 * the fractures with a row of zeros, but on the diagonal, where a constant value @c is applied.
 * The same dofs on the right and side are replaced by @c multiplied by the pressure imposed on the fractures.
 * This procedure is done only if the system is then solved through a direct method and so the system matrix is 
 * stored as a simple SpMat.
 */
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
    void apply(const Real pressure) throw();

    //! Default destructor
    ~FixPressureDofs() = default;

private:

    //! Problem
    ProblemType * M_problem;

}; // class FixPressureDofs

template <typename ProblemType>
void FixPressureDofs<ProblemType>::apply(const Real pressure) throw()
{
	auto Algebry  = M_problem->getAlgebry();
	auto & A    = std::get<0>(Algebry);
	auto & b    = std::get<2>(Algebry);

	if(A != nullptr)
	{
    const UInt dofTotal = A->rows();
    const UInt dofFractures = M_problem->getMesh().getFractureFacetsIdsVector().size();
    const UInt dofPorousMatrix = dofTotal - dofFractures;

    // Clear the fractures rows
    A->bottomRows( dofFractures ) *= 0.;

    const Real weight = A->diagonal().sum() / dofPorousMatrix;

    b->tail( dofFractures ).setConstant( weight * pressure );

    for(UInt index = dofPorousMatrix; index < dofTotal; ++index )
        A->coeffRef( index, index ) = weight;

    A->prune(0.);
	}
	else
	{
		std::stringstream error;	
		error << "Error: fix pressure only if the solver is Direct"<<std::endl;
		throw std::runtime_error(error.str());	
	}
}

} // namespace FVCode3D
#endif /* FIXPRESSUREDOFS_HPP_ */
