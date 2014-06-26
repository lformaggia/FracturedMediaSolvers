/*!
 *	@file BC.cpp
 *	@brief This class handles the boundary conditions of the Darcy problem.
 */

#include <FVCode3D/mesh/Mesh3D.hpp>
#include <FVCode3D/boundaryCondition/BC.hpp>

namespace FVCode3D
{

void BoundaryConditions::setBoundaryConditions(std::vector<BorderBC> & borderBC)
{
	for(std::vector<BorderBC>::const_iterator it = borderBC.begin(); it != borderBC.end(); ++it)
		M_bordersBCMap.emplace(std::piecewise_construct, std::forward_as_tuple(it->getId()), std::forward_as_tuple(*it));

	for(std::map<UInt,BorderBC>::iterator it = M_bordersBCMap.begin(); it != M_bordersBCMap.end(); ++it)
		it->second.M_bcContainer = this;
}

} // namespace FVCode3D
