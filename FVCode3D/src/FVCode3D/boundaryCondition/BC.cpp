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

UInt BoundaryConditions::selectBC_onFractureEdge(const Rigid_Mesh::Border_Tip_Edge & edge_it) const
{
	bool isD = false;
	UInt borderId = 0;
        
	// select which BC to apply
	for(auto border_it : edge_it.getBorderIds())
	{
		// BC = D > N && the one with greatest id
		if(M_bordersBCMap.at(border_it).getBCType() == Dirichlet)
		{	
			if(!isD)
			{
				isD = true;
				borderId = border_it;
			}
			else
				borderId = (border_it > borderId) ? border_it : borderId;
		}
            else if(!isD && M_bordersBCMap.at(border_it).getBCType() == Neumann)
            {	
                borderId = (border_it > borderId) ? border_it : borderId;
			}
     }
     return borderId;
}

} // namespace FVCode3D
