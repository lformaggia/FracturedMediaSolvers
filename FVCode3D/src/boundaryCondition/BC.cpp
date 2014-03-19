/*!
 *	@file BC.cpp
 *	@brief This class handles the boundary conditions of the Darcy problem.
 */

#include "mesh/Mesh3D.hpp"
#include "boundaryCondition/BC.hpp"

BoundaryConditions::BoundaryConditions(std::vector<BorderBC> & borderbc)
{
	for(std::vector<BorderBC>::const_iterator it = borderbc.begin(); it != borderbc.end(); ++it)
		BordersBCMap.emplace(std::piecewise_construct, std::forward_as_tuple(it->getId()), std::forward_as_tuple(*it));

	for(std::map<UInt,BorderBC>::iterator it = BordersBCMap.begin(); it != BordersBCMap.end(); ++it)
		it->second.m_bcContainer = this;
}
