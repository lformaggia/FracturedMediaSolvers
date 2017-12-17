/*!
 * @file SaddlePoint.cpp
 * @brief Classe that assemble the saddle point blocks.
 */

#include <cmath>
#include <FVCode3D/core/TypeDefinition.hpp>
#include <FVCode3D/mesh/RigidMesh.hpp>
#include <FVCode3D/assembler/global_operator.hpp>
#include <FVCode3D/assembler/SaddlePoint.hpp>

namespace FVCode3D
{
	
void SaddlePoint_StiffMatHandler::assemble()
{
	//Define the block matrices
	auto & M = getM();
	auto & B = getB();
	auto & T = getT();
	
	// Define the degrees of fredoom
	const UInt numFracture  = M_mesh.getFractureFacetsIdsVector().size();
	const UInt numFacetsTot = M_mesh.getFacetsVector().size() + numFracture;
	const UInt numCell      = M_mesh.getCellsVector().size();
	
	// Define the global inner product
	global_InnerProduct gIP(M_mesh, dType::dFacet, dType::dFacet, numFacetsTot, numFacetsTot);
	gIP.ShowMe();
	
	// Define the global divergence operator
	global_Div gDIV(M_mesh, dType::dCell, dType::dFacet, numCell, numFacetsTot);
	gDIV.ShowMe();
	
	// Define the coupling conditions
	CouplingConditions coupling(M_mesh, dType::dFracture, dType::dFacet, numFracture, numFacetsTot, gIP.getMatrix());
	coupling.ShowMe();
	
	// Define the flux operator
	FluxOperator fluxOP(M_mesh, dType::dFracture, dType::dFracture, numFracture, numFracture);
	fluxOP.ShowMe();
	
	// Define the global bulk builder
	global_BulkBuilder gBulkBuilder(M_mesh, gIP, gDIV);
	// Reserve space for the bulk matrices
	gBulkBuilder.reserve_space(M, B);
	// Build the bulk matrices
	gBulkBuilder.build(M, B);
	
	// Define the fracture builder
	FractureBuilder FBuilder(M_mesh, coupling, fluxOP);
	// Reserve space for the fracture matrices
	FBuilder.reserve_space(M, B, T);
	// Build the fracture matrices
	FBuilder.build(M, B, T);

	// Define the BCs imposition object
	BCimposition BCimp(M_mesh, M_bc);
	// Impose BCs on bulk
	BCimp.ImposeBConBulk(M, B, M_b);
	// Impose BCs in fractures
	BCimp.ImposeBConFracture_onT(T, M_b, fluxOP);
	
	//Comrpess the saddle point matrix
	M_SP.makeCompressed();		
}
	
		
} // FVCode3D
