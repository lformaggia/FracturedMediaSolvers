/*!
 *	@file Quadrature.cpp
 *	@brief Class for integrating functions on a Rigid_Mesh (definitions).
 */

#include "mesh/Rigid_Mesh.hpp"
#include "mesh/Properties.hpp"
#include "quadrature/Quadrature.hpp"

namespace Darcy
{

Quadrature::Quadrature (const Geometry::Rigid_Mesh & rigid_mesh):
	 M_mesh (rigid_mesh), M_size (rigid_mesh.getCellsVector().size()+ rigid_mesh.getFractureFacetsIdsVector().size()),
	 M_quadrature(std::move(std::unique_ptr< QuadratureRule >(new CentroidQuadrature))),
	 M_fractureQuadrature(std::move(std::unique_ptr< QuadratureRule >(new CentroidQuadrature))),
	 M_properties(rigid_mesh.getPropertiesMap())
{}

Quadrature::Quadrature (const Geometry::Rigid_Mesh & rigid_mesh, const QuadratureRule & quadrature):
	M_mesh (rigid_mesh), M_size (rigid_mesh.getCellsVector().size()+ rigid_mesh.getFractureFacetsIdsVector().size()),
	M_quadrature(std::move(quadrature.clone())),
	M_fractureQuadrature(std::move(std::unique_ptr<QuadratureRule >(new CentroidQuadrature))),
	M_properties(rigid_mesh.getPropertiesMap())
{}

Quadrature::Quadrature (const Geometry::Rigid_Mesh & rigid_mesh, const QuadratureRule & quadrature, const QuadratureRule & fracturequadrature):
	M_mesh (rigid_mesh), M_size (rigid_mesh.getCellsVector().size()+ rigid_mesh.getFractureFacetsIdsVector().size()),
	M_quadrature(std::move(quadrature.clone())),
	M_fractureQuadrature(std::move(fracturequadrature.clone())),
	M_properties(rigid_mesh.getPropertiesMap()){}

Real Quadrature::Integrate(const Vector & Integrand)
{
	UInt IntSize = Integrand.size();
	if(IntSize != M_size)
		std::cerr << "ERROR: DIMENSION OF INTEGRAND FUNCTION DIFFERS FROM DIMENSION OF MESH" << std::endl;
	Real integral = 0;
	for (auto cell_it : M_mesh.getCellsVector())
	{
		integral += cell_it.getVolume()*Integrand(cell_it.getId());
	}

	for (auto facet_it : M_mesh.getFractureFacetsIdsVector())
	{
		Real _volume = M_properties.getProperties(facet_it.getZoneCode()).M_aperture * facet_it.getFacet().size();
		integral += _volume*Integrand(facet_it.getIdasCell());
		integral -= _volume/2.*Integrand(facet_it.getSeparated()[0]);
		integral -= _volume/2.*Integrand(facet_it.getSeparated()[1]);
	}

	return integral;
}

Real Quadrature::L2Norm(const Vector& Integrand)
{
	Real integral = 0;
	for (auto cell_it : M_mesh.getCellsVector())
	{
		integral += cell_it.getVolume()*Integrand(cell_it.getId())*Integrand(cell_it.getId());
	}

	for (auto facet_it : M_mesh.getFractureFacetsIdsVector())
	{
		Real _volume = M_properties.getProperties(facet_it.getZoneCode()).M_aperture * facet_it.getFacet().size();
		integral += _volume*Integrand(facet_it.getIdasCell())*Integrand(facet_it.getIdasCell());
		integral -= _volume/2.*Integrand(facet_it.getSeparated()[0])*Integrand(facet_it.getSeparated()[0]);
		integral -= _volume/2.*Integrand(facet_it.getSeparated()[1])*Integrand(facet_it.getSeparated()[1]);
	}

	return sqrt(integral);
}

Vector Quadrature::CellIntegrate (const std::function<Real(Generic_Point)>& func)
{
	UInt N = M_mesh.getCellsVector().size() + M_mesh.getFractureFacetsIdsVector().size();
	UInt counter = 0;
	Vector result(N);
	Real partRes;

	UInt NeighboursId;
	std::vector<Generic_Point> v_Nodes;

	for (auto cell_it : M_mesh.getCellsVector())
	{
		for (auto vertex_it : cell_it.getVertexesIds())
			v_Nodes.push_back (M_mesh.getNodesVector()[vertex_it]);

		result (counter) = M_quadrature->apply(v_Nodes, cell_it.getVolume(), func);

		v_Nodes.clear();
		++counter;
	}

	for (auto facet_it : M_mesh.getFractureFacetsIdsVector())
	{
		// contribution from the fracture
		Real _volume = M_properties.getProperties(facet_it.getZoneCode()).M_aperture * facet_it.getFacet().size();

		for (auto vertex_it : M_mesh.getFacetsVector()[facet_it.getId()].getVertexesIds())
			v_Nodes.push_back (M_mesh.getNodesVector()[vertex_it]);

		result (counter) = M_fractureQuadrature->apply(v_Nodes, _volume, func);
		v_Nodes.clear();

		// remove the contribution of the fracture from the first separated cell
		NeighboursId = M_mesh.getFacetsVector()[facet_it.getFacetId()].getSeparatedCellsIds()[0];

		for (auto vertex_it : M_mesh.getCellsVector()[NeighboursId].getVertexesIds())
			v_Nodes.push_back (M_mesh.getNodesVector()[vertex_it]);

		partRes = M_quadrature->apply(v_Nodes, M_mesh.getCellsVector()[NeighboursId].getVolume(), func);
		v_Nodes.clear();

		result (NeighboursId) -= _volume/2. * partRes / M_mesh.getCellsVector()[NeighboursId].getVolume();

		// remove the contribution of the fracture from the second separated cell
		NeighboursId = M_mesh.getFacetsVector()[facet_it.getFacetId()].getSeparatedCellsIds()[1];

		for (auto vertex_it : M_mesh.getCellsVector()[NeighboursId].getVertexesIds())
			v_Nodes.push_back (M_mesh.getNodesVector()[vertex_it]);

		partRes = M_quadrature->apply(v_Nodes, M_mesh.getCellsVector()[NeighboursId].getVolume(), func);
		v_Nodes.clear();

		result (NeighboursId) -= _volume/2. * partRes / M_mesh.getCellsVector()[NeighboursId].getVolume();

		++counter;
	}

	return result;
}

Real Quadrature::Integrate(const std::function<Real(Generic_Point)>& Integrand)
{
	UInt N = M_mesh.getCellsVector().size() + M_mesh.getFractureFacetsIdsVector().size();
	Vector result(N);
	Real integral = 0;
	result = CellIntegrate(Integrand);

	for (UInt it = 0; it < N; ++it)
		integral += result (it);
	return integral;
}

} // namespace Darcy
