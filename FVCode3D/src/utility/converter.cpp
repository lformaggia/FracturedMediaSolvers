/*!
 *	@file converter.cpp
 *	@brief Methods to convert format files (definitions).
 */

#include "mesh/Mesh3D.hpp"
#include "mesh/Properties.hpp"
#include "utility/converter.hpp"

void saveAsSolverFormat(const std::string filename, Geometry::Mesh3D & mesh, Geometry::PropertiesMap & properties)
{
	std::fstream file;

	file.open (filename.c_str(), std::ios_base::out);

	if (file.is_open())
	{
		std::cout << std::endl << " File: " << filename << ", successfully opened";
	}
	else
	{
		std::cerr << std::endl << " *** Error: file not opened *** " << std::endl << std::endl;
		return;
	}

	std::cout << std::endl << " Converting from TPFA to solver format... " << std::endl;

	UInt nNodes, nFacets, nCells, nFractures;
	UInt nodesFacet, facetsCell, facetsFracture;
	UInt i, j;

	std::vector<Geometry::Point3D> & nodesRef = mesh.getNodesVector();
	std::map<UInt, Geometry::Mesh3D::Facet3D> & facetsRef = mesh.getFacetsMap();
	std::map<UInt, Geometry::Mesh3D::Cell3D> & cellsRef = mesh.getCellsMap();
	Geometry::FractureNetwork3D & FN = mesh.getFn();

	nNodes = nodesRef.size();
	nFacets = facetsRef.size();
	nCells = cellsRef.size();
	nFractures = FN.size();

	file << std::scientific << std::setprecision(0);

	// Save header
	file << "# Mesh3D SolverDataFile" << std::endl << std::endl;

	// Save nodes
	file << "POINTS" << std::endl;
	file << nNodes << std::endl;
	file << std::scientific << std::setprecision(10);
	for(i=0; i < nNodes; ++i)
	{
		file << nodesRef[i].x() << " ";
		file << nodesRef[i].y() << " ";
		file << nodesRef[i].z() << std::endl;
	}
	file << std::endl;
	file << std::scientific << std::setprecision(0);

	// Save facets
	file << "FACETS" << std::endl;
	file << nFacets << std::endl;
	for(i=0; i < nFacets; ++i)
	{
		file << std::scientific << std::setprecision(0);
		nodesFacet = facetsRef[i].getNumberOfPoints();
		file << nodesFacet << " ";

		for(j=0; j < nodesFacet; ++j)
			file << facetsRef[i].getIdVertex(j) << " ";

		file << facetsRef[i].getSeparatedCells().size() << " ";
		for(std::set<UInt>::const_iterator it = facetsRef[i].getSeparatedCells().begin(); it != facetsRef[i].getSeparatedCells().end(); ++it)
			file << *it << " ";

		file << facetsRef[i].getBorderId() << " ";
		file << facetsRef[i].isFracture() << " ";
		file << std::scientific << std::setprecision(10);
		file << properties.getProperties(facetsRef[i].getZoneCode()).M_aperture << " ";
		file << properties.getProperties(facetsRef[i].getZoneCode()).M_porosity << " ";
		file << properties.getProperties(facetsRef[i].getZoneCode()).M_permeability << " ";
		for(j=0; j < 4; ++j)
			file << "0" << " ";
		file << "0" << std::endl;
	}
	file << std::endl;
	file << std::scientific << std::setprecision(0);

	// Save cells
	file << "CELLS" << std::endl;
	file << nCells << std::endl;
	for(i=0; i < nCells; ++i)
	{
		file << std::scientific << std::setprecision(0);
		facetsCell = cellsRef[i].getNumberOfFacets();
		file << facetsCell << " ";

		for(std::set<UInt>::const_iterator it = cellsRef[i].getFacetsSet().begin(); it != cellsRef[i].getFacetsSet().end(); ++it)
			file << *it << " ";

		file << "0" << " ";
		file << std::scientific << std::setprecision(10);
		file << properties.getProperties(cellsRef[i].getZoneCode()).M_porosity << " ";
		file << properties.getProperties(cellsRef[i].getZoneCode()).M_permeability << " ";
		for(j=0; j < 4; ++j)
			file << "0" << " ";
		file << "0" << std::endl;
	}
	file << std::endl;
	file << std::scientific << std::setprecision(0);

	// Save fracture network
	file << "FRACTURE_NETWORK" << std::endl;
	file << nFractures << std::endl;
	for(i=0; i < nFractures; ++i)
	{
		facetsFracture = FN.getFracture(i).getNumberOfFractureFacets();
		file << facetsFracture << " ";

		for(j=0; j < facetsFracture-1; ++j)
			file << FN.getFracture(i).getFractureFacetsId()[j] << " ";
		file << FN.getFracture(i).getFractureFacetsId()[facetsFracture-1] << std::endl;
	}

	file.close();
}
