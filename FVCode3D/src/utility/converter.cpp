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

	std::cout << std::endl << " Converting TPFA Standard in TPFA with BC... " << std::endl;

	UInt nNodes, nFacets, nCells, nZones, nFractures, volumeCorrection=1;
	UInt nodesFacet, facetsCell, facetsFracture;
	UInt i, j;

	std::vector<Geometry::Point3D> & nodesRef = mesh.getNodesVector();
	std::map<UInt, Geometry::Mesh3D::Facet3D> & facetsRef = mesh.getFacetsMap();
	std::map<UInt, Geometry::Mesh3D::Cell3D> & cellsRef = mesh.getCellsMap();
	Geometry::FractureNetwork3D & FN = mesh.getFn();

	nNodes = nodesRef.size();
	nFacets = facetsRef.size();
	nCells = cellsRef.size();
	nZones = properties.getNumberOfZone();
	nFractures = FN.size();

	// Save header
	file << nNodes << " ";
	file << nFacets << " ";
	file << nCells << " ";
	file << nZones << " ";
	file << nFractures << " ";
	file << volumeCorrection << std::endl;

	// Save nodes
	for(i=0; i < nNodes; ++i)
	{
		file << nodesRef[i].x() << " ";
		file << nodesRef[i].y() << " ";
		file << nodesRef[i].z() << std::endl;
	}

	// Save facets
	for(i=0; i < nFacets; ++i)
	{
		nodesFacet = facetsRef[i].getNumberOfPoints();
		file << nodesFacet << " ";

		for(UInt j=0; j < nodesFacet; ++j)
			file << facetsRef[i].getIdVertex(j) << " ";

		file << facetsRef[i].getZoneCode() << " ";
		file << facetsRef[i].getBorderId() << std::endl;
	}

	// Save cells
	for(i=0; i < nCells; ++i)
	{
		facetsCell = cellsRef[i].getNumberOfFacets();
		file << facetsCell << " ";

		for(std::set<UInt>::const_iterator it = cellsRef[i].getFacetsSet().begin(); it != cellsRef[i].getFacetsSet().end(); ++it)
			file << *it << " ";

		file << cellsRef[i].getZoneCode() << std::endl;
	}

	// Save properties
	for(std::map<UInt,Geometry::Properties>::const_iterator it = properties.getProperties().begin(); it != properties.getProperties().end(); ++it)
	{
		file << it->first << " ";
		file << it->second.M_aperture << " ";
		file << it->second.M_porosity << " ";
		file << it->second.M_permeability << std::endl;
	}

	// Save fracture network
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
