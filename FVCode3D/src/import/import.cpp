/*!
 *	@file import.cpp
 *	@brief Classes for loading files (definitions).
 */

#include "mesh/Mesh3D.hpp"
#include "mesh/Properties.hpp"
#include "import/import.hpp"

#include <utility>
#include <map>

void ImporterTPFAStandard::import()
{
	std::ifstream file;
	file.open(M_filename.c_str(), std::ios_base::in);
	if(!file)
	{
		std::cerr << "File not opened!" << std::endl;
		exit(0);
	}

	UInt nNodes, nFacets, nCells, nZones, nFractures, volumeCorrection;
	UInt nodesFacet, facetsCell, facetsFracture, bcId=0;
	Int zone;
	UInt i, j;
	Real x, y, z;
	Real aperture, porosity, permeability;
	std::vector<UInt> tmp;

	std::vector<Geometry::Point3D> & nodesRef = M_mesh.getNodesVector();
	std::map<UInt, Geometry::Mesh3D::Facet3D> & facetsRef = M_mesh.getFacetsMap();
	std::map<UInt, Geometry::Mesh3D::Cell3D> & cellsRef = M_mesh.getCellsMap();
	Geometry::FractureNetwork3D & FN = M_mesh.getFn();

	Geometry::Properties prop;

	// Read header
	file >> nNodes;
	file >> nFacets;
	file >> nCells;
	file >> nZones;
	nFractures = 0;
	file >> volumeCorrection;

	// Read nodes
	nodesRef.resize(nNodes);

	for(i=0; i < nNodes; ++i)
	{
		file >> x; file >> y; file >> z;
		nodesRef[i] = Geometry::Point3D(x,y,z);
	}

	// Read facets
	for(i=0; i < nFacets; ++i)
	{
		file >> nodesFacet;
		tmp.resize(nodesFacet);

		for(UInt j=0; j < nodesFacet; ++j)
			file >> tmp[j];
		file >> zone;
		//facetsRef.insert( std::pair<UInt, Geometry::Mesh3D::Facet3D>( i, Geometry::Mesh3D::Facet3D(&M_mesh, tmp, zone+1 , 0)) );
		facetsRef.emplace( std::piecewise_construct, std::forward_as_tuple(i), std::forward_as_tuple(&M_mesh, tmp, zone+1, bcId) );
	}

	// Read cells
	for(i=0; i < nCells; ++i)
	{
		file >> facetsCell;
		tmp.resize(facetsCell);

		for(UInt j=0; j < facetsCell; ++j)
			file >> tmp[j];
		file >> zone;
		//cellsRef.insert( std::pair<UInt,Geometry::Mesh3D::Cell3D>( i, Geometry::Mesh3D::Cell3D(&M_mesh, tmp, zone+1 )) );
		cellsRef.emplace( std::piecewise_construct, std::forward_as_tuple(i), std::forward_as_tuple(&M_mesh, tmp, zone+1) );
	}
	tmp.clear();

	prop.setProperties(1., 0.25, 100);
	M_properties.setZone(0, prop);

	// Read properties
	for(i=0; i < nZones; ++i)
	{
		file >> zone;
		file >> aperture;
		file >> porosity;
		file >> permeability; // dummy parameter
		file >> permeability;
		prop.setProperties(aperture, porosity, permeability);
		M_properties.setZone(zone+1, prop);
	}

	// Read fracture network
	for(i=0; i < nFractures; ++i)
	{
		file >> facetsFracture;
		tmp.resize(facetsFracture);
		for(j=0; j < facetsFracture; ++j)
			file >> tmp[j];
		FN.emplace_back(Geometry::Fracture3D(M_mesh, tmp, i));
	}

	file.close();
}

void ImporterTPFAWithBC::import()
{

	std::ifstream file;
	file.open(M_filename.c_str(), std::ios_base::in);
	if(!file)
	{
		std::cerr << "File not opened!" << std::endl;
		exit(0);
	}

	UInt nNodes, nFacets, nCells, nZones, nFractures, volumeCorrection;
	UInt nodesFacet, facetsCell, facetsFracture, zone, bcId;
	UInt i, j;
	Real x, y, z;
	Real aperture, porosity, permeability;
	std::vector<UInt> tmp;

	std::vector<Geometry::Point3D> & nodesRef = M_mesh.getNodesVector();
	std::map<UInt, Geometry::Mesh3D::Facet3D> & facetsRef = M_mesh.getFacetsMap();
	std::map<UInt, Geometry::Mesh3D::Cell3D> & cellsRef = M_mesh.getCellsMap();
	Geometry::FractureNetwork3D & FN = M_mesh.getFn();

	Geometry::Properties prop;

	// Read header
	file >> nNodes;
	file >> nFacets;
	file >> nCells;
	file >> nZones;
	file >> nFractures;
	file >> volumeCorrection;

	// Read nodes
	nodesRef.resize(nNodes);

	for(i=0; i < nNodes; ++i)
	{
		file >> x; file >> y; file >> z;
		nodesRef[i] = Geometry::Point3D(x,y,z);
	}

	// Read facets
	for(i=0; i < nFacets; ++i)
	{
		file >> nodesFacet;
		tmp.resize(nodesFacet);

		for(UInt j=0; j < nodesFacet; ++j)
			file >> tmp[j];
		file >> zone;
		file >> bcId;
		//facetsRef.insert( std::pair<UInt, Geometry::Mesh3D::Facet3D>( i, Geometry::Mesh3D::Facet3D(&M_mesh, tmp, zone, bcId)) );
		facetsRef.emplace( std::piecewise_construct, std::forward_as_tuple(i), std::forward_as_tuple(&M_mesh, tmp, zone, bcId) );
	}

	// Read cells
	for(i=0; i < nCells; ++i)
	{
		file >> facetsCell;
		tmp.resize(facetsCell);

		for(UInt j=0; j < facetsCell; ++j)
			file >> tmp[j];
		file >> zone;

		//cellsRef.insert( std::move(std::pair<UInt,Geometry::Mesh3D::Cell3D>( i, Geometry::Mesh3D::Cell3D(&M_mesh, tmp, zone))) );
		cellsRef.emplace( std::piecewise_construct, std::forward_as_tuple(i), std::forward_as_tuple(&M_mesh, tmp, zone) );
	}
	tmp.clear();

	// Read properties
	for(i=0; i < nZones; ++i)
	{
		file >> zone;
		file >> aperture;
		file >> porosity;
		file >> permeability;
		prop.setProperties(aperture, porosity, permeability);
		M_properties.setZone(zone, prop);
	}

	// Read fracture network
	for(i=0; i < nFractures; ++i)
	{
		file >> facetsFracture;
		tmp.resize(facetsFracture);
		for(j=0; j < facetsFracture; ++j)
			file >> tmp[j];
		FN.emplace_back(Geometry::Fracture3D(M_mesh, tmp, i));
	}

	file.close();

}
