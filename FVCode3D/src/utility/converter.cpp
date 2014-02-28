/*!
 *	@file converter.cpp
 *	@brief Methods to convert format files (definitions).
 */

#include "mesh/Mesh3D.hpp"
#include "mesh/Properties.hpp"
#include "utility/converter.hpp"

void readTPFAStandardAsTPFAWithBC(Geometry::Mesh3D & mesh, Geometry::PropertiesMap & properties, const Real theta)
{
	Geometry::FractureNetwork3D FN(mesh);
	std::vector<Geometry::Fracture3D> fracturesVector;
	std::map<UInt, Geometry::Fracture3D> fracturesMap;
	std::map<UInt, Geometry::Fracture3D>::iterator itF;
	Geometry::Point3D normal, center, centerFace;
	Real max;
	UInt compMax;

	extractBC(mesh,properties, theta);

	for(std::map<UInt, Geometry::Mesh3D::Facet3D>::iterator it = mesh.getFacetsMap().begin(); it != mesh.getFacetsMap().end(); ++it)
	{
		if (it->second.getZoneCode() > 0 && it->second.getBorderId()==0)
		{
			itF = fracturesMap.find(it->second.getZoneCode());
			if (itF != fracturesMap.end())
				itF->second.push_back(it->first);
			else
			{
				fracturesMap.emplace(std::piecewise_construct, std::forward_as_tuple(it->second.getZoneCode()), std::forward_as_tuple(mesh) );
				fracturesMap.at(it->second.getZoneCode()).push_back(it->first);
				fracturesMap.at(it->second.getZoneCode()).getId() = it->second.getZoneCode();
			}
		}
	}

	fracturesVector.reserve(fracturesMap.size());
	for(itF = fracturesMap.begin();  itF != fracturesMap.end(); ++itF)
		fracturesVector.push_back(itF->second);

	FN.addFractures(fracturesVector);

	mesh.addFractureNetwork(FN);
}

void extractBC(Geometry::Mesh3D & mesh, Geometry::PropertiesMap & properties, const Real theta)
{
	Geometry::Point3D normal, center, centerFace;
	Real max;
	UInt compMax;

	for(std::map<UInt, Geometry::Mesh3D::Facet3D>::iterator it = mesh.getFacetsMap().begin(); it != mesh.getFacetsMap().end(); ++it)
	{
		if (it->second.getSeparatedCells().size() == 1)
		{
			mesh.getCellsMap().at(*(it->second.getSeparatedCells().begin())).computeCentroid();
			center = mesh.getCellsMap().at(*(it->second.getSeparatedCells().begin())).getCentroid();

			it->second.computeCentroid();
			centerFace = it->second.getCentroid();
			normal = it->second.computeNormal();
			center -= centerFace;
			center.normalize();

			if(normal * center > 0.)
				normal = -normal;

			normal = Geometry::rotateOf(normal, theta);

			max = std::max(std::fabs(normal.x()), std::fabs(normal.y()));
			compMax = std::fabs(normal.x()) > std::fabs(normal.y()) ? 0 : 1;

			max = std::max(max, std::fabs(normal.z()));
			compMax = max > std::fabs(normal.z()) ? compMax : 2;

			if(normal[compMax]>0.)
				it->second.setBorderID(2*compMax+1);
			else
				it->second.setBorderID(2*compMax+2);
		}
	}
}

void convertFromTPFAStandardToTPFAWithBC(const std::string filename, Geometry::Mesh3D & mesh, Geometry::PropertiesMap & properties)
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
	Real aperture, porosity, permeability;

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
