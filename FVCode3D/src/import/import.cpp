	/*!
 *	@file import.cpp
 *	@brief Classes for loading files (definitions).
 */

#include "mesh/Mesh3D.hpp"
#include "mesh/Properties.hpp"
#include "import/import.hpp"

#include <utility>
#include <map>

void Importer::extractBC(const Real theta)
{
	Geometry::Point3D normal, center, centerFace;
	Real max;
	UInt compMax;

	for(std::map<UInt, Geometry::Mesh3D::Facet3D>::iterator it = M_mesh.getFacetsMap().begin(); it != M_mesh.getFacetsMap().end(); ++it)
	{
		if (it->second.getSeparatedCells().size() == 1)
		{
			center = M_mesh.getCellsMap().at(*(it->second.getSeparatedCells().begin())).getCentroid();
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

void Importer::addBCAndFractures(const Real theta)
{
	extractBC(theta);
	addFractures();
}

void ImporterMedit::import(bool fracturesOn)
{
	std::ifstream file;
	file.open(M_filename.c_str(), std::ios_base::in);
	if(!file)
	{
		std::cerr << "File not opened!" << std::endl;
		exit(0);
	}

	UInt nNodes, nFacets, nCells;
	UInt zone, maxZone, bcId;
	UInt i, j;
	Real x, y, z;
	std::vector<UInt> tmp, tmpNodes(3), tmpFacets(4);
	std::set<UInt> zones;
	std::string buffer = "";

	std::vector<Geometry::Point3D> & nodesRef = M_mesh.getNodesVector();
	std::map<UInt, Geometry::Mesh3D::Facet3D> & facetsRef = M_mesh.getFacetsMap();
	std::map<UInt, Geometry::Mesh3D::Cell3D> & cellsRef = M_mesh.getCellsMap();

	Geometry::Properties prop;

	// Read nodes
	const std::string s2findN = "Vertices";
	while(buffer!=s2findN)
		getline(file, buffer);

	file >> nNodes;
	nodesRef.resize(nNodes);

	for(i=0; i < nNodes; ++i)
	{
		file >> x; file >> y; file >> z;
		file >> buffer;
		nodesRef[i] = Geometry::Point3D(x,y,z);
	}
	
	// Read facets
	const std::string s2findT = "Triangles";
	while(buffer!=s2findT)
		getline(file, buffer);
	
	file >> nFacets;

	tmp.resize(3);
	for(i=0; i < nFacets; ++i)
	{
		for(UInt j=0; j < 3; ++j)
		{
			file >> tmp[j];
			tmp[j]--;
		}
		file >> zone;
		bcId = zone <= 1000 ? 1 : 0;
		zone = (zone > 1000) && (zone <= 2000) ? zone : 0;
		facetsRef.emplace( std::piecewise_construct, std::forward_as_tuple(i), std::forward_as_tuple(&M_mesh, tmp, (zone)*static_cast<UInt>(fracturesOn), bcId) );
		if(zone > 1 && fracturesOn && zones.find(zone) == zones.end())
		{
			prop.setProperties(1., 1., 1.);
			M_properties.setZone(zone, prop);
		}
	}
	tmp.clear();

	// Read cells
	const std::string s2findE = "Tetrahedra";
	while(buffer!=s2findE)
		getline(file, buffer);
	
	file >> nCells;
	
	maxZone = *std::max_element(zones.begin(), zones.end());
	prop.setProperties(1., 1., 1.);
	M_properties.setZone(maxZone+1, prop);
	
	tmp.resize(4);
	for(i=0; i < nCells; ++i)
	{		
		// get the nodes that define the cell
		for(UInt j=0; j < 4; ++j)
		{
			file >> tmp[j];
			tmp[j]--;
		}
		file >> buffer;
		
		// costruisco le facce della cella partendo dai nodi
		for(j=0; j < 4; ++j)
		{
			tmpNodes[0] = tmp[j % 4];
			tmpNodes[1] = tmp[(j+1) % 4];
			tmpNodes[2] = tmp[(j+2) % 4];
			tmpFacets[j] = M_mesh.getFacetFromNodes(tmpNodes);
		}
		
		cellsRef.emplace( std::piecewise_construct, std::forward_as_tuple(i), std::forward_as_tuple(&M_mesh, tmpFacets, maxZone+1) );
	}
	tmp.clear();

	file.close();
}

void ImporterMedit::addFractures()
{
	Geometry::FractureNetwork3D FN(M_mesh);
	std::vector<Geometry::Fracture3D> fracturesVector;
	std::map<UInt, Geometry::Fracture3D> fracturesMap;
	std::map<UInt, Geometry::Fracture3D>::iterator itF;
	Geometry::Point3D normal, center, centerFace;

	for(std::map<UInt, Geometry::Mesh3D::Facet3D>::iterator it = M_mesh.getFacetsMap().begin(); it != M_mesh.getFacetsMap().end(); ++it)
	{
		if (it->second.getZoneCode() > 0 && it->second.getBorderId()==0)
		{
			itF = fracturesMap.find(it->second.getZoneCode());
			if (itF != fracturesMap.end())
				itF->second.push_back(it->first);
			else
			{
				fracturesMap.emplace(std::piecewise_construct, std::forward_as_tuple(it->second.getZoneCode()), std::forward_as_tuple(M_mesh) );
				fracturesMap.at(it->second.getZoneCode()).push_back(it->first);
				fracturesMap.at(it->second.getZoneCode()).getId() = it->second.getZoneCode();
			}
		}
	}

	fracturesVector.reserve(fracturesMap.size());
	for(itF = fracturesMap.begin();  itF != fracturesMap.end(); ++itF)
		fracturesVector.push_back(itF->second);

	FN.addFractures(fracturesVector);

	M_mesh.addFractureNetwork(FN);
}

void ImporterTPFA::import(bool fracturesOn)
{
	std::ifstream file;
	file.open(M_filename.c_str(), std::ios_base::in);
	if(!file)
	{
		std::cerr << "File not opened!" << std::endl;
		exit(0);
	}

	UInt nNodes, nFacets, nCells, nZones, volumeCorrection;
	UInt nodesFacet, facetsCell, bcId=0;
	Int zone;
	UInt i, j;
	Real x, y, z;
	Real aperture, porosity, permeability;
	std::vector<UInt> tmp;

	std::vector<Geometry::Point3D> & nodesRef = M_mesh.getNodesVector();
	std::map<UInt, Geometry::Mesh3D::Facet3D> & facetsRef = M_mesh.getFacetsMap();
	std::map<UInt, Geometry::Mesh3D::Cell3D> & cellsRef = M_mesh.getCellsMap();

	Geometry::Properties prop;

	// Read header
	file >> nNodes;
	file >> nFacets;
	file >> nCells;
	file >> nZones;
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

		for(j=0; j < nodesFacet; ++j)
			file >> tmp[j];
		file >> zone;
		facetsRef.emplace( std::piecewise_construct, std::forward_as_tuple(i), std::forward_as_tuple(&M_mesh, tmp, (zone+1)*static_cast<UInt>(fracturesOn), bcId) );
	}

	// Read cells
	for(i=0; i < nCells; ++i)
	{
		file >> facetsCell;
		tmp.resize(facetsCell);

		for(j=0; j < facetsCell; ++j)
			file >> tmp[j];
		file >> zone;
		cellsRef.emplace( std::piecewise_construct, std::forward_as_tuple(i), std::forward_as_tuple(&M_mesh, tmp, zone+1) );
	}
	tmp.clear();

	//prop.setProperties(1., 0.25, 100);
	//M_properties.setZone(0, prop);

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

	file.close();
}

void ImporterTPFA::addFractures()
{
	Geometry::FractureNetwork3D FN(M_mesh);
	std::vector<Geometry::Fracture3D> fracturesVector;
	std::map<UInt, Geometry::Fracture3D> fracturesMap;
	std::map<UInt, Geometry::Fracture3D>::iterator itF;
	Geometry::Point3D normal, center, centerFace;

	for(std::map<UInt, Geometry::Mesh3D::Facet3D>::iterator it = M_mesh.getFacetsMap().begin(); it != M_mesh.getFacetsMap().end(); ++it)
	{
		if (it->second.getZoneCode() > 0 && it->second.getBorderId()==0)
		{
			itF = fracturesMap.find(it->second.getZoneCode());
			if (itF != fracturesMap.end())
				itF->second.push_back(it->first);
			else
			{
				fracturesMap.emplace(std::piecewise_construct, std::forward_as_tuple(it->second.getZoneCode()), std::forward_as_tuple(M_mesh) );
				fracturesMap.at(it->second.getZoneCode()).push_back(it->first);
				fracturesMap.at(it->second.getZoneCode()).getId() = it->second.getZoneCode();
			}
			//M_properties.getProperties(it->second.getZoneCode()).M_porosity = 1;
			//M_properties.getProperties(it->second.getZoneCode()).M_permeability = 1000;
			//M_properties.getProperties(it->second.getZoneCode()).M_aperture = 1;
		}
	}

	fracturesVector.reserve(fracturesMap.size());
	for(itF = fracturesMap.begin();  itF != fracturesMap.end(); ++itF)
		fracturesVector.push_back(itF->second);

	FN.addFractures(fracturesVector);

	M_mesh.addFractureNetwork(FN);
}

void ImporterForSolver::import(bool fracturesOn)
{
	std::ifstream file;
	file.open(M_filename.c_str(), std::ios_base::in);
	if(!file)
	{
		std::cerr << "File not opened!" << std::endl;
		exit(0);
	}

	UInt nNodes, nFacets, nCells, nFractures;
	UInt nodesFacet, facetsCell, facetsFracture, nSepCells, bcId, isFrac;
	UInt zone;
	UInt i, j;
	Real x, y, z;
	Real aperture, porosity, permeability;
	std::vector<UInt> tmp;
	std::string buffer = "";
	const UInt permSize = 1;

	std::vector<Geometry::Point3D> & nodesRef = M_mesh.getNodesVector();
	std::map<UInt, Geometry::Mesh3D::Facet3D> & facetsRef = M_mesh.getFacetsMap();
	std::map<UInt, Geometry::Mesh3D::Cell3D> & cellsRef = M_mesh.getCellsMap();
	Geometry::FractureNetwork3D & FN = M_mesh.getFn();
	Geometry::Properties prop;

	// Read nodes
	const std::string s2findN = "POINTS";
	while(buffer!=s2findN)
		file >> buffer;

	file >> nNodes;
	nodesRef.resize(nNodes);

	for(i=0; i < nNodes; ++i)
	{
		file >> x; file >> y; file >> z;
		nodesRef[i] = Geometry::Point3D(x,y,z);
	}

	// Read facets with properties
	const std::string s2findF = "FACETS";
	while(buffer!=s2findF)
		file >> buffer;

	file >> nFacets;

	zone = 1;
	for(i=0; i < nFacets; ++i)
	{
		file >> nodesFacet;
		tmp.resize(nodesFacet);

		for(UInt j=0; j < nodesFacet; ++j)
			file >> tmp[j];
		file >> nSepCells;
		for(UInt j=0; j<nSepCells; ++j)
			file >> buffer;
		file >> bcId;
		file >> isFrac;

		if(isFrac)
		{
			file >> aperture;
			file >> porosity;
			for(UInt j=0; j<6; ++j)
				if (j<permSize)
					file >> permeability;
				else
					file >> buffer;
		}
		facetsRef.emplace( std::piecewise_construct, std::forward_as_tuple(i), std::forward_as_tuple(&M_mesh, tmp, isFrac*static_cast<UInt>(fracturesOn)*zone, bcId) );
		if (isFrac*static_cast<UInt>(fracturesOn))
		{
			prop.setProperties(aperture, porosity, permeability);
			M_properties.setZone(zone, prop);
			zone++;
		}
	}

	// Read cells with properties
	const std::string s2findC = "CELLS";
	while(buffer!=s2findC)
		file >> buffer;

	file >> nCells;

	aperture = 1.;
	for(i=0; i < nCells; ++i)
	{
		file >> facetsCell;
		tmp.resize(facetsCell);

		for(UInt j=0; j < facetsCell; ++j)
			file >> tmp[j];

		file >> porosity;
		for(UInt j=0; j<6; ++j)
			if (j<permSize)
				file >> permeability;
			else
				file >> buffer;
		cellsRef.emplace( std::piecewise_construct, std::forward_as_tuple(i), std::forward_as_tuple(&M_mesh, tmp, zone) );
		prop.setProperties(aperture, porosity, permeability);
		M_properties.setZone(zone, prop);
		zone++;
	}
	tmp.clear();

	// Read fracture network
	const std::string s2findFN= "FRACTURE_NETWORK";
	while(buffer!=s2findFN)
		file >> buffer;

	file >> nFractures;

	for(i=0; i < nFractures; ++i)
	{
		file >> facetsFracture;
		tmp.resize(facetsFracture);
		for(j=0; j < facetsFracture; ++j)
			file >> tmp[j];
		if(facetsFracture>0)
			FN.emplace_back(Geometry::Fracture3D(M_mesh, tmp, i));
		//M_properties.getProperties(it->second.getZoneCode()).M_porosity = 1;
		for(j=0; j<facetsFracture; ++j)
		{
//			M_properties.getProperties(facetsRef[tmp[j]].getZoneCode()).M_permeability = 1000;
//			M_properties.getProperties(facetsRef[tmp[j]].getZoneCode()).M_aperture = 1;
		}
	}

	file.close();
}
