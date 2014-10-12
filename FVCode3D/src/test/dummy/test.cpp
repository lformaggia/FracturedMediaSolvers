//#define NDEBUG // add this macro to disable asserts
#include <cassert>
#include <iostream>
#include <chrono>
#include <vector>
#include <fstream>

#include "core/TypeDefinition.hpp"
#include "core/data.hpp"
#include "mesh/Mesh3D.hpp"
#include "property/Properties.hpp"
#include "mesh/cartesianGrid.hpp"
#include "import/import.hpp"
#include "export/exportVTU.hpp"
#include "utility/converter.hpp"
#include "geometry/operations.hpp"

int main(int argc, char * argv[])
{
	GetPot command_line(argc,argv);
	const std::string dataFileName = command_line.follow("data.txt", 2, "-f", "--file");

	std::cout << "Read Data..." << std::flush;
	Data data(dataFileName);
	data.setMeshType(Data::MeshFormatType::Medit);
	data.fractureOn(true);
	data.verbose(true);
	std::cout << " done." << std::endl;

	std::cout << std::endl;
	data.showMe();
	std::cout << std::endl;

	std::cout << "Define Mesh and Properties..." << std::flush;
	Geometry::Mesh3D mesh;
	Geometry::PropertiesMap propMap(data.getMobility());
	std::cout << " done." << std::endl;

	std::ifstream file;
	file.open((data.getMeshDir() + data.getMeshFile()).c_str(), std::ios_base::in);
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

	std::vector<Geometry::Point3D> & nodesRef = mesh.getNodesVector();
	std::map<UInt, Geometry::Mesh3D::Facet3D> & facetsRef = mesh.getFacetsMap();
	std::map<UInt, Geometry::Mesh3D::Cell3D> & cellsRef = mesh.getCellsMap();

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

		Geometry::Point3D centroid;
		UInt N = tmp.size();
		Real tmpArea, totArea = 0.;
		for (UInt j = 1; j < N-1; ++j)
		{
			tmpArea = Geometry::triangleArea(	mesh.getNodesVector()[tmp[0]],
												mesh.getNodesVector()[tmp[j]],
												mesh.getNodesVector()[tmp[j+1]]
											);
			totArea += tmpArea;
			centroid += tmpArea * (	mesh.getNodesVector()[tmp[0]] +
									mesh.getNodesVector()[tmp[j]] +
									mesh.getNodesVector()[tmp[j+1]]
									) / 3.;
		}
		centroid /= totArea;

		if( centroid[0] >= 0.5 && centroid[0] <= 1.5 &&
			centroid[2] > 0.5 - 0.01 && centroid[2] < 0.5 + 0.01)
			zone = 1001;

		bcId = zone <= 1000 ? 1 : 0;
		zone = (zone > 1000) && (zone <= 2000) ? zone : 0;
		if (zone == 1001)
			std::cout<<"happy"<<std::endl;
		facetsRef.emplace( std::piecewise_construct, std::forward_as_tuple(i), std::forward_as_tuple(&mesh, tmp, (zone)*static_cast<UInt>(data.fractureOn()), bcId) );
		if(zone > 1 && data.fractureOn() && zones.find(zone) == zones.end())
		{
			prop.setProperties(1., 1., 1.);
			propMap.setZone(zone, prop);
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
	propMap.setZone(maxZone+1, prop);

	mesh.buildNodesToFacetMap();

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
			tmpFacets[j] = mesh.getFacetFromNodes(tmpNodes);
		}

		cellsRef.emplace( std::piecewise_construct, std::forward_as_tuple(i), std::forward_as_tuple(&mesh, tmpFacets, maxZone+1) );
	}
	tmp.clear();

	file.close();

	std::cout << "Add fractures" << std::endl;

	Geometry::FractureNetwork3D FN(mesh);
	std::vector<Geometry::Fracture3D> fracturesVector;
	std::map<UInt, Geometry::Fracture3D> fracturesMap;
	std::map<UInt, Geometry::Fracture3D>::iterator itF;
	Geometry::Point3D normal, center, centerFace;

	for(std::map<UInt, Geometry::Mesh3D::Facet3D>::iterator it = mesh.getFacetsMap().begin(); it != mesh.getFacetsMap().end(); ++it)
	{
		if (it->second.getZoneCode() > 0 && it->second.getBorderId()==0)
		{
			std::cout<<"zone: "<<it->second.getZoneCode()<<std::endl;
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

	std::cout << "Compute facet ids of the fractures..." << std::flush;
	mesh.updateFacetsWithFractures();
	std::cout << " done." << std::endl << std::endl;

	std::cout << "Export..." << std::flush;
	saveAsMeditFormat(data.getOutputDir() + data.getOutputFile() + ".mesh", mesh);
	std::cout << " done." << std::endl;

	return 0;
}
