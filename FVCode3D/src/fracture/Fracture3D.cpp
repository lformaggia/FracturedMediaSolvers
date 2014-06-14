/*!
 * @file Fracture3D.cpp
 * @brief Class that handles the fractures (definition).
 */

#include "fracture/Fracture3D.hpp"
#include "mesh/Mesh3D.hpp"

namespace FVCode3D
{

Fracture3D::Fracture3D(const Geometry::Mesh3D & mesh):
	M_id(0), M_mesh(mesh) {}

Fracture3D::Fracture3D(const Geometry::Fracture3D & f):
	M_id(f.getId()), M_fractureFacets(f.getFractureFacetsId()), M_mesh(f.getMesh()) {}

Fracture3D::Fracture3D(const Geometry::Mesh3D & mesh, const std::vector<UInt> & fractureFacets, const UInt id):
	M_id(id), M_fractureFacets(fractureFacets), M_mesh(mesh) {}

bool Fracture3D::exportVTK(const std::string & filename) const
{
	// TODO Fare una bool per scegliere se fare out o append?

	const UInt nFaces = M_fractureFacets.size();
	UInt totalNodes = 0;
	std::set<UInt> nodes;

	for(UInt i=0; i < nFaces; ++i)
	{
		totalNodes += (M_mesh.getFacetsMap()).at(M_fractureFacets[i]).getNumberOfVertices();
		nodes.insert( M_mesh.getFacetsMap().at(M_fractureFacets[i]).getVerticesVector().begin(), M_mesh.getFacetsMap().at(M_fractureFacets[i]).getVerticesVector().end() );
		//std::copy( M_mesh.getFacetsMap().at(M_fractureFacets[i]).getVerticesVector().begin(), M_mesh.getFacetsMap().at(M_fractureFacets[i]).getVerticesVector().end(),
		//		std::inserter( nodes, nodes.end() ) );
	}

	std::ofstream filestr;

	filestr.open (filename.c_str(), std::ios_base::out | std::ios_base::app);

	if (filestr.is_open())
	{
		std::cout << std::endl << " File: " << filename << ", successfully opened";
	}
	else
	{
		std::cerr << std::endl
				  << " *** Error: file not opened *** " << std::endl << std::endl;
		return  0;
	}

	std::cout << std::endl
			  << " Exporting Fracture3D discretization in Vtk format... " << std::endl;

	const UInt nPoints = nodes.size();
	const UInt CellType = 7; // for VTK_POLYGON
	UInt nNodesFacet;

	// Header
	filestr << "# vtk DataFile Version 3.1" << std::endl;
	filestr << "this is a file created for Paraview" << std::endl;
	filestr << "ASCII" << std::endl;
	filestr << "DATASET UNSTRUCTURED_GRID" << std::endl;
	filestr << std::endl;	// The fifth line is empty.

	filestr << std::scientific << std::setprecision(10);

	// Pointdata
	filestr << "POINTS " << nPoints << " double" << std::endl;
	for(UInt i=0; i < nPoints; ++i)
		filestr << M_mesh.getNodesVector()[i].x() << " " << M_mesh.getNodesVector()[i].y() << " " << M_mesh.getNodesVector()[i].z() << std::endl;
	filestr << std::endl;

	// Celldata
	filestr << "CELLS " << nFaces << " " << totalNodes + nFaces << std::endl;
	for(UInt i=0; i<nFaces; ++i){
		nNodesFacet = M_mesh.getFacetsMap().at(M_fractureFacets[i]).getNumberOfVertices();
		filestr << nNodesFacet;
		for(UInt j=0; j < nNodesFacet; ++j)
			filestr << M_mesh.getFacetsMap().at(M_fractureFacets[i]).getVerticesVector()[j];
		filestr << std::endl;
	}
	filestr << std::endl;

	// Celltype
	filestr << "CELL_TYPES " << nFaces << std::endl;
	for(UInt i=0; i < nFaces; ++i)
		filestr << CellType << std::endl;
	filestr << std::endl;

	filestr.close();

	return true;
}

void Fracture3D::showMe(std::ostream & out) const
{
	out << "Type = Fracture3D " << std::endl;
	out << " M_fractureFacets :" << std::endl;
	out << " [ ";
	for(UInt i=0; i < M_fractureFacets.size()-1 ; ++i)
		out << M_fractureFacets[i] << " , ";
	out << M_fractureFacets[M_fractureFacets.size() - 1] << " ] " << std::endl;
	out << " Physical Properties : " << std::endl;
}

} // namespace FVCode3D
