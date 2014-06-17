/*!
 * @file FractureNetwork3D.cpp
 * @brief Class that handles the fracture network (definitions).
 */

#include "fracture/FractureNetwork3D.hpp"
#include "mesh/Mesh3D.hpp"
#include "fracture/Fracture3D.hpp"

#include <fstream>
#include <limits>
#include <cstring>

namespace FVCode3D
{

FractureNetwork3D::FractureNetwork3D(const Mesh3D & mesh):
		M_mesh(mesh) {}

FractureNetwork3D::FractureNetwork3D(const FractureNetwork3D & fn):
		M_fractureNetwork(fn.getNetwork()), M_mesh(fn.getMesh()) {}

FractureNetwork3D::FractureNetwork3D(const Mesh3D & mesh, const std::vector<Fracture3D> & fractures):
		M_fractureNetwork(fractures), M_mesh(mesh) {}

void FractureNetwork3D::addFractures(std::vector<Fracture3D> & fractures)
{
	M_fractureNetwork.clear();
	std::move(fractures.begin(), fractures.end(), std::back_inserter(M_fractureNetwork));
	fractures.erase(fractures.begin(), fractures.end());
}

bool FractureNetwork3D::exportVTK(const std::string & prefixFileName) const
{
	const UInt nFN = M_fractureNetwork.size();
	bool status = true;

	for( UInt i=0; i < nFN && status; ++i)
    {
        std::stringstream vtkFileName;
        vtkFileName << prefixFileName << i << ".vtk";
		status = M_fractureNetwork[i].exportVTK(vtkFileName.str());
    }

	return status;
}

bool FractureNetwork3D::exportNetworkVTK(const std::string & filename) const
{
	const UInt nFN = M_fractureNetwork.size();
	bool status = true;

	for(UInt i=0; i < nFN && status; ++i)
		status = M_fractureNetwork[i].exportVTK(filename);

	return status;
}

void FractureNetwork3D::showMe(std::ostream & out) const
{
	const UInt nFN = M_fractureNetwork.size();
	out << "Type = FractureNetwork3D " << std::endl;
	out << std::endl;
	for( UInt i = 0; i < nFN; ++i )
	{
		out << " f[" << i << "] : ";
		M_fractureNetwork[i].showMe(out);
		out << std::endl;
	}

	out << std::endl;
}

} // namespace FVCode3D
