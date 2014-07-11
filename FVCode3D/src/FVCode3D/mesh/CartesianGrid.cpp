 /*!
 * @file cartesianGrid.cpp
 * @brief Class that generate a hexahedral structured (Cartesian) grid (definitions).
 */

#include <FVCode3D/mesh/Mesh3D.hpp>
#include <FVCode3D/property/Properties.hpp>
#include <FVCode3D/mesh/CartesianGrid.hpp>

namespace FVCode3D
{

void CartesianGrid::generate(bool fracturesOn, const Real Lx, const Real Ly, const Real Lz, const UInt Nx, const UInt Ny, const UInt Nz)
{
	Real hx = Lx/Nx;
	Real hy = Ly/Ny;
	Real hz = Lz/Nz;
	UInt nNodes = (Nx+1)*(Ny+1)*(Nz+1);
	UInt i,j,k;

	UInt bcId, zone, maxZone = 0, count = 0;

	std::vector<UInt> tmp, tmpFacets;

	std::vector<Point3D> & nodesRef = M_mesh.getNodesVector();
	std::map<UInt, Mesh3D::Facet3D> & facetsRef = M_mesh.getFacetsMap();
	std::map<UInt, Mesh3D::Cell3D> & cellsRef = M_mesh.getCellsMap();

	nodesRef.reserve(nNodes);

	// create nodes
	for(k=0; k <= Nz; ++k)
	{
		for(j=0; j <= Ny; ++j)
		{
			for(i=0; i <= Nx; ++i)
			{
				nodesRef.emplace_back(hx*i, hy*j, hz*k); // Point3D
			}
		}
	}

	tmp.resize(4);

	// create facets parallel to x axis
	for(i=0; i <= Nx; ++i)
	{
		for(k=0; k < Nz; ++k)
		{
			for(j=0; j < Ny; ++j)
			{
				tmp[0] = i + (Nx+1) * j + (Nx+1)*(Ny+1) * k;
				tmp[1] = i + (Nx+1) * j + (Nx+1)*(Ny+1) * (k+1);
				tmp[2] = i + (Nx+1) * (j+1) + (Nx+1)*(Ny+1) * (k+1);
				tmp[3] = i + (Nx+1) * (j+1) + (Nx+1)*(Ny+1) * k;

				bcId = ((i == 0) || (i == Nx)) ? 1 : 0;
				zone = 0;

				facetsRef.emplace( std::piecewise_construct, std::forward_as_tuple(count), std::forward_as_tuple(&M_mesh, tmp, (zone)*static_cast<UInt>(fracturesOn), bcId) );
				++count;
			}
		}
	}

	// create facets parallel to y axis
	for(j=0; j <= Ny; ++j)
	{
		for(i=0; i < Nx; ++i)
		{
			for(k=0; k < Nz; ++k)
			{
				tmp[0] = i + (Nx+1) * j + (Nx+1)*(Ny+1) * k;
				tmp[1] = i + (Nx+1) * j + (Nx+1)*(Ny+1) * (k+1);
				tmp[2] = i+1 + (Nx+1) * j + (Nx+1)*(Ny+1) * (k+1);
				tmp[3] = i+1 + (Nx+1) * j + (Nx+1)*(Ny+1) * k;

				bcId = (( j==0 ) || (j == Ny)) ? 1 : 0;
				zone = 0;

				facetsRef.emplace( std::piecewise_construct, std::forward_as_tuple(count), std::forward_as_tuple(&M_mesh, tmp, (zone)*static_cast<UInt>(fracturesOn), bcId) );
				++count;
			}
		}
	}

	// create facets parallel to z axis
	for(k=0; k <= Nz; ++k)
	{
		for(j=0; j < Ny; ++j)
		{
			for(i=0; i < Nx; ++i)
			{
				tmp[0] = i + (Nx+1) * j + (Nx+1)*(Ny+1) * k;
				tmp[1] = i+1 + (Nx+1) * j + (Nx+1)*(Ny+1) * k;
				tmp[2] = i+1 + (Nx+1) * (j+1) + (Nx+1)*(Ny+1) * k;
				tmp[3] = i + (Nx+1) * (j+1) + (Nx+1)*(Ny+1) * k;

				bcId = ((k == 0) || (k == Nz)) ? 1 : 0;
				zone = 0;

				facetsRef.emplace( std::piecewise_construct, std::forward_as_tuple(count), std::forward_as_tuple(&M_mesh, tmp, (zone)*static_cast<UInt>(fracturesOn), bcId) );
				++count;
			}
		}
	}

	tmpFacets.resize(6);

	count = 0;

	M_mesh.buildNodesToFacetMap();

	// create cells
	for(k=0; k < Nz; ++k)
	{
		for(j=0; j < Ny; ++j)
		{
			for(i=0; i < Nx; ++i)
			{
				// bottom
				tmp[0] = i + (Nx+1) * j + (Nx+1)*(Ny+1) * k;
				tmp[1] = i+1 + (Nx+1) * j + (Nx+1)*(Ny+1) * k;
				tmp[2] = i+1 + (Nx+1) * (j+1) + (Nx+1)*(Ny+1) * k;
				tmp[3] = i + (Nx+1) * (j+1) + (Nx+1)*(Ny+1) * k;

				tmpFacets[0] = M_mesh.getFacetFromNodes(tmp);

				// top
				tmp[0] = i + (Nx+1) * j + (Nx+1)*(Ny+1) * (k+1);
				tmp[1] = i+1 + (Nx+1) * j + (Nx+1)*(Ny+1) * (k+1);
				tmp[2] = i+1 + (Nx+1) * (j+1) + (Nx+1)*(Ny+1) * (k+1);
				tmp[3] = i + (Nx+1) * (j+1) + (Nx+1)*(Ny+1) * (k+1);

				tmpFacets[1] = M_mesh.getFacetFromNodes(tmp);

				// left
				tmp[0] = i + (Nx+1) * j + (Nx+1)*(Ny+1) * k;
				tmp[1] = i + (Nx+1) * j + (Nx+1)*(Ny+1) * (k+1);
				tmp[2] = i+1 + (Nx+1) * j + (Nx+1)*(Ny+1) * (k+1);
				tmp[3] = i+1 + (Nx+1) * j + (Nx+1)*(Ny+1) * k;

				tmpFacets[2] = M_mesh.getFacetFromNodes(tmp);

				// right
				tmp[0] = i + (Nx+1) * (j+1) + (Nx+1)*(Ny+1) * k;
				tmp[1] = i + (Nx+1) * (j+1) + (Nx+1)*(Ny+1) * (k+1);
				tmp[2] = i+1 + (Nx+1) * (j+1) + (Nx+1)*(Ny+1) * (k+1);
				tmp[3] = i+1 + (Nx+1) * (j+1) + (Nx+1)*(Ny+1) * k;

				tmpFacets[3] = M_mesh.getFacetFromNodes(tmp);

				// back
				tmp[0] = i + (Nx+1) * j + (Nx+1)*(Ny+1) * k;
				tmp[1] = i + (Nx+1) * j + (Nx+1)*(Ny+1) * (k+1);
				tmp[2] = i + (Nx+1) * (j+1) + (Nx+1)*(Ny+1) * (k+1);
				tmp[3] = i + (Nx+1) * (j+1) + (Nx+1)*(Ny+1) * k;

				tmpFacets[4] = M_mesh.getFacetFromNodes(tmp);

				// front
				tmp[0] = i+1 + (Nx+1) * j + (Nx+1)*(Ny+1) * k;
				tmp[1] = i+1 + (Nx+1) * j + (Nx+1)*(Ny+1) * (k+1);
				tmp[2] = i+1 + (Nx+1) * (j+1) + (Nx+1)*(Ny+1) * (k+1);
				tmp[3] = i+1 + (Nx+1) * (j+1) + (Nx+1)*(Ny+1) * k;

				tmpFacets[5] = M_mesh.getFacetFromNodes(tmp);

				cellsRef.emplace( std::piecewise_construct, std::forward_as_tuple(count), std::forward_as_tuple(&M_mesh, tmpFacets, maxZone+1) );

				++count;
			}
		}
	}
}

void CartesianGrid::extractBC(const Real theta)
{
	Point3D normal, center, centerFace;
	Real max;
	UInt compMax;

	for(std::map<UInt, Mesh3D::Facet3D>::iterator it = M_mesh.getFacetsMap().begin(); it != M_mesh.getFacetsMap().end(); ++it)
	{
		if (it->second.getSeparatedCells().size() == 1)
		{
			center = M_mesh.getCellsMap().at(*(it->second.getSeparatedCells().begin())).getCentroid();
			centerFace = it->second.getCentroid();

			normal = it->second.getUnsignedNormal();
			center -= centerFace;
			center.normalize();

			if(normal * center > 0.)
				normal = -normal;

			normal = rotateOf(normal, theta);

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

void CartesianGrid::addFractures(const std::map<UInt,UInt> & facetIdToZone)
{
	std::map<UInt, Mesh3D::Facet3D> & facetsRef = M_mesh.getFacetsMap();

	FractureNetwork3D FN(M_mesh);
	std::vector<Fracture3D> fracturesVector;
	std::map<UInt, Fracture3D> fracturesMap;
	std::map<UInt, Fracture3D>::iterator itF;
	Point3D normal, center, centerFace;

	for(std::map<UInt,UInt>::const_iterator it = facetIdToZone.begin(); it!= facetIdToZone.end(); ++it)
		if(facetsRef.find(it->first) != facetsRef.end())
			facetsRef[it->first].setZoneCode(it->second);

	for(std::map<UInt, Mesh3D::Facet3D>::iterator it = facetsRef.begin(); it != facetsRef.end(); ++it)
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

void CartesianGrid::addBCAndFractures(const std::map<UInt,UInt> & facetIdToZone, const Real theta)
{
	extractBC(theta);
	addFractures(facetIdToZone);
}

} // namespace FVCode3D
