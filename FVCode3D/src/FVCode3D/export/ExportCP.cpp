/*
 * @file ExportCP.cpp
 * @brief Classes for saving files in the CP format (GRDECL) (definitions).
 */

#include <algorithm>
#include <functional>

#include <FVCode3D/export/ExportCP.hpp>

namespace FVCode3D
{

void ExporterCP::exportMesh(const Mesh3D & mesh, const std::string filename) throw()
{
    std::ostringstream specGridSS, coordSS, zCornSS, actnumSS;

    const std::vector<Point3D> & pointsRef = mesh.getNodesVector();

    auto conf = [](Point3D p1, Point3D p2) { return p1.z() < p2.z(); };

    Real zmax, zmin;

    zmax = (std::max_element(std::begin(pointsRef), std::end(pointsRef), conf))->z();
    zmin = (std::min_element(std::begin(pointsRef), std::end(pointsRef), conf))->z();

    typedef std::pair<Real, Real> basePillar;
    std::function<bool(basePillar, basePillar)> compPillar = [](basePillar b1, basePillar b2)
    {
        if(b1.second == b2.second)
            return b1.first < b2.first;
        else
            return b1.second < b2.second;
    };

    std::set<basePillar, std::function<bool(basePillar, basePillar)> > pillarSet(compPillar);

    for(auto it : pointsRef)
        pillarSet.insert( std::make_pair(it.x(), it.y()) );

    coordSS << "COORD" << std::endl;
    UInt Nx = 0, Ny = 0, Nz = 0;
    Real prevY = std::begin(pillarSet)->second;

    coordSS.precision(15);

    for(auto it : pillarSet)
    {
        if(prevY == it.second)
        {
            if(Ny == 0)
                Nx++;
        }
        else
            Ny++;

        coordSS << it.first << " " << it.second << " " << zmin << " "
                << it.first << " " << it.second << " " << zmax << std::endl;

        prevY = it.second;
    }
    coordSS << "/" << std::endl << std::endl;

    Nx--;

    basePillar prevPillar = *std::begin(pillarSet);

    Nz = std::count_if(std::begin(pointsRef), std::end(pointsRef),
            [&prevPillar](Point3D p){ return (p.x() == prevPillar.first) && (p.y() == prevPillar.second); }) - 1;

    zCornSS << "ZCORN" << std::endl;

    zCornSS << 4*(Nx+1)*(Ny+1) - 4*(Nx+1) - 4*(Ny-1) - 4 << "*" << zmin << " ";

    for(UInt i=1; i < Nz; ++i)
        zCornSS << 2*(4*(Nx+1)*(Ny+1) - 4*(Nx+1) - 4*(Ny-1) - 4) << "*" << (zmax-zmin)/Nz*i+zmin << " ";

    zCornSS << 4*(Nx+1)*(Ny+1) - 4*(Nx+1) - 4*(Ny-1) - 4 << "*" << zmax << std::endl << "/" << std::endl << std::endl;

    actnumSS << "ACTNUM" << std::endl;

    actnumSS << Nx*Ny*Nz << "*1 " << "/" << std::endl;

    std::fstream cpFile;

    cpFile.open (filename.c_str(), std::ios_base::out);

    cpFile.precision(15);

    if (cpFile.is_open())
    {
        std::cout << std::endl << " File: " << filename << ", successfully opened";
    }
    else
    {
    	throw std::runtime_error("Error: file " + filename + " not opened.");
    }

    std::cout << std::endl << " Exporting Mesh in GRDECL format... " << std::endl;

    specGridSS << "SPECGRID" << std::endl;
    specGridSS << Nx << " " << Ny << " " << Nz << " /" << std::endl << std::endl;

    cpFile << specGridSS.str() << coordSS.str() << zCornSS.str() << actnumSS.str();

    cpFile.close();
}

} // namespace FVCode3D
