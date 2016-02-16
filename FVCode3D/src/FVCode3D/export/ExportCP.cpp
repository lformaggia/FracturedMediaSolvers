/*
 * @file ExportCP.cpp
 * @brief Classes for saving files in the CP format (GRDECL) (definitions).
 */

#include <algorithm>
#include <functional>

#include <FVCode3D/export/ExportCP.hpp>

namespace FVCode3D
{

void ExporterCP::exportMesh(const Mesh3D & mesh, const std::string filename) const throw()
{

    std::ostringstream specGridSS, coordSS, zCornSS, actnumSS;

    const std::vector<Point3D> & pointsRef = mesh.getNodesVector();

    typedef std::pair<Real, Real> basePillar;
    typedef std::pair<Real, Real> baseIJPillar;
    typedef std::function<bool(basePillar, basePillar)> compPillarFct;

    compPillarFct compPillar = [](basePillar b1, basePillar b2)
    {
        if(b1.second == b2.second)
            return b1.first < b2.first;
        else
            return b1.second < b2.second;
    };

    std::map< basePillar, std::vector<Real>, compPillarFct > coordPillar( compPillar);
    std::map< baseIJPillar, basePillar, compPillarFct > ijPillar( compPillar );

    std::set<basePillar, compPillarFct > pillarSet(compPillar);

    for(auto it : pointsRef)
    {
        basePillar pairPillar{ it.x(), it.y() };
        pillarSet.insert( pairPillar );

        auto cIt = coordPillar.find( pairPillar );
        if ( cIt == std::end( coordPillar ) )
        {
            std::vector<Real> value ( 1, it.z() );
            coordPillar.insert( std::make_pair( pairPillar, value ) );
        } // if
        else
        {
            cIt->second.push_back( it.z() );
        } // else
    } // for

    for( auto it: coordPillar )
    {
        std::sort( std::begin( it.second ), std::end( it.second ) );
    } // for

    coordSS << "COORD" << std::endl;
    UInt Nx = 0, Ny = 0, Nz = 0;
    Real prevY = std::begin( coordPillar )->first.second;

    coordSS.precision(15);

    Int iPillar(-1), jPillar(0);
    for(auto it : coordPillar )
    {
        if(prevY == it.first.second)
        {
            if(Ny == 0) Nx++;
            iPillar++;
        }
        else
        {
            Ny++;
            jPillar++;
            iPillar = 0;
        }

        const Real x = it.first.first;
        const Real y = it.first.second;

        const Real zmin = *std::min_element( std::begin( it.second ),
                                             std::end( it.second ) );

        const Real zmax = *std::max_element( std::begin( it.second ),
                                             std::end( it.second ) );

        coordSS << x << " " << y << " " << zmin << " "
                << x << " " << y << " " << zmax
                << std::endl;

        baseIJPillar ij{ iPillar, jPillar };
        ijPillar.insert( std::make_pair( ij, it.first ) );

        prevY = it.first.second;
    } // for
    coordSS << "/" << std::endl << std::endl;

    Nx--;
    Nz = std::begin( coordPillar )->second.size() - 1;

    zCornSS.precision(15);
    zCornSS << "ZCORN" << std::endl;

    for ( UInt k = 0; k < Nz + 1; ++k )
    {
        UInt numJ(1);
        if( k != 0 && k != Nz ) numJ = 2;

        for( UInt jk = 0; jk < numJ; ++jk )
        {
            for( UInt j = 0; j < Ny + 1; ++j )
            {
                UInt numI(1);
                if( j != 0 && j != Ny ) numI = 2;

                for( UInt ij = 0; ij < numI; ++ij )
                {
                    for( UInt i = 0; i < Nx; ++i )
                    {
                        auto first = ijPillar.find( { i, j } )->second;
                        auto second = ijPillar.find( { i + 1, j } )->second;

                        const Real z1 = coordPillar.find( first )->second[ k ];
                        const Real z2 = coordPillar.find( second )->second[ k ];

                        zCornSS << z1 << " " << z2 << std::endl;
                    } // for
                } // for
            } // for
        } // for
    } // for

    zCornSS << std::endl << "/" << std::endl << std::endl;

    actnumSS << "ACTNUM" << std::endl
             << Nx*Ny*Nz << "*1 " << std::endl << "/" << std::endl;

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

} // exportMesh

} // namespace FVCode3D
