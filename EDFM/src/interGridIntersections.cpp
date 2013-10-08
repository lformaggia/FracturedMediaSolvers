/*!
*  @file interGridIntersections.cpp
*  @brief Class to store intersections between a grid and a fault (definition).
*
*  @author Luca Turconi <lturconi@gmail.com>
*  @date 31-08-2012
*
*/

#include<fstream>
#include<limits>
#include <iomanip>

#include "interGridIntersections.hpp"


namespace Intersect
{
  // --------------------   Class GridIntersections   --------------------

  // ==================================================
  // Constructors & Destructor
  // ==================================================
  GridIntersections::GridIntersections (const Geometry::CPgrid& g) :
    M_CellIntersectionsMap(), M_Nx (g.Nx() ), M_Ny (g.Ny() ), M_Nz (g.Nz() ), M_Nintersections (0),
    M_meanError(), M_maxError(), M_missedCells(), M_additionalCells() {}

  GridIntersections::GridIntersections (const UInt& Nx, const UInt& Ny, const UInt& Nz) :
    M_CellIntersectionsMap(), M_Nx (Nx), M_Ny (Ny), M_Nz (Nz), M_Nintersections (0),
    M_meanError(), M_maxError(), M_missedCells(), M_additionalCells() {}

  GridIntersections::~GridIntersections() {}

  // ==================================================
  // Get Methods
  // ==================================================
  CellIntersections GridIntersections::getCellIntersections (const UInt& id) const
  {
    if ( M_CellIntersectionsMap.find (id) == M_CellIntersectionsMap.end() )
    {
      std::cerr << " *** Error: cell with id = " << id
                << " has no stored intersections ***" << std::endl;
      CellIntersections err;
      err.i() = 0;
      err.j() = 0;
      err.k() = 0;
      return err;
    }

    return M_CellIntersectionsMap.find (id)->second;
  }

  CellIntersections GridIntersections::getCellIntersections (const UInt& i, const UInt& j, const UInt& k) const
  {
    if ( i > M_Nx )
    {
      std::cerr << " *** Error: index i out of range: i = " << i << " and Nx = "
                << M_Nx << " ***" << std::endl;
      CellIntersections err;
      err.i() = 0;
      err.j() = 0;
      err.k() = 0;
      return err;
    }
    if ( j > M_Ny )
    {
      std::cerr << " *** Error: index j out of range: j = " << j << " and Ny = "
                << M_Ny << " ***" << std::endl;
      CellIntersections err;
      err.i() = 0;
      err.j() = 0;
      err.k() = 0;
      return err;
    }
    if ( i > M_Nz )
    {
      std::cerr << " *** Error: index k out of range: k = " << k << " and Nz = "
                << M_Nz << " ***" << std::endl;
      CellIntersections err;
      err.i() = 0;
      err.j() = 0;
      err.k() = 0;
      return err;
    }

    UInt id = i + (j - 1) * M_Nx + (k - 1) * M_Nx * M_Ny;

    return this->getCellIntersections (id);
  }

  // ==================================================
  // Methods
  // ==================================================
  bool GridIntersections::insert (const CellIntersections& inter, bool printError)
  {
    UInt id = inter.i() + (inter.j() - 1) * M_Nx + (inter.k() - 1) * M_Nx * M_Ny;
    std::pair<GridIntersections_Iterator_Type, bool> ret;

    ret = M_CellIntersectionsMap.insert ( std::pair<UInt, CellIntersections> (id, inter) );
    if (ret.second == false)
    {
      if (printError)
        std::cerr << " *** Error: insertion failed: element already existed! *** " << std::endl;
      return 0;
    }

    if (ret.second == true)
      M_Nintersections += inter.size();

    return 1;
  }

  void GridIntersections::clearAll( )
  {
    M_Nintersections = 0;
    M_meanError = 0;
    M_missedCells = 0;
    M_additionalCells = 0;
    M_maxError.first = 0;
    M_maxError.second = 0;
    M_CellIntersectionsMap.clear();
  }

  void GridIntersections::computeErrors (const GridIntersections& exactInter)
  {
    M_maxError.first = 0;
    M_maxError.second = 0;
    M_missedCells = 0;
    M_additionalCells = 0;

    Real cumulError (0);
    UInt n (0);

    GridIntersections_Iterator_Type it;
    GridIntersections_Const_Iterator_Type itEx;

    for (UInt id = 1; id <= M_Nx * M_Ny * M_Nz; ++id)
    {
      it = this->find ( id );
      itEx = exactInter.find ( id );

      if (it != this->end() && itEx != exactInter.end() )
      {
        it->second.computeErrors (itEx->second);

        cumulError += it->second.meanError();
        ++n;

        if (it->second.maxError().second > M_maxError.second)
        {
          M_maxError.second = it->second.maxError().second;
          M_maxError.first = id;
        }
      }

      if (it == this->end() && itEx != exactInter.end() )
      {
        ++M_missedCells;
      }

      if (it != this->end() && itEx == exactInter.end() )
      {
        ++M_additionalCells;
      }
    }

    M_meanError = cumulError / n;
  }

  void GridIntersections::errorSummary (std::ostream&   out) const
  {
    out << " GridIntersections: ERROR SUMMARY" << std::endl;
    out << " --------------------------------------- " << std::endl;
    out << " gridDimensions: " << M_Nx << " x " << M_Ny << " x " << M_Nz << std::endl;
    out << " Intersections Found: " << M_Nintersections << std::endl;
    out << " meanError = " << M_meanError << std::endl;
    out << " maxError = " << M_maxError.second << "  in Cell (id) "
        << M_maxError.first << " --> (" << i (M_maxError.first) << ","
        << j (M_maxError.first) << "," << k (M_maxError.first) << ")" << std::endl;
    out << " missedCells = " << M_missedCells << std::endl;
    out << " additionalCells = " << M_additionalCells << std::endl;
    out << " --------------------------------------- " << std::endl;
  }

  bool GridIntersections::exportVtk (const std::string& fileName) const
  {
    std::fstream filestr;

    filestr.open (fileName.c_str(), std::ios_base::out);

    if (filestr.is_open() )
    {
      std::cout << std::endl << " File: " << fileName << ", successfully opened";
    }
    else
    {
      std::cerr << std::endl << " *** Error: file not opened *** " << std::endl << std::endl;
      return  0;
    }

    std::cout << std::endl << " Exporting grid in Vtk format... " << std::endl;

    UInt nPoints = M_Nintersections;
    UInt CellType = 1; // for VTK_POINT

    // Header
    filestr << "# vtk DataFile Version 3.1" << std::endl;
    filestr << "this is a file created for Paraview" << std::endl;
    filestr << "ASCII" << std::endl;
    filestr << "DATASET UNSTRUCTURED_GRID" << std::endl;
    filestr << std::endl; // The fifth line is empty.

    filestr << std::scientific << std::setprecision (10);

    // Pointdata
    filestr << "POINTS " << nPoints << " double" << std::endl;

    for (GridIntersections_Const_Iterator_Type it = this->begin();
         it != this->end(); ++it)
    {
      for (CellIntersections_Const_Iterator_Type jt = (*it).second.begin();
           jt != (*it).second.end(); ++jt)
      {
        filestr << jt->second.x << " "
                << jt->second.y << " "
                << jt->second.z << std::endl;
      }
    }

    filestr << std::endl;

    // Celldata
    filestr << "CELLS " << nPoints << " " << 2 * nPoints << std::endl;

    for (UInt i = 0; i < M_Nintersections; ++i)
    {
      filestr << "1 " << i << std::endl;
    }

    filestr << std::endl;

    filestr << "CELL_TYPES " << nPoints << std::endl;

    for (UInt i = 0; i < nPoints; ++i)
      filestr << CellType << std::endl;
    filestr << std::endl;

    filestr.close();

    return 1;
  }

  void GridIntersections::showMe (std::ostream&   out) const
  {
    out << " GRID INTERSECTIONS: " << std::endl;
    out << "  " << M_Nintersections << " Points" << std::endl;
    GridIntersections_Const_Iterator_Type it;

    for (it = this->begin(); it != this->end(); ++it)
    {
      (*it).second.showMe();
      out << " -- o -- o -- o -- o -- o -- o -- o -- " << std::endl;
    }
  }

} // namespace Intersect