/*!
 *  @file interCellIntersections.cpp
 *  @brief Class to store intersections between a cell and a fault (definition).
 *
 *  @author Luca Turconi <lturconi@gmail.com>
 *  @date 31-08-2012
 *
 */

#include<fstream>
#include<limits>
#include <iomanip>

#include "interCellIntersections.hpp"

namespace Intersect
{
  // --------------------   Class CellIntersections   --------------------

  // ==================================================
  // Constructors & Destructor
  // ==================================================
  CellIntersections::CellIntersections() : M_i(), M_j(), M_k(), M_intersections(),
    M_intersectionsError(), M_meanError(), M_maxError(),
    M_missedIntersections(), M_additionalIntersections() {}

  CellIntersections::CellIntersections (const Geometry::CPcell& c) :
    M_i (c.i() ), M_j (c.j() ), M_k (c.k() ), M_intersections(), M_intersectionsError(),
    M_meanError(), M_maxError(), M_missedIntersections(), M_additionalIntersections() {}

  CellIntersections::CellIntersections (const UInt& i, const UInt& j, const UInt& k) :
    M_i (i), M_j (j), M_k (k), M_intersections(), M_intersectionsError(),
    M_meanError(), M_maxError(), M_missedIntersections(), M_additionalIntersections() {}

  CellIntersections::~CellIntersections() {}

  // ==================================================
  // Methods
  // ==================================================
  bool CellIntersections::insert (const Intersection& inter, bool printError)
  {
    std::pair<CellIntersections_Iterator_Type, bool> ret;

    ret = M_intersections.insert ( inter );
    if (ret.second == false)
    {
      if (printError)
        std::cerr << " *** Error: insertion failed: element already existed! *** " << std::endl;
      return 0;
    }

    return 1;
  }

  bool CellIntersections::importIntersections (const CellIntersections& cellInter)
  {
    if ( cellInter.i() != M_i || cellInter.j() != M_j || cellInter.k() != M_k )
    {
      std::cerr << "Error: CellIntersections objects refer to different cells!" << std::endl;
      return 0;
    }

    for (CellIntersections_Const_Iterator_Type it = cellInter.begin();
         it != cellInter.end(); ++it)
    {
      this->insert (*it, 0);
    }
    return 1;
  }

  void CellIntersections::computeErrors (const CellIntersections& exactInter)
  {
    M_intersectionsError.clear();
    M_maxError.second = 0;
    M_maxError.first = 0;
    M_missedIntersections = 0;
    M_additionalIntersections = 0;

    Geometry::Point3D p (std::numeric_limits<Real>::quiet_NaN(),
                         std::numeric_limits<Real>::quiet_NaN(),
                         std::numeric_limits<Real>::quiet_NaN() );

    Real cumulError (0);
    UInt n (0);
    Intersection inter;

    CellIntersections_Const_Iterator_Type itEx;
    CellIntersections_Const_Iterator_Type it;

    for (UInt i = 1; i <= 12; ++i)
    {
      it = this->find ( i );
      itEx = exactInter.find ( i );

      if (it != this->end() && itEx != exactInter.end() )
      {
        inter.first = i;
        inter.second = itEx->second - it->second;
        M_intersectionsError.insert (inter);
        cumulError += inter.second.norm();
        ++n;
        if (inter.second.norm() > M_maxError.second)
        {
          M_maxError.second = inter.second.norm();
          M_maxError.first = inter.first;
        }
      }

      if (it == this->end() && itEx != exactInter.end() )
      {
        // only in exactInter --> -1
        p.x = -1;
        inter.first = i;
        inter.second = p;
        ++M_missedIntersections;
        M_intersectionsError.insert (inter);
      }

      if (it != this->end() && itEx == exactInter.end() )
      {
        // only in M_intersections --> +1
        p.x = +1;
        inter.first = i;
        inter.second = p;
        ++M_additionalIntersections;
        M_intersectionsError.insert (inter);
      }

    }

    if (n > 0)
      M_meanError = cumulError / n;
  }

  void CellIntersections::showErrors (std::ostream&   out) const
  {
    out << " Cell: (" << M_i << "," << M_j << "," << M_k << ")" << std::endl;
    out << " Errors: " << std::endl;
    for (std::map<UInt, Geometry::Point3D>::const_iterator it = M_intersectionsError.begin();
         it != M_intersectionsError.end(); ++it)
    {
      out << "    Edge : " << (*it).first
          << "  Error : " << (*it).second
          << "  ErrorNorm = " << (*it).second.norm() << std::endl;
    }
  }

  void CellIntersections::errorSummary (std::ostream&   out) const
  {
    out << " Cell (" << M_i << "," << M_j << "," << M_k << "): ERROR SUMMARY" << std::endl;
    out << " --------------------------------------- " << std::endl;
    out << " meanError = " << M_meanError << std::endl;
    out << " maxError = " << M_maxError.second << "  on Edge " << M_maxError.first << std::endl;
    out << " missedIntersections = " << M_missedIntersections << std::endl;
    out << " additionalIntersections = " << M_additionalIntersections << std::endl;
    out << " --------------------------------------- " << std::endl;
  }

  bool CellIntersections::exportVtk (const std::string& fileName) const
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

    UInt nPoints = this->size();
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

    for (CellIntersections_Const_Iterator_Type it = this->begin();
         it != this->end(); ++it)
    {
      filestr << it->second.x << " "
              << it->second.y << " "
              << it->second.z << std::endl;
    }


    filestr << std::endl;

    // Celldata
    filestr << "CELLS " << nPoints << " " << 2 * nPoints << std::endl;

    for (UInt i = 0; i < nPoints; ++i)
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

  void CellIntersections::showMe (std::ostream&   out) const
  {
    out << " Cell: (" << M_i << "," << M_j << "," << M_k << ")" << std::endl;
    out << " Intersections: " << std::endl;
    for (std::map<UInt, Geometry::Point3D>::const_iterator it = M_intersections.begin();
         it != M_intersections.end(); ++it)
    {
      out << "    Edge : " << (*it).first
          << "  Point : " << (*it).second << std::endl;
    }
  }

} // namespace Intersect

