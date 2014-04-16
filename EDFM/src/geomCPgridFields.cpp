/*!
*  @file geomCPgridFields.cpp
*  @brief Classes for permeability and transmissibility field on Corner Point Grid (definition).
*
*  @author Luca Turconi <lturconi@gmail.com>
*  @date 22-09-2012
*
*/

#include "EclipseFileTranslator.hpp"
#include "geomCPgridFields.hpp"
#include "geomCPgrid.hpp"


namespace Geometry
{

  // --------------------   Class PermeabilityField   --------------------

  // ==================================================
  // Constructors & Destructor
  // ==================================================
  PermeabilityField::PermeabilityField() :
    M_Nx (0), M_Ny (0), M_Nz (0), M_isotropic (1), M_permx(), M_permy(), M_permz() {}

  PermeabilityField::PermeabilityField (const UInt& nx, const UInt& ny, const UInt& nz, const bool& iso) :
    M_Nx (nx), M_Ny (ny), M_Nz (nz), M_isotropic (iso), M_permx(), M_permy(), M_permz() {}

  PermeabilityField::PermeabilityField (const CPgrid& grid, const bool& iso) :
    M_Nx (grid.Nx() ), M_Ny (grid.Ny() ), M_Nz (grid.Nz() ), M_isotropic (iso) {}

  PermeabilityField::~PermeabilityField() {}

  // ==================================================
  // Get Methods
  // ==================================================
  Vector3D PermeabilityField::getCellPermVec (const UInt& i, const UInt& j, const UInt& k) const
  {
    Vector3D v (getCellPermX (i, j, k), getCellPermY (i, j, k), getCellPermZ (i, j, k) );
    return v;
  }

  // ==================================================
  // Set Methods
  // ==================================================
  void PermeabilityField::setPermx (const std::vector<Real>& p)
  {
    UInt len = M_Nx * M_Ny * M_Nz ;
    if (p.size() != len)
      std::cerr << " *** Error: in PermX dimension *** "
                << std::endl << std::endl;
    else
      M_permx = p;
  }

  void PermeabilityField::setPermy (const std::vector<Real>& p)
  {
    UInt len = M_Nx * M_Ny * M_Nz ;
    if (p.size() != len)
    {
      std::cerr << " *** Error: in PermY dimension *** "
                << std::endl << std::endl;
    }
    else
    {
      M_isotropic = 0;
      M_permy = p;
    }
  }

  void PermeabilityField::setPermz (const std::vector<Real>& p)
  {
    UInt len = M_Nx * M_Ny * M_Nz ;
    if (p.size() != len)
    {
      std::cerr << " *** Error: in PermZ dimension *** "
                << std::endl << std::endl;
    }
    else
    {
      M_isotropic = 0;
      M_permz = p;
    }
  }

  // ==================================================
  // Methods
  // ==================================================
  bool PermeabilityField::importFromFile (const std::string& filename, const bool& iso,
                                          const bool& scriptWithSpecialChar)
  {
    UInt nval = M_Nx * M_Ny * M_Nz;
    M_isotropic = iso;

    std::cout << std::endl << "Reading section PERMX..." << std::endl;
    if ( !EclipseFile::readSection (filename, "PERMX", M_permx, scriptWithSpecialChar, nval) )
    {
      std::cerr << " *** Error: in reading section PERMX! *** " << std::endl;
      return 0;
    }

    if (!M_isotropic)
    {
      std::cout << std::endl << "Reading section PERMY..." << std::endl;
      if ( !EclipseFile::readSection (filename, "PERMY", M_permy, scriptWithSpecialChar, nval) )
      {
        std::cerr << " *** Error: in reading section PERMY! *** " << std::endl;
        return 0;
      }

      std::cout << std::endl << "Reading section PERMZ..." << std::endl;
      if ( !EclipseFile::readSection (filename, "PERMZ", M_permz, scriptWithSpecialChar, nval) )
      {
        std::cerr << " *** Error: in reading section PERMZ! *** " << std::endl;
        return 0;
      }
    }

    return 1;
  }

  bool PermeabilityField::exportToFile (const std::string& filename,
                                        const std::ios_base::openmode& mode)
  {
    if (M_permx.size() == 0)
    {
      std::cerr << std::endl << " *** Error: Permeability Field not defined *** "
                << std::endl << std::endl;
      return  0;
    }

    std::cout << std::endl << "Writing section PERMX..." << std::endl;
    if ( !EclipseFile::writeSection (M_permx, "PERMX", filename, mode) )
    {
      std::cerr << std::endl << " *** Error: in writing section PERMX *** "
                << std::endl << std::endl;
      return  0;
    }

    if (!M_isotropic)
    {
      std::cout << std::endl << "Writing section PERMY..." << std::endl;
      if ( !EclipseFile::writeSection (M_permy, "PERMY", filename) )
      {
        std::cerr << std::endl << " *** Error: in writing section PERMY *** "
                  << std::endl << std::endl;
        return  0;
      }

      std::cout << std::endl << "Writing section PERMZ..." << std::endl;
      if ( !EclipseFile::writeSection (M_permz, "PERMZ", filename) )
      {
        std::cerr << std::endl << " *** Error: in writing section PERMZ *** "
                  << std::endl << std::endl;
        return  0;
      }
    }

    return 1;
  }

  void PermeabilityField::showMe (std::ostream&   out) const
  {
    out << "Type = PermeabilityField : " << std::endl;
    out << " Dim = (Nx,Ny,Nz) : Type = (UInt,UInt,UInt) : ( "
        << M_Nx << " , " << M_Ny << " , " << M_Nz << " )" << std::endl;
    out << " Isotropic Field : ";
    if (M_isotropic)
    {
      out << "TRUE" << std::endl;
    }
    else
    {
      out << "FALSE" << std::endl;
    }
    out << " PermX : Type = vector<Real> : size = " << M_permx.size() << std::endl;
    if (!M_isotropic)
    {
      out << " PermY : Type = vector<Real> : size = " << M_permy.size() << std::endl;
      out << " PermZ : Type = vector<Real> : size = " << M_permz.size() << std::endl;
    }
  }


  // --------------------   Class TransmissibilityField   --------------------

  // ==================================================
  // Constructors & Destructor
  // ==================================================
  TransmissibilityField::TransmissibilityField() :
    M_Nx (0), M_Ny (0), M_Nz (0), M_tranx(), M_trany(), M_tranz() {}

  TransmissibilityField::TransmissibilityField (const UInt& nx, const UInt& ny, const UInt& nz) :
    M_Nx (nx), M_Ny (ny), M_Nz (nz), M_tranx(), M_trany(), M_tranz() {}

  TransmissibilityField::TransmissibilityField (const CPgrid& grid) :
    M_Nx (grid.Nx() ), M_Ny (grid.Ny() ), M_Nz (grid.Nz() ) {}

  TransmissibilityField::~TransmissibilityField() {}

  // ==================================================
  // Get Methods
  // ==================================================
  Vector3D TransmissibilityField::getCellTranVec (const UInt& i, const UInt& j, const UInt& k) const
  {
    Vector3D v (getCellTranX (i, j, k), getCellTranY (i, j, k), getCellTranZ (i, j, k) );
    return v;
  }

  // ==================================================
  // Set Methods
  // ==================================================
  void TransmissibilityField::setTranx (const std::vector<Real>& p)
  {
    UInt len = M_Nx * M_Ny * M_Nz - 1;
    if (p.size() != len)
      std::cerr << " *** Error: in TranX dimension *** "
                << std::endl << std::endl;
    else
      M_tranx = p;
  }

  void TransmissibilityField::setTrany (const std::vector<Real>& p)
  {
    UInt len = M_Nx * M_Ny * M_Nz - 1;
    if (p.size() != len)
    {
      std::cerr << " *** Error: in TranY dimension *** "
                << std::endl << std::endl;
    }
    else
    {
      M_trany = p;
    }
  }

  void TransmissibilityField::setTranz (const std::vector<Real>& p)
  {
    UInt len = M_Nx * M_Ny * M_Nz - 1;
    if (p.size() != len)
    {
      std::cerr << " *** Error: in TranZ dimension *** "
                << std::endl << std::endl;
    }
    else
    {
      M_tranz = p;
    }
  }

  // ==================================================
  // Methods
  // ==================================================
  bool TransmissibilityField::importFromFile (const std::string& filename,
                                              const bool& scriptWithSpecialChar)
  {
    UInt nval = M_Nx * M_Ny * M_Nz;

    std::cout << std::endl << "Reading section TRANX..." << std::endl;
    if ( !EclipseFile::readSection (filename, "TRANX", M_tranx, scriptWithSpecialChar, nval) )
    {
      std::cerr << " *** Error: in reading section TRANX! *** " << std::endl;
      return 0;
    }

    std::cout << std::endl << "Reading section TRANY..." << std::endl;
    if ( !EclipseFile::readSection (filename, "TRANY", M_trany, scriptWithSpecialChar, nval) )
    {
      std::cerr << " *** Error: in reading section TRANY! *** " << std::endl;
      return 0;
    }

    std::cout << std::endl << "Reading section TRANZ..." << std::endl;
    if ( !EclipseFile::readSection (filename, "TRANZ", M_tranz, scriptWithSpecialChar, nval) )
    {
      std::cerr << " *** Error: in reading section TRANZ! *** " << std::endl;
      return 0;
    }

    return 1;
  }

  bool TransmissibilityField::exportToFile (const std::string& filename,
                                            const std::ios_base::openmode& mode)
  {
    if (M_tranx.size() == 0)
    {
      std::cerr << std::endl << " *** Error: Permeability Field not defined *** "
                << std::endl << std::endl;
      return  0;
    }

    std::cout << std::endl << "Writing section TRANX..." << std::endl;
    if ( !EclipseFile::writeSection (M_tranx, "TRANX", filename, mode) )
    {
      std::cerr << std::endl << " *** Error: in writing section TRANX *** "
                << std::endl << std::endl;
      return  0;
    }


    std::cout << std::endl << "Writing section TRANY..." << std::endl;
    if ( !EclipseFile::writeSection (M_trany, "TRANY", filename) )
    {
      std::cerr << std::endl << " *** Error: in writing section TRANY *** "
                << std::endl << std::endl;
      return  0;
    }

    std::cout << std::endl << "Writing section TRANZ..." << std::endl;
    if ( !EclipseFile::writeSection (M_tranz, "TRANZ", filename) )
    {
      std::cerr << std::endl << " *** Error: in writing section TRANZ *** "
                << std::endl << std::endl;
      return  0;
    }

    return 1;
  }

  void TransmissibilityField::showMe (std::ostream&   out) const
  {
    out << "Type = TransmissibilityField : " << std::endl;
    out << " Dim = (Nx,Ny,Nz) : Type = (UInt,UInt,UInt) : ( "
        << M_Nx << " , " << M_Ny << " , " << M_Nz << " )" << std::endl;

    out << " TranX : Type = vector<Real> : size = " << M_tranx.size() << std::endl;
    out << " TranY : Type = vector<Real> : size = " << M_trany.size() << std::endl;
    out << " TranZ : Type = vector<Real> : size = " << M_tranz.size() << std::endl;

  }

} // namespace Geometry
