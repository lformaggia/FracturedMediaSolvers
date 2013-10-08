/*!
*  @file geomCPgrid.hpp
*  @brief Class for Corner Point Grid management.
*
*  @author Luca Turconi <lturconi@gmail.com>
*  @date 22-09-2012
*
*/

#ifndef GEOMCPGRID_HPP_
#define GEOMCPGRID_HPP_

#include<vector>
#include<iostream>
#include<string>
#include<fstream>
#include<memory>
#include "TypeDefinition.hpp"
#include "geomCPgridElements.hpp"
#include "geomCPgridFields.hpp"
#include "adtree.hpp"
#include <iterator>

namespace Geometry
{
  /*!
    @class CPgrid

    @author Luca Turconi <lturconi@gmail.com>, Anna Scotti, Luca Formaggia

      This class implements the concept of Corner Point grid.

      Geometric informations about the grid are saved in the two vector:
      M_coord and M_zcorn.
    M_coord contains the coordinates of the points generating the pillars.
    M_zcorn stores the coordinates z of the cell vertices.

    It is possible to define on the grid a permeability and
    a transmissibility field.

      Tools for import-export of the grid and its fileds
      in Eclipse data file format are available.

      The class provides methods to extract cells and pillars.

      The method showMe allows to print the class informations.

      For 3D visualization use the exportVtk method.

    */
  class CPgrid
  {
  public:
    //! @name Constructor & Destructor
    //@{

    //! Constructor, reading from file
    /*!
     * @param filename The Eclipse data file with grid informations
     * @param scriptWithSpecialChar Set TRUE if the Eclipse data file is generated on windows
     */
    CPgrid (const std::string& filename, const bool& scriptWithSpecialChar = 1, std::string direzione = "x", Real angle = 0, bool rotate_z = false);

    //! Destructor
    ~CPgrid();

    //@}

    //! @name Get Methods
    //@{

    //! Get grid dimension in x direction
    /*!
     * @return The number of cells in x direction
     */
    inline UInt Nx() const
    {
      return M_Nx;
    }

    //! Get grid dimension in y direction
    /*!
     * @return The number of cells in y direction
     */
    inline UInt Ny() const
    {
      return M_Ny;
    }

    //! Get grid dimension in z direction
    /*!
     * @return The number of cells in z direction
     */
    inline UInt Nz() const
    {
      return M_Nz;
    }

    //! Get zcorn data
    /*!
     * @return The vector zcorn
     */
    inline const std::vector<Real>& zcorn() const
    {
      return M_zcorn;
    }

    inline const std::vector<Real>& coord() const
    {
      return M_coord;
    }

    //! Get pillar
    /*!
     * The pillars are numbered from (1,1) to (M_Nx+1,M_Ny+1).
     * The convention adopted is the following:
     *
     *     1,2  2,2  3,2            Nx,2 Nx+1,2
     *   1,1 | 2,1| 3,1|           Nx,1|Nx+1,1
     *     | |--|-|--|-|   --   --   | |--|-|
     *     |/|  |/|  |/|             |/|  |/|
     *     |----|-|--| | --   --     |----|-|
     *     | |--|-|--|-|   --   --   | |--|-|
     *     |/   |/   |/              |/   |/
     *     |----|----|  --   --      |----|
     *
     * @param i Pillar index i
     * @param j Pillar index j
     * @return The pillar (i,j)
     */
    CPpillar pillar (const UInt& i, const UInt& j) const;

    //! Get cell
    /*!
     * The cells of the grid are numbered from (1,1,1) to (M_Nx,M_Ny,M_Nz).
     * @param i Cell index i
     * @param j Cell index j
     * @param k Cell index k
     * @return The cell (i,j,k)
     */
    CPcell cell (const UInt& i, const UInt& j, const UInt& k) const;

    void whereIs (Point3D const&, std::vector<UInt>&) const;

    void buildBB (Point3D , Point3D, Point3D, Point3D, std::vector<UInt>&) const;


    //! Get permeability field
    /*!
     * @return The permeability field
     */
    const PermeabilityField& getPermeabilityField() const;

    //! Get transmissibility field
    /*!
     * @return The transmissibility field
     */
    const TransmissibilityField& getTransmissibilityField() const;

    //@}

    //! @name Methods
    //@{

    void rescaleZ (Real, Real);

    //! Import permeability from file
    /*!
     * It reads permeability from Eclipse data file (section PERMX, PERMY, PERMZ)
     * @param filename The name of Eclipse data file
     * @param iso TRUE if the permeability field is isotropic (read only PERMX section)
     * @param scriptWithSpecialChar Set TRUE if the Eclipse data file is generated on windows
     * @return TRUE -> operation ended correctly
          FALSE -> an error occurred
     */

    bool importPermFromFile (const std::string& filename,  const bool& iso = 1,
                             const bool& scriptWithSpecialChar = 1);

    //! Export permeability to file
    /*!
     * It writes permeability in an Eclipse data file (section PERMX, PERMY, PERMZ)
     * @param filename The name of Eclipse data file
     * @param mode The stream opening mode flags
     * @return TRUE -> operation ended correctly
          FALSE -> an error occurred
     */
    bool exportPermToFile (const std::string& filename,
                           const std::ios_base::openmode& mode = std::ios_base::out | std::ios_base::trunc);

    //! Import transmissibility from file
    /*!
     * It reads transmissibility from Eclipse data file (section TRANX, TRANY, TRANZ)
     * @param filename The name of Eclipse data file
     * @param scriptWithSpecialChar Set TRUE if the Eclipse data file is generated on windows
     * @return TRUE -> operation ended correctly
          FALSE -> an error occurred
     */
    bool importTranFromFile (const std::string& filename,
                             const bool& scriptWithSpecialChar = 1);

    //! Export transmissibility to file
    /*!
     * It writes transmissibility in an Eclipse data file (section TRANX, TRANY, TRANZ)
     * @param filename The name of Eclipse data file
     * @param mode The stream opening mode flags
     * @return TRUE -> operation ended correctly
          FALSE -> an error occurred
     */
    bool exportTranToFile (const std::string& filename,
                           const std::ios_base::openmode& mode = std::ios_base::out | std::ios_base::trunc);

    //! Export grid to file
    /*!
     * It writes grid in an Eclipse data file (section SPECGRID, COORD, ZCORN, ACTNUM)
     * @param filename The name of Eclipse data file
     * @param mode The stream opening mode flags
     * @return TRUE -> operation ended correctly
          FALSE -> an error occurred
     */
    bool exportToFile (const std::string& filename,
                       const std::ios_base::openmode& mode = std::ios_base::out | std::ios_base::trunc);

    //! Export in vtk format
    /*!
     * It generates the vtk file, for the 3D visualization of the grid.
     * (Use paraview to open the vtk file)
     * @param filename The name of the vtk file created by this method
     * @return TRUE -> operation ended correctly
          FALSE -> an error occurred
     */
    bool exportVtk (const std::string& filename) const;
    bool exportVtk (const std::string& filename, std::vector<double>&, const std::string&, int) const;
    //! Display general information about the content of the class ) const;
    //! Display general information about the content of the class ) const;
    //! Display general information about the content of the class
    /*!
     * List of things displayed in the class
     * @param out Specify the output format (std::cout by default)
     */
    void showMe (std::ostream&   out = std::cout) const;

    bool readTextLine (const std::string fileName, const std::string sectionName,
                       std::string&  data, const bool specialChar = 1, int skip = 0);

    //@}
    //! To get the search tree
    /*!
      It is a 'logical const' method. In fact it may change the state by  operating
      on the mutable member G_tree_ptr.
      @note This method can be called ONLY after the construction of the CPgrid has been completed.
      @todo If you want to make thing supersafe you should add a G_tree_ptr.reset() statement at the end of all
      methods that change the state of the CPgrid.
     */
    ADT::ADTree const* searchTree() const
    {
      // Create the tree on the fly
      if (G_tree_ptr.get() == 0) G_tree_ptr.reset (new ADT::ADTree (*this) );
      return G_tree_ptr.get();
    }
  private:
    UInt M_Nx, M_Ny, M_Nz;
    std::vector<Real> M_gridBB;
    std::vector<Real> M_coord;
    std::vector<Real> M_zcorn;
    std::vector<UInt> M_actnum;
    bool M_permDefined;
    PermeabilityField M_permeability;
    bool M_tranDefined;
    TransmissibilityField M_transmissibility;
    //!For spatial search.
    /*!
      Mutable because it might be generated also on constant objects.
    */
    mutable std::auto_ptr<ADT::ADTree> G_tree_ptr;

  };


} // namespace Geometry

#endif /* GEOMCPGRID_HPP_ */
