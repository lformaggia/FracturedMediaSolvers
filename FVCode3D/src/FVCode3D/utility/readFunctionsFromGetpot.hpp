/*
 * readFunctionsFromGetpot.hpp
 *
 *  Created on: Oct 25, 2018
 *      Author: forma
 */

#ifndef SUBPROJECTS__FVCODE3D_SRC_FVCODE3D_UTILITY_READFUNCTIONSFROMGETPOT_HPP_
#define SUBPROJECTS__FVCODE3D_SRC_FVCODE3D_UTILITY_READFUNCTIONSFROMGETPOT_HPP_

#include <FVCode3D/boundaryCondition/BC.hpp>
#include <GetPot>
#include <vector>
#include <string>
#include <FVCode3D/utility/functionLoader.hpp>
namespace FVCode3D
{
  /*
   * enum BCType
{
    Dirichlet=0,
    Neumann=1
};

enum BorderLabel : UInt
{
    Internal    = 0,
    Left        = 1,
    Right       = 2,
    Front       = 3,
    Back        = 4,
    Bottom      = 5,
    Top         = 6,
    FractureTip = std::numeric_limits<UInt>::max()
};
   */
  class ReadFunctionsFromGetpot
  {
  public:
    //! The class is default constructible
    ReadFunctionsFromGetpot()=default;

    //! We take directly the info from the getpot data associated to a getpot file
    /*!
     *  The file should contain the following fields
     *  functionLibrary=libraryName  (example libfunction.so)
     *  [bc]
     *  [bc/Left]
     *  type= 0 or 1 (0 Dirichlet , 1 Neumann)
     *  f=FunctioName
     *  the same for Right,Front, Back,Bottom,Top
     *  If type is missing Dirichlet is assumed
     *
     * If the library is missing is not an error, but only the in-built
     * functions fZero and fOne will be available
     *  @par gp A GetPot data object
     *   */
    ReadFunctionsFromGetpot(GetPot const & gp);
    //! I can give the GetPot file name directly
    /*!
     * @par fileName a valid GetPot file
     *
     */
    ReadFunctionsFromGetpot(std::string fileName):ReadFunctionsFromGetpot(GetPot(fileName)){};
    //! Sets the getpot data and loads the library
    void setGetpotData(GetPot const & gp);

    //! Loads the source function
    /*
     * The GetPot data should contain
     * [problem]
     * sourceFunction=functionName
     *
     * if not present fZero is assumed
     *
     * @return a function wrapped in a std::function object
     */
    Func getSourceFunction() const;

    //! Loads library parsing the getpot data
    /*
     *  The getpot data  should contain the following fields
     *  functionLibrary=libraryName  (example libfunction.so)
     *  [bc]
     *  [bc/Left]
     *  type= 0 or 1 (0 Dirichlet , 1 Neumann)
     *  f=FunctioName
     *  the same for Right,Front, Back,Bottom,Top
     *  If type is missing Dirichlet is assumed
     *
     *
     * the field functionLibrary=libraryname
     * should be present. If not only fZero and fOne will be available
     *
     * @par gp a A getPot data object
     */
    std::vector<BoundaryConditions::BorderBC> getBCForStandardDomain() const;
    /*
    //! Releases possible shared library attached
    void release(){this->M_functionLoader.release();}
    //! The destructor releases possible library
    ~ReadFunctionsFromGetpot()
    {
      this->release();
    }*/
 private:
    //! Loads library with function
    /*!
     * The getpot data should contain
     * functionLibrary=libname
     * @par gp the getopot data object
     */
    void loadLib(GetPot const & gp);
   //! Loads library given name
    /*!
     * @par libFileName the name of the library
     */
    void loadLib(std::string const & libFileName);
    FunctionLoader M_functionLoader;
    GetPot M_gp;
    bool invalid=true;
    
    const static std::vector<std::pair<std::string,BorderLabel> > bcList;
    //     {"Internal",BorderLabel::Internal}, // not needed here
    const static std::vector<std::pair<std::string,BCType> > bcType;
    static constexpr UInt noType=99999;
};
  

}

#endif /* SUBPROJECTS__FVCODE3D_SRC_FVCODE3D_UTILITY_READFUNCTIONSFROMGETPOT_HPP_ */
