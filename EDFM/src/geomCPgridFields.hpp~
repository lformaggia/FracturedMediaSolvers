 /*!
 *	@file geomCPgridFields.hpp
 *	@brief Classes for permeability and transmissibility field on Corner Point Grid.
 *
 *	@author Luca Turconi <lturconi@gmail.com>
 *  @date 22-09-2012
 *
 */ 

#ifndef GEOMCPGRIDFIELDS_HPP_
#define GEOMCPGRIDFIELDS_HPP_

#include<vector>
#include<iostream>
#include<string>
#include<fstream>

#include "TypeDefinition.hpp"
#include "geomPoint3D.hpp"


namespace Geometry
{
class CPgrid;	// Forward declaration

/*!
		@class PermeabilityField
		
		@author Luca Turconi <lturconi@gmail.com>
		
    	This class implements the concept of permeability field.
    	
    	Field informations are stored in vectors: M_permx, M_permy, M_permz.
    	
    	The constructor build an empty field. It sets only the grid dimensions.
    	To fill an object with value of permeability use the method importFromFile.
    	
    	The bool variable M_isotropic is TRUE for isotropic fields.
    	
    	A tools of different methods to extract permeability data are available.
    	
    	The method showMe allows to print the class informations.
    	
    */
class PermeabilityField{
public:
	//! @name Constructor & Destructor
	//@{
		
	//! Empty Constructor
	PermeabilityField();
	
	//! Constructor, getting field dimensions
	/*!
	 * @param nx The x-dimension
	 * @param ny The y-dimension
	 * @param nz The z-dimension
	 * @param iso TRUE if the permeability field is isotropic
	 */
	PermeabilityField(const UInt & nx, const UInt & ny, const UInt & nz, const bool & iso=1);
	
	//! Constructor, getting associated grid
	/*!
	 * @param grid The associated Corner Point grid
	 * @param iso TRUE if the permeability field is isotropic
	 */
	PermeabilityField(const CPgrid & grid, const bool & iso=1);
	
	//! Destructor
	~PermeabilityField();
	
	//@}
	
	//! @name Get Methods
	//@{
		
	//! Get value from permx vector
	/*!
	 * @param i The id of the value in permx vector
	 * @return The permeability value
	 */
	inline Real & getPermxElement(const UInt & i)
		{ return M_permx[i]; };
	
	//! Get value from permx vector (const)
	/*!
	 * @param i The id of the value in permx vector
	 * @return The permeability value
	 */
	inline Real getPermxElement(const UInt & i) const
		{ return M_permx[i]; };
		
	//! Get value from permy vector
	/*!
	 * @param i The id of the value in permy vector
	 * @return The permeability value
	 */
	inline Real & getPermyElement(const UInt & i)
		{ return M_isotropic ? M_permx[i] : M_permy[i]; };
	
	//! Get value from permy vector (const)
	/*!
	 * @param i The id of the value in permy vector
	 * @return The permeability value
	 */
	inline Real getPermyElement(const UInt & i) const
		{ return M_isotropic ? M_permx[i] : M_permy[i]; };
		
	//! Get value from permz vector
	/*!
	 * @param i The id of the value in permz vector
	 * @return The permeability value
	 */
	inline Real & getPermzElement(const UInt & i)
		{ return M_isotropic ? M_permx[i] : M_permz[i]; };
	
	//! Get value from permz vector (const)
	/*!
	 * @param i The id of the value in permz vector
	 * @return The permeability value
	 */
	inline Real getPermzElement(const UInt & i) const
		{ return M_isotropic ? M_permx[i] : M_permz[i]; };
	
	//! Get permx vector
	/*!
	 * @return A reference to permx vector
	 */
	inline std::vector<Real> & getPermxVec()
		{ return M_permx; }
	
	//! Get permx vector (const)
	/*!
	 * @return A reference to permx vector
	 */
	inline const std::vector<Real> & getPermxVec() const
		{ return M_permx; }
		
	//! Get permy vector
	/*!
	 * @return A reference to permy vector
	 */
	inline std::vector<Real> & getPermyVec()
		{ return M_isotropic ? M_permx : M_permy; }
		
	//! Get permy vector (const)
	/*!
	 * @return A reference to permy vector
	 */
	inline const std::vector<Real> & getPermyVec() const
		{ return M_isotropic ? M_permx : M_permy; }
		
	//! Get permz vector
	/*!
	 * @return A reference to permz vector
	 */
	inline std::vector<Real> & getPermzVec()
		{ return M_isotropic ? M_permx : M_permz; }
	
	//! Get permz vector (const)
	/*!
	 * @return A reference to permz vector
	 */
	inline const std::vector<Real> & getPermzVec() const
		{ return M_isotropic ? M_permx : M_permz; }
    
    //! Get permx of given cell
	/*!
	 * @param i The cell index i
	 * @param j The cell index j
	 * @param k The cell index k
	 * @return A reference to permx value of given cell
	 */
    inline Real & getCellPermX(const UInt & i, const UInt & j, const UInt & k)
		{ return M_permx[i+j*M_Nx+k*M_Nx*M_Ny]; }
	
	//! Get permx of given cell (const)
	/*!
	 * @param i The cell index i
	 * @param j The cell index j
	 * @param k The cell index k
	 * @return The permx value of given cell
	 */
	inline Real getCellPermX(const UInt & i, const UInt & j, const UInt & k) const
		{ return M_permx[i+j*M_Nx+k*M_Nx*M_Ny]; }
	
	//! Get permy of given cell
	/*!
	 * @param i The cell index i
	 * @param j The cell index j
	 * @param k The cell index k
	 * @return A reference to permy value of given cell
	 */
    inline Real & getCellPermY(const UInt & i, const UInt & j, const UInt & k)
		{ return M_isotropic ? M_permx[i+j*M_Nx+k*M_Nx*M_Ny] : M_permy[i+j*M_Nx+k*M_Nx*M_Ny]; }
	
	//! Get permy of given cell (const)
	/*!
	 * @param i The cell index i
	 * @param j The cell index j
	 * @param k The cell index k
	 * @return The permy value of given cell
	 */
	inline Real getCellPermY(const UInt & i, const UInt & j, const UInt & k) const
		{ return M_isotropic ? M_permx[i+j*M_Nx+k*M_Nx*M_Ny] : M_permy[i+j*M_Nx+k*M_Nx*M_Ny]; }
    
    //! Get permz of given cell
	/*!
	 * @param i The cell index i
	 * @param j The cell index j
	 * @param k The cell index k
	 * @return A reference to permz value of given cell
	 */
    inline Real & getCellPermZ(const UInt & i, const UInt & j, const UInt & k)
		{ return M_isotropic ? M_permx[i+j*M_Nx+k*M_Nx*M_Ny] : M_permz[i+j*M_Nx+k*M_Nx*M_Ny]; }
	
	//! Get permz of given cell (const)
	/*!
	 * @param i The cell index i
	 * @param j The cell index j
	 * @param k The cell index k
	 * @return The permz value of given cell
	 */
	inline Real getCellPermZ(const UInt & i, const UInt & j, const UInt & k) const
		{ return M_isotropic ? M_permx[i+j*M_Nx+k*M_Nx*M_Ny] : M_permz[i+j*M_Nx+k*M_Nx*M_Ny]; }
    
    //! Get permeability values of given cell
	/*!
	 * @param i The cell index i
	 * @param j The cell index j
	 * @param k The cell index k
	 * @return A 3D vector containing the permeability values (permx,permy,permz) of given cell
	 */
	Vector3D getCellPermVec(const UInt & i, const UInt & j, const UInt & k) const;

	//@}
	
	//! @name Set Methods
	//@{
	
	//! Set permx vector
	/*!
	 * @param p The new permx vector
	 */
	void setPermx(const std::vector<Real> & p);
	
	//! Set permy vector
	/*!
	 * @param p The new permy vector
	 */
	void setPermy(const std::vector<Real> & p);
	
	//! Set permz vector
	/*!
	 * @param p The new permz vector
	 */
	void setPermz(const std::vector<Real> & p);
	
	//! Set field dimensions
	/*!
	 * @param nx The x-dimension
	 * @param ny The y-dimension
	 * @param nz The z-dimension
	 */
	inline void setDim(const UInt nx, const UInt ny, const UInt nz)
		{ M_Nx=nx; M_Ny=ny; M_Nz=nz; }
		
	//! Set field isotropy
	/*!
	 * @param iso TRUE if the permeability field is isotropic
	 */
	inline void setIsotropy(const bool & iso=1)
		{ M_isotropic=iso; }
	
	//@}
	
	//! @name Methods
	//@{
	
	//! Import permeability from file
	/*!
	 * It reads permeability from Eclipse data file (section PERMX, PERMY, PERMZ)
	 * @param filename The name of Eclipse data file
	 * @param iso TRUE if the permeability field is isotropic (read only PERMX section)
	 * @param scriptWithSpecialChar Set TRUE if the Eclipse data file is generated on windows
	 * @return TRUE -> operation ended correctly
				FALSE -> an error occurred
	 */
	bool importFromFile(const std::string & filename, const bool & iso=1,
						const bool & scriptWithSpecialChar=1);
	
	//! Export permeability to file
	/*!
	 * It writes permeability in an Eclipse data file (section PERMX, PERMY, PERMZ)
	 * @param filename The name of Eclipse data file
	 * @param mode The stream opening mode flags
	 * @return TRUE -> operation ended correctly
				FALSE -> an error occurred
	 */
	bool exportToFile(const std::string & filename,
			const std::ios_base::openmode & mode = std::ios_base::out | std::ios_base::trunc);
	
	//! Display general information about the content of the class
	/*!
	 * List of things displayed in the class
	 * @param out Specify the output format (std::cout by default)
	 */
	void showMe(std::ostream & out=std::cout) const;
	
	//@}
	
private:
	UInt M_Nx,M_Ny,M_Nz;
	bool M_isotropic;
    std::vector<Real> M_permx;
    std::vector<Real> M_permy;
    std::vector<Real> M_permz;
};
	

/*!
		@class TransmissibilityField
		
		@author Luca Turconi <lturconi@gmail.com>
		
    	This class implements the concept of transmissibility field.
    	
    	Field informations are stored in vectors: M_tranx, M_trany, M_tranz.
    	
    	The constructor build an empty field. It sets only the grid dimensions.
    	To fill an object with value of transmissibility use the method importFromFile.
    	
    	A tools of different methods to extract transmissibility data are available.
    	
    	The method showMe allows to print the class informations.
    	
    */
class TransmissibilityField{
public:
	//! @name Constructor & Destructor
	//@{
		
	//! Empty Constructor
	TransmissibilityField();
	
	//! Constructor, getting field dimensions
	/*!
	 * @param nx The x-dimension
	 * @param ny The y-dimension
	 * @param nz The z-dimension
	 */
	TransmissibilityField(const UInt & nx, const UInt & ny, const UInt & nz);
	
	//! Constructor, getting associated grid
	/*!
	 * @param grid The associated Corner Point grid
	 */
	TransmissibilityField(const CPgrid & grid);
	
	//! Destructor
	~TransmissibilityField();
	
	//@}
	
	//! @name Get Methods
	//@{
		
	//! Get value from tranx vector
	/*!
	 * @param i The id of the value in tranx vector
	 * @return The transmissibility value
	 */
	inline Real & getTranxElement(const UInt & i)
		{ return M_tranx[i]; };
		
	//! Get value from tranx vector (const)
	/*!
	 * @param i The id of the value in tranx vector
	 * @return The transmissibility value
	 */
	inline Real getTranxElement(const UInt & i) const
		{ return M_tranx[i]; };
		
	//! Get value from trany vector
	/*!
	 * @param i The id of the value in trany vector
	 * @return The transmissibility value
	 */
	inline Real & getTranyElement(const UInt & i)
		{ return M_trany[i]; };
		
	//! Get value from trany vector (const)
	/*!
	 * @param i The id of the value in trany vector
	 * @return The transmissibility value
	 */
	inline Real getTranyElement(const UInt & i) const
		{ return M_trany[i]; };
	
	//! Get value from tranz vector
	/*!
	 * @param i The id of the value in tranz vector
	 * @return The transmissibility value
	 */
	inline Real & getTranzElement(const UInt & i)
		{ return M_tranz[i]; };
		
	//! Get value from tranz vector (const)
	/*!
	 * @param i The id of the value in tranz vector
	 * @return The transmissibility value
	 */
	inline Real getTranzElement(const UInt & i) const
		{ return M_tranz[i]; };
	
	//! Get tranx vector
	/*!
	 * @return A reference to tranx vector
	 */
	inline std::vector<Real> & getTranxVec()
		{ return M_tranx; }
	
	//! Get tranx vector (const)
	/*!
	 * @return A reference to tranx vector
	 */
	inline const std::vector<Real> & getTranxVec() const
		{ return M_tranx; }
		
	//! Get trany vector
	/*!
	 * @return A reference to trany vector
	 */
	inline std::vector<Real> & getTranyVec()
		{ return M_trany; }
	
	//! Get trany vector (const)
	/*!
	 * @return A reference to trany vector
	 */
	inline const std::vector<Real> & getTranyVec() const
		{ return M_trany; }
	
	//! Get tranz vector
	/*!
	 * @return A reference to tranz vector
	 */
	inline std::vector<Real> & getTranzVec()
		{ return M_tranz; }
	
	//! Get tranz vector (const)
	/*!
	 * @return A reference to tranz vector
	 */
	inline const std::vector<Real> & getTranzVec() const
		{ return M_tranz; }
	
	//! Get tranx of given cell
	/*!
	 * @param i The cell index i
	 * @param j The cell index j
	 * @param k The cell index k
	 * @return A reference to tranx value of given cell
	 */
	inline Real & getCellTranX(const UInt & i, const UInt & j, const UInt & k)
		{ return M_tranx[i+j*M_Nx+k*M_Nx*M_Ny]; }
	
	//! Get tranx of given cell (const)
	/*!
	 * @param i The cell index i
	 * @param j The cell index j
	 * @param k The cell index k
	 * @return The tranx value of given cell
	 */
	inline Real getCellTranX(const UInt & i, const UInt & j, const UInt & k) const
		{ return M_tranx[i+j*M_Nx+k*M_Nx*M_Ny]; }
	
	//! Get trany of given cell
	/*!
	 * @param i The cell index i
	 * @param j The cell index j
	 * @param k The cell index k
	 * @return A reference to trany value of given cell
	 */
	inline Real & getCellTranY(const UInt & i, const UInt & j, const UInt & k)
		{ return M_trany[i+j*M_Nx+k*M_Nx*M_Ny]; }
	
	//! Get trany of given cell (const)
	/*!
	 * @param i The cell index i
	 * @param j The cell index j
	 * @param k The cell index k
	 * @return The trany value of given cell
	 */
	inline Real getCellTranY(const UInt & i, const UInt & j, const UInt & k) const
		{ return M_trany[i+j*M_Nx+k*M_Nx*M_Ny]; }
	
	//! Get tranz of given cell
	/*!
	 * @param i The cell index i
	 * @param j The cell index j
	 * @param k The cell index k
	 * @return A reference to tranz value of given cell
	 */
	inline Real & getCellTranZ(const UInt & i, const UInt & j, const UInt & k)
		{ return M_tranz[i+j*M_Nx+k*M_Nx*M_Ny]; }
	
	//! Get tranz of given cell (const)
	/*!
	 * @param i The cell index i
	 * @param j The cell index j
	 * @param k The cell index k
	 * @return The tranz value of given cell
	 */
	inline Real getCellTranZ(const UInt & i, const UInt & j, const UInt & k) const
		{ return M_tranz[i+j*M_Nx+k*M_Nx*M_Ny]; }
	
	//! Get transmissibility values of given cell
	/*!
	 * @param i The cell index i
	 * @param j The cell index j
	 * @param k The cell index k
	 * @return A 3D vector containing the transmissibility values (tranx,trany,tranz) of given cell
	 */
	Vector3D getCellTranVec(const UInt & i, const UInt & j, const UInt & k) const;
	
	//@}
	
	//! @name Set Methods
	//@{
	
	//! Set tranx vector
	/*!
	 * @param p The new tranx vector
	 */
	void setTranx(const std::vector<Real> & p);
	
	//! Set trany vector
	/*!
	 * @param p The new trany vector
	 */
	void setTrany(const std::vector<Real> & p);
	
	//! Set tranz vector
	/*!
	 * @param p The new tranz vector
	 */
	void setTranz(const std::vector<Real> & p);
	
	//! Set field dimensions
	/*!
	 * @param nx The x-dimension
	 * @param ny The y-dimension
	 * @param nz The z-dimension
	 */
	inline void setDim(const UInt & nx, const UInt & ny, const UInt & nz)
		{ M_Nx=nx; M_Ny=ny; M_Nz=nz; }
	
	//@}
	
	//! @name Methods
	//@{
	
	//! Import transmissibility from file
	/*!
	 * It reads transmissibility from Eclipse data file (section TRANX, TRANY, TRANZ)
	 * @param filename The name of Eclipse data file
	 * @param scriptWithSpecialChar Set TRUE if the Eclipse data file is generated on windows
	 * @return TRUE -> operation ended correctly
				FALSE -> an error occurred
	 */
	bool importFromFile(const std::string & filename,
						const bool & scriptWithSpecialChar=1);
	
	//! Export transmissibility to file
	/*!
	 * It writes transmissibility in an Eclipse data file (section TRANX, TRANY, TRANZ)
	 * @param filename The name of Eclipse data file
	 * @param mode The stream opening mode flags
	 * @return TRUE -> operation ended correctly
				FALSE -> an error occurred
	 */
	bool exportToFile(const std::string & filename,
			const std::ios_base::openmode & mode = std::ios_base::out | std::ios_base::trunc);
	
	//! Display general information about the content of the class
	/*!
	 * List of things displayed in the class
	 * @param out Specify the output format (std::cout by default)
	 */
	void showMe(std::ostream & out=std::cout) const;
	
	//@}
	
private:
	UInt M_Nx,M_Ny,M_Nz;
    std::vector<Real> M_tranx;
    std::vector<Real> M_trany;
    std::vector<Real> M_tranz;
};

} // namespace Geometry

#endif /* GEOMCPGRIDFIELDS_HPP_ */