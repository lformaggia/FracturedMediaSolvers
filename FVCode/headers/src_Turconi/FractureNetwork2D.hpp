 /*!
 *	@file Fracture2D.hpp
 *	@brief Class for fracture network in 2D space.
 *
 *	@author Luca Turconi <luca.turconi@moxoff.com>
 *
 */ 
 
 
#ifndef FRACTURENETWORK2D_HPP_
#define FRACTURENETWORK2D_HPP_

#include "TypeDefinition.hpp"
#include "Fracture2D.hpp"

namespace Geometry{

	/*!
		@class FractureNetwork2D
    	This class implements the concept of fracture network in 2D space.
    */
class FractureNetwork2D{
public:
	//! @name Constructor & Destructor
	//@{
	
	//! Empty constructor
	FractureNetwork2D();
	
	//! Copy constructor
	FractureNetwork2D(const Geometry::FractureNetwork2D & fn);

	//! Constructor, getting two points
	/*!
	 * @param fabFileName The name of the fab file containing the network description
	 */
	FractureNetwork2D(const std::string & fabFileName);
	
	//! Destructor
	~FractureNetwork2D();
	
	//@}
	
	//! @name Get Methods
	//@{
	
	//! Get the vector storing the fractures of the network (const)
	/*!
	 * @return A constant reference to M_fractureNetwork
	 */
	const std::vector<Geometry::Fracture2D> & getNetwork() const
		{ return M_fractureNetwork; }

	//! Get the vector storing the fractures of the network
	/*!
	 * @return A reference to M_fractureNetwork
	 */
	std::vector<Geometry::Fracture2D> & getNetwork()
		{ return M_fractureNetwork; }
	
	//@}
	
	//! @name Methods
	//@{
	
	//! Test if the fractures are defined using meter as unit of measurement
	/*!
	 * @return TRUE -> Fractures defined using meter
	 		   FALSE -> Fractures defined using another unit
	 */
	bool isMetric() const
		{ return M_isMetric; }
	
	//! Insert a Fracture2D in the Fracture Network
	/*!
	 * @param f The fracture to be inserted
	 */
	void push_back(const Geometry::Fracture2D & f)
		{ M_fractureNetwork.push_back(f); }
	
	//! The number of fractures stored in the network
	/*!
	 * @return The number of fractures stored in the network
	 */
	UInt size() const
		{ return M_fractureNetwork.size(); }
	
	//! Compute the distance between the fracture network and the point p
	/*!
	 * @param p The point
	 * @return The distance between this and p
	 */
	Real distance(const Geometry::Point2D & p) const;
	
	//! Build the discretization of all the fractures in the network using elements of dimension as similar as possible to h
	/*!
	 * @param h The desidered elements dimension
	 */
	void buildFractureNetworkDiscretization(const Real & h);
	
	//! Build the discretization of all the fractures in the network using a fixed number of elements
	/*!
	 * @param Nelements The numeber of elements to be used for the discretization of each fracture
	 */
	void buildFractureNetworkDiscretization(const UInt & Nelements);
	
	//! Export in different vtk file the segments representing each fracture in the network
	/*!
	 * (Use paraview to open the vtk file)
	 * @param filename The prefix name of the vtk files created by this method
	 * @return TRUE -> operation ended correctly
			   FALSE -> an error occurred
	 */
	bool exportVtk(const std::string & filename="FractureNetwork") const;
	
	//! Export in a single vtk file all the segments representing each fracture in the network
	/*!
	 * (Use paraview to open the vtk file)
	 * @param filename The name of the vtk file created by this method
	 * @return TRUE -> operation ended correctly
			   FALSE -> an error occurred
	 */
	bool exportNetworkVtk(const std::string & fileName="FractureNetwork.vtk") const;
	
	//! Export in different vtk file the set of the nodes used for the discretization of each fracture in the network
	/*!
	 * (Use paraview to open the vtk file)
	 * @param filename The prefix name of the vtk files created by this method
	 * @return TRUE -> operation ended correctly
			   FALSE -> an error occurred
	 */
	bool exportFractureDiscretizationVtk(const std::string & filename="FractureNetwork") const;
	
	//! Display general information about the content of the class
	/*!
	 * @param out Specify the output format (std::cout by default)
	 */
	void showMe(std::ostream & out=std::cout) const;
	
private:

	//! It stores the the fractures of the network
	std::vector<Geometry::Fracture2D> M_fractureNetwork;
	//! Unit of measurement: TRUE for meter
	bool M_isMetric;
};


} // namespace Geometry


namespace FabFile{

/*!
	@fn readProperties
   	This function allows to read the fracture properties contained in a string extracted from a .fab data file.
   	@param line The string containing the data
   	@param n_points The number of points used to define the fracture (in 3D space)
   	@param perm A reference to the variable in which the function store the fracture permeability
   	@param compr A reference to the variable in which the function store the fracture compressibility
   	@param aperture A reference to the variable in which the function store the fracture aperture
    */
void readProperties( std::string line, UInt & n_punti,
					 Real & perm, Real & compr, Real & aperture );
/*!
	@fn readPoints
   	This function allows to read the coordinates of a single fracture point contained in a string extracted from a .fab data file.
   	@param line The string containing the data
   	@param px A reference to the variable in which the function store the x coordinate of the cosidered point
   	@param py A reference to the variable in which the function store the y coordinate of the cosidered point
   	@param pz A reference to the variable in which the function store the z coordinate of the cosidered point
    */		 
void readPoints( std::string line, Real & px, Real & py, Real & pz );

} // namespace FabFile


#endif /* FRACTURENETWORK2D_HPP_ */
