/*!
 * @file FractureNetwork3D.hpp
 * @brief Class that handles the fracture network.
 */

#ifndef FRACTURENETWORK3D_HPP_
#define FRACTURENETWORK3D_HPP_

#include "core/TypeDefinition.hpp"

namespace FVCode3D
{

class Mesh3D;
class Fracture3D;

//! Class that handles the fracture network
/*!
 * @class FractureNetwork3D
 * This class implements the concept of fracture network as a vector of Fracture3D.
 */
class FractureNetwork3D
{
public:
	//! @name Constructor & Destructor
	//@{

	//! Constructor
	/*!
	 * @param mesh reference to a Mesh3D
	 */
	FractureNetwork3D(const Mesh3D & mesh);

	//! Copy constructor
	/*!
	 * @param fn reference to a FractureNetwork3D
	 */
	FractureNetwork3D(const FractureNetwork3D & fn);

	//! Constructor from a vector of fractures
	/*!
	 * @param mesh reference to a Mesh3D
	 * @param fractures vector of fractures (Fracture3D)
	 */
	FractureNetwork3D(const Mesh3D & mesh, const std::vector<Fracture3D> & fractures);

	//! Destructor
	~FractureNetwork3D() = default;

	//@}

	//! @name Get Methods
	//@{

	//! Get the vector that store the fractures of the network (const)
	/*!
	 * @return A constant reference to the vector of fractures
	 */
	const std::vector<Fracture3D> & getNetwork() const
		{ return M_fractureNetwork; }

	//! Get the vector that store the fractures of the network
	/*!
	 * @return A reference to the vector of fractures
	 */
	std::vector<Fracture3D> & getNetwork()
		{ return M_fractureNetwork; }

	//! Get the i-th fracture (const)
	/*!
	 * @param i id of the fracture
	 * @return A the i-th fracture
	 */
	const Fracture3D & getFracture(const UInt i) const
		{ return M_fractureNetwork[i]; }

	//! Get the Mesh3D
	/*!
	 * @return The 3D mesh
	 */
	const Mesh3D & getMesh() const
		{ return M_mesh; }

	//@}

	//! @name Methods
	//@{

	//! Insert a vector of Fracture3D in the Fracture Network
	/*!
	 * @param fractures The fractures to insert
	 */
	void addFractures(std::vector<Fracture3D> & fractures);

	//! The number of fractures stored in the network
	/*!
	 * @return The number of fractures stored in the network
	 */
	UInt size() const
		{ return M_fractureNetwork.size(); }

	//! Export in different vtk file the segments representing each fracture in the network
	/*!
	 * (Use paraview to open the vtk file)
	 * @param filename The prefix name of the vtk files created by this method
	 * @return TRUE  -> operation ended correctly
	 *		   FALSE -> an error occurred
	 */
	bool exportVTK(const std::string & prefixFileName="FractureNetwork") const;

	//! Export in a single vtk file all the segments representing each fracture in the network
	/*!
	 * (Use paraview to open the vtk file)
	 * @param filename The name of the vtk file created by this method
	 * @return TRUE  -> operation ended correctly
	 *		   FALSE -> an error occurred
	 */
	bool exportNetworkVTK(const std::string & filename="FractureNetwork.vtk") const;

	//! Display general information about the content of the class
	/*!
	 * @param out Specify the output format (std::cout by default)
	 */
	void showMe(std::ostream & out=std::cout) const;

private:

	//! It stores the the fractures of the network
	std::vector<Fracture3D> M_fractureNetwork;
	//! Reference to the mesh
	const Mesh3D & M_mesh;

};

} // namespace FVCode3D

#endif /* FRACTURENETWORK3D_HPP_ */
