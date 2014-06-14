/*!
 * @file Fracture3D.hpp
 * @brief Class that handles the fractures.
 */

#ifndef FRACTURE3D_HPP_
#define FRACTURE3D_HPP_

#include "core/TypeDefinition.hpp"
#include <fstream>

namespace FVCode3D
{

class Mesh3D;

//! Class that implements a fracture.
/*!
	@class Fracture3D
    This class implements the concept of a fracture as a vector of facets of a Mesh3D.
 */
class Fracture3D{
public:
	//! @name Constructor & Destructor
	//@{

	//! Constructor
	/*!
	 * @param mesh reference of a Mesh3D
	 */
	Fracture3D(const Mesh3D & mesh);

	//! Copy constructor
	/*!
	 * @param f Fracture3D
	 */
	Fracture3D(const Fracture3D & f);

	//! Constructor, getting the vector of fracture facets ids
	/*!
	 * @param mesh reference of a Mesh3D
	 * @param fractureFacets The vector of facets ids that define the fracture
	 * @param id identifier of the fracture
	 */
	Fracture3D(const Mesh3D & mesh, const std::vector<UInt> & fractureFacets, const UInt id);

	//! Destructor
	~Fracture3D() = default;

	//@}

	//! @name Operators
	//@{

	//! Assignment operator
	/*!
	 * @param f Fracture3D
	 * @return this Fracture3D
	 */
	Fracture3D & operator=(Fracture3D f)
		{ M_id = f.getId(); M_fractureFacets = f.getFractureFacetsId(); return *this; };

	//@}

	//! @name Get Methods
	//@{

	//! Get the id of the fracture (const)
	/*!
	 * @return the id of the fracture
	 */
	const UInt & getId() const
		{ return M_id; }

	//! Get the id of the fracture
	/*!
	 * @return reference to the id of the fracture
	 */
	UInt & getId()
		{ return M_id; }

	//! Get the vector of the fracture facets ids (const)
	/*!
	 * @return A constant reference to the vector of fracture facets ids
	 */
	const std::vector<UInt> & getFractureFacetsId() const
		{ return M_fractureFacets; }

	//! Get the vector of the fracture facets ids
	/*!
	 * @return A reference to the vector of fracture facets ids
	 */
	std::vector<UInt> & getFractureFacetsId()
		{ return M_fractureFacets; }

	//! Get the number of facets used for the fracture discretization (const)
	/*!
	 * @return The number of facets used for the fracture discretization
	 */
	UInt getNumberOfFractureFacets() const
		{ return M_fractureFacets.size(); }

	//! Get a constant reference to the 3D mesh (const)
	/*!
	 * @return The 3D mesh
	 */
	const Mesh3D & getMesh() const
		{ return M_mesh; }

	//@}

	//! @name Set Methods
	//@{

	//! Insert the id of a fracture facet
	/*!
	 * @param idFacet id of a fracture facet
	 */
	void push_back(const UInt idFacet)
		{ return M_fractureFacets.push_back(idFacet); }

	//@}

	//! @name Methods
	//@{

	//! Export in vtk format the facets representing the fracture (as VTK_POLYGON)
	/*!
	 * (Use paraview to open the vtk file)
	 * @param filename The name of the vtk file created by this method
	 * @return TRUE  -> operation ended correctly
			   FALSE -> an error occurred
	 */
	bool exportVTK(const std::string & filename) const;

	//! Display general information about the content of the class
	/*!
	 * @param out Specify the output format (std::cout by default)
	 */
	void showMe(std::ostream & out=std::cout) const;

	//@}

private:

	//! Id of the fracture
	UInt M_id;
	//! It stores the facets id used for the discretization
	std::vector<UInt> M_fractureFacets;
	//! Reference to the mesh
	const Mesh3D & M_mesh;

};

} // namespace FVCode3D

#endif /* FRACTURE3D_HPP_ */
