/*!
 *	@file import.hpp
 *	@brief Classes for loading files.
 */

#ifndef IMPORT_HPP_
#define IMPORT_HPP_

#include "core/TypeDefinition.hpp"

class Mesh3D;
class PropertiesMap;

//! Class used to read files.
/*!
 * @class Importer
 * This is a base abstract class that allows to import files.
 */
class Importer
{
public:

	//! Constructor
	/*!
	 * Constructor from a grid file
	 * @param filename filename of the file
	 * @param mesh reference to a Geometry::Mesh3D
	 * @param properties reference to a Geometry::PropertiesMap
	 */
	Importer(const std::string filename, Geometry::Mesh3D & mesh, Geometry::PropertiesMap & properties):
		M_filename(filename), M_mesh(mesh), M_properties(properties) {}

	//! Import from a grid file
	/*!
	 * @param fracturesOn if true, imports the fractures, else the fractures are disabled
	 */
	virtual void import(bool fracturesOn) = 0;

	//! Generate the BC ids from a standard TPFA file format
	/*!
	 * Add the BCs ids to the boundary facets.
	 * The id is set by considering the maximum component and the sign of the normal of a boundary facet.
	 * @param theta rotation angle along z-axis. It is used only to compute the BC ids. Default = 0
	 * @pre import the file
	 */
	virtual void extractBC(const Real theta = 0.);
	
	//! Add the fractures network
	/*!
	 * Add the fracture network
	 * @pre import the file
	 */
	virtual void addFractures() = 0;
	
	//! Add the lacking data
	/*!
	 * Add the BCs and the fracture network
	 * @param theta rotation angle along z-axis. It is used only to compute the BC ids. Default = 0
	 * @pre import the file
	 */
	virtual void addBCAndFractures(const Real theta = 0.);
	
	//! Get filename
	/*!
	 * @return the filename
	 */
	const std::string getFileName() const
		{ return M_filename; }

	//! Get mesh
	/*!
	 * @return reference to the mesh
	 */
	const Geometry::Mesh3D & getMesh() const
		{ return M_mesh; }

	//! Get properties
	/*!
	 * @return reference to the properties
	 */
	const Geometry::PropertiesMap & getProperties() const
		{ return M_properties; }

	//! Destructor
	virtual ~Importer() = default;

protected:

	//! Filename
	std::string M_filename;
	//! Reference to a Geometry::Mesh3D
	Geometry::Mesh3D & M_mesh;
	//! Reference to a Geometry::PropertiesMap
	Geometry::PropertiesMap & M_properties;

private:

	//! No default constructor
	Importer();

	//! No copy-constructor
	Importer(const Importer &);

	//! No assignment operator
	Importer & operator=(const Importer &);

};

//! Class used to read a Medit format file (.mesh).
/*!
 * @class ImporterMedit
 * This class allows to read a Medit format file (.mesh).
 */
class ImporterMedit : public Importer
{
public:

	//! Constructor
	/*!
	 * Constructor from a medit file
	 * @param filename filename of the .mesh file
	 * @param mesh reference to a Geometry::Mesh3D
	 * @param properties reference to a Geometry::PropertiesMap
	 */
	ImporterMedit(const std::string filename, Geometry::Mesh3D & mesh, Geometry::PropertiesMap & properties):
		Importer(filename, mesh, properties) {}

	//! Import from a .mesh file
	/*!
	 * Read points, polygons, polyhedra
	 * @param fracturesOn if true, imports the fractures, else the fractures are disabled
	 */
	virtual void import(bool fracturesOn = true);
	
	//! Add the fractures network
	/*!
	 * Add the fracture network
	 * @pre import the file
	 */
	virtual void addFractures();
	
	//! Destructor
	virtual ~ImporterMedit() = default;

private:

	//! No default constructor
	ImporterMedit();

	//! No copy-constructor
	ImporterMedit(const ImporterMedit &);

	//! No assignment operator
	ImporterMedit & operator=(const ImporterMedit &);
};

//! Class used to read a standard TPFA format file (.grid).
/*!
 * @class ImporterTPFA
 * This class allows to read a standard TPFA format file (.grid).
 */
class ImporterTPFA : public Importer
{
public:

	//! Constructor
	/*!
	 * Constructor from a grid file
	 * @param filename filename of the .grid file
	 * @param mesh reference to a Geometry::Mesh3D
	 * @param properties reference to a Geometry::PropertiesMap
	 */
	ImporterTPFA(const std::string filename, Geometry::Mesh3D & mesh, Geometry::PropertiesMap & properties):
		Importer(filename, mesh, properties) {}

	//! Import from a standard grid file
	/*!
	 * Read points, polygons, polyhedra, zone properties
	 * @param fracturesOn if true, imports the fractures, else the fractures are disabled
	 */
	virtual void import(bool fracturesOn = true);

	//! Add the fractures network
	/*!
	 * Add the fracture network
	 * @pre import the file
	 */
	virtual void addFractures();
	
	//! Destructor
	virtual ~ImporterTPFA() = default;

private:

	//! No default constructor
	ImporterTPFA();

	//! No copy-constructor
	ImporterTPFA(const ImporterTPFA &);

	//! No assignment operator
	ImporterTPFA & operator=(const ImporterTPFA &);
};

//! Class used to read files optimized for the solver (.fvg).
/*!
 * @class ImporterForSolver
 * This class allows to read files optimized for the solver (.fvg).
 * This format contains also the BC id for each polygon and the fracture network.
 */
class ImporterForSolver : public Importer
{
public:

	//! Constructor
	/*!
	 * Constructor from a file
	 * @param filename filename of the file
	 * @param mesh reference to a Geometry::Mesh3D
	 * @param properties reference to a Geometry::PropertiesMap
	 */
	ImporterForSolver(const std::string filename, Geometry::Mesh3D & mesh, Geometry::PropertiesMap & properties):
		Importer(filename, mesh, properties) {}

	//! Import from a file with boundary conditions and fracture network
	/*!
	 * Read points, polygons, polyhedra, zone properties, BC ids, fracture network
	 * @param fracturesOn if true, imports the fractures, else the fractures are disabled
	 */
	virtual void import(bool fracturesOn = true);
	
	//! Destructor
	virtual ~ImporterForSolver() = default;

private:
	//! No default constructor
	ImporterForSolver();

	//! No copy-constructor
	ImporterForSolver(const ImporterForSolver &);

	//! No assignment operator
	ImporterForSolver & operator=(const ImporterForSolver &);
	
	//! Hidden
	void addFractures(){};
	//! Hidden
	using Importer::addBCAndFractures;

};


#endif /* IMPORT_HPP_ */
