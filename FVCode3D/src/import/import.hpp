/*!
 *	@file import.hpp
 *	@brief Classes for loading files.
 */

#ifndef IMPORT_HPP_
#define IMPORT_HPP_

#include "core/TypeDefinition.hpp"

class Mesh3D;
class PropertiesMap;

//! Class used to read TPFA format files.
/*!
 * @class ImporterTPFA
 * This is a base abstract class that allows to read TPFA format files.
 */
class ImporterTPFA
{
public:

	//! Constructor
	/*!
	 * Constructor from a grid file
	 * @param filename filename of the .grid file
	 * @param mesh reference to a Geometry::Mesh3D
	 * @param proporties reference to a Geometry::PropertiesMap
	 */
	ImporterTPFA(const std::string filename, Geometry::Mesh3D & mesh, Geometry::PropertiesMap & properties):
		M_filename(filename), M_mesh(mesh), M_properties(properties) {}

	//! Import from a grid file
	virtual void import() = 0;

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
	virtual ~ImporterTPFA(){};

protected:

	//! Filename
	std::string M_filename;
	//! Reference to a Geometry::Mesh3D
	Geometry::Mesh3D & M_mesh;
	//! Reference to a Geometry::PropertiesMap
	Geometry::PropertiesMap & M_properties;

private:

	//! No default constructor
	ImporterTPFA();

	//! No copy-constructor
	ImporterTPFA(const ImporterTPFA &);

	//! No assignment operator
	ImporterTPFA & operator=(const ImporterTPFA &);

};

//! Class used to read a standard TPFA format file.
/*!
 * @class ImporterTPFAStandard
 * This class allows to read a standard TPFA format files.
 */
class ImporterTPFAStandard : public ImporterTPFA
{
public:

	//! Constructor
	/*!
	 * Constructor from a grid file
	 * @param filename filename of the .grid file
	 * @param mesh reference to a Geometry::Mesh3D
	 * @param proporties reference to a Geometry::PropertiesMap
	 */
	ImporterTPFAStandard(const std::string filename, Geometry::Mesh3D & mesh, Geometry::PropertiesMap & properties):
		ImporterTPFA(filename, mesh, properties) {}

	//! Import from a standard grid file
	/*!
	 * Read points, polygons, polyhedra, zone properties
	 */
	virtual void import();

	//! Destructor
	virtual ~ImporterTPFAStandard(){};

private:

	//! No default constructor
	ImporterTPFAStandard();

	//! No copy-constructor
	ImporterTPFAStandard(const ImporterTPFAStandard &);

	//! No assignment operator
	ImporterTPFAStandard & operator=(const ImporterTPFAStandard &);
};

//! Class used to read an improved version of a TPFA format file.
/*!
 * @class ImporterTPFAWithBC
 * This class allows to read an improved version TPFA format files.
 * In fact, this format contains also the BC id for each polygon and the fracture network.
 */
class ImporterTPFAWithBC : public ImporterTPFA
{
public:

	//! Constructor
	/*!
	 * Constructor from a grid file
	 * @param filename filename of the .grid file
	 * @param mesh reference to a Geometry::Mesh3D
	 * @param proporties reference to a Geometry::PropertiesMap
	 */
	ImporterTPFAWithBC(const std::string filename, Geometry::Mesh3D & mesh, Geometry::PropertiesMap & properties):
		ImporterTPFA(filename, mesh, properties) {}

	//! Import from a grid file with boundary conditions and fracture network
	/*!
	 * Read points, polygons, polyhedra, zone properties, BC ids, fracture network
	 */
	virtual void import();

	//! Destructor
	virtual ~ImporterTPFAWithBC(){};

private:
	//! No default constructor
	ImporterTPFAWithBC();

	//! No copy-constructor
	ImporterTPFAWithBC(const ImporterTPFAWithBC &);

	//! No assignment operator
	ImporterTPFAWithBC & operator=(const ImporterTPFAWithBC &);

};


#endif /* IMPORT_HPP_ */
