/*!
 *	@file export.hpp
 *	@brief Classes for saving files.
 */

#ifndef EXPORT_HPP_
#define EXPORT_HPP_

#include "core/TypeDefinition.hpp"

class Mesh3D;
class Rigid_Mesh;

//! Flag for exporter
/*!
 * @enum PropertiesToExport
 * This enum permits to select which properties to export.
 */
enum PropertiesToExport
{
	ZoneCode 		= 0x01,
	Aperture		= 0x02,
	Permeability	= 0x04,
	Porosity		= 0x08,
	BorderID		= 0x10,
	Pressure		= 0x20,
	Velocity		= 0x40,
	Other			= 0x80,
};

//! Class to export mesh, fractures, solution and properties.
/*!
 * @class Exporter
 * This is a base abstact class that allows to export mesh, fractures, solution, and properties.
 * Each derived class permits to export in a specific file format.
 */
class Exporter
{
public:

	typedef Geometry::Mesh3D Mesh3D;
	typedef Geometry::Mesh3D::Cell3D Cell3D;
	typedef Geometry::Mesh3D::Facet3D Facet3D;

	typedef Geometry::Rigid_Mesh Rigid_Mesh;
	typedef Geometry::Rigid_Mesh::Cell Cell;
	typedef Geometry::Rigid_Mesh::Facet Facet;
	typedef Geometry::Rigid_Mesh::Fracture_Facet Fracture_Facet;
	typedef Geometry::Rigid_Mesh::Fracture_Juncture Fracture_Juncture;

	//! Contructor
	Exporter(){}

	//! Export the mesh (only cells)
	/*!
	 * @param mesh reference of a Geometry::Mesh3D
	 * @param filename name of the file
	 */
	virtual void exportMesh(const Mesh3D & mesh, const std::string filename) = 0;

	//! Export the fracture facets
	/*!
	 * @param mesh reference of a Geometry::Mesh3D
	 * @param filename name of the file
	 */
	virtual void exportFractures(const Mesh3D & mesh, const std::string filename) = 0;

	//! Export the mesh, cells and fracture facets, in a single file
	/*!
	 * @param mesh reference of a Geometry::Mesh3D
	 * @param filename name of the file
	 */
	virtual void exportMeshWithFractures(const Mesh3D & mesh, const std::string filename) = 0;

	//! Export the fracture junctures
	/*!
	 * @param mesh reference of a Geometry::Rigid_Mesh
	 * @param filename name of the file
	 */
	virtual void exportFractureJunctures(const Rigid_Mesh & mesh, const std::string filename) = 0;

	//! Export the solution on cells and fracture facets in a single file
	/*!
	 * @param mesh reference of a Geometry::Rigid_Mesh
	 * @param filename name of the file
	 * @param sol Eigen vector that contain the solution (cells + fracture facets)
	 */
	virtual void exportSolution(const Rigid_Mesh & mesh, const std::string filename, const Eigen::VectorXd & sol) = 0;

	//! Export the solution on fracture facets
	/*!
	 * @param mesh reference of a Geometry::Rigid_Mesh
	 * @param filename name of the file
	 * @param sol Eigen vector that contain the solution (cells + fracture facets)
	 */
	virtual void exportSolutionOnFractures(const Rigid_Mesh & mesh, const std::string filename, const Eigen::VectorXd & sol) = 0;

	//! Export the a specific property on cells and fracture facets
	/*!
	 * @param mesh reference of a Geometry::Rigid_Mesh
	 * @param filename name of the file
	 * @param propertiesType flag used to select which properties to export
	 * @param property pointer to a generic property
	 */
	virtual void exportWithProperties(const Rigid_Mesh & mesh, const std::string filename, const Flag8bit propertiesType, const std::vector<Real> * property = NULL ) = 0;

	//! Export the a specific property on cells and fracture facets
	/*!
	 * @param mesh reference of a Geometry::Mesh3D
	 * @param filename name of the file
	 * @param propertiesType flag used to select which properties to export
	 * @param property pointer to a generic property
	 */
	virtual void exportWithProperties(const Mesh3D & mesh, const std::string filename, const Flag8bit propertiesType, const std::vector<Real> * property = NULL ) = 0;

	//! Destructor
	virtual ~Exporter(){};

};

//! Class to export mesh, fractures, solution and properties in VTU format
/*!
 * @class ExporterVTU
 * This class allows to export mesh, fractures, solution, and properties in VTU format.
 */
class ExporterVTU : public Exporter
{
public:

	//! Constructor
	ExporterVTU(){}

	//! Export the mesh (only cells)
	/*!
	 * @param mesh reference of a Geometry::Mesh3D
	 * @param filename name of the file
	 */
	virtual void exportMesh(const Mesh3D & mesh, const std::string filename);

	//! Export the fracture facets
	/*!
	 * @param mesh reference of a Geometry::Mesh3D
	 * @param filename name of the file
	 */
	virtual void exportFractures(const Mesh3D & mesh, const std::string filename);

	//! Export the mesh, cells and fracture facets, in a single file
	/*!
	 * @param mesh reference of a Geometry::Mesh3D
	 * @param filename name of the file
	 */
	virtual void exportMeshWithFractures(const Mesh3D & mesh, const std::string filename);

	//! Export the fracture junctures
	/*!
	 * @param mesh reference of a Geometry::Rigid_Mesh
	 * @param filename name of the file
	 */
	virtual void exportFractureJunctures(const Rigid_Mesh & mesh, const std::string filename);

	//! Export the solution on cells and fracture facets in a single file
	/*!
	 * @param mesh reference of a Geometry::Rigid_Mesh
	 * @param filename name of the file
	 * @param sol Eigen vector that contain the solution (cells + fracture facets)
	 */
	virtual void exportSolution(const Rigid_Mesh & mesh, const std::string filename, const Eigen::VectorXd & sol);

	//! Export the solution on fracture facets
	/*!
	 * @param mesh reference of a Geometry::Rigid_Mesh
	 * @param filename name of the file
	 * @param sol Eigen vector that contain the solution (cells + fracture facets)
	 */
	virtual void exportSolutionOnFractures(const Rigid_Mesh & mesh, const std::string filename, const Eigen::VectorXd & sol);

	//! Export the solution on fracture facets
	/*!
	 * @param mesh reference of a Geometry::Rigid_Mesh
	 * @param filename name of the file
	 * @param sol Eigen vector that contain the solution (cells + fracture facets)
	 */
	virtual void exportWithProperties(const Rigid_Mesh & mesh, const std::string filename, const Flag8bit propertiesType, const std::vector<Real> * property = NULL ) ;

	//! Export the a specific property on cells and fracture facets
	/*!
	 * @param mesh reference of a Geometry::Mesh3D
	 * @param filename name of the file
	 * @param propertiesType flag used to select which properties to export
	 * @param property pointer to a generic property
	 */
	virtual void exportWithProperties(const Mesh3D & mesh, const std::string filename, const Flag8bit propertiesType, const std::vector<Real> * property = NULL);

	//! Destrcutor
	virtual ~ExporterVTU() {};
};

//! Class to export mesh, fractures, solution and properties in VTK format
/*!
 * @class ExporterVTK
 * This class allows to export mesh, fractures, solution, and properties in VTK format.
 * @warning Not implemented.
 */
class ExporterVTK : public Exporter
{
public:

	ExporterVTK(){}

	virtual void exportMesh(const Mesh3D & mesh, const std::string filename);

	virtual void exportFractures(const Mesh3D & mesh, const std::string filename);

	virtual void exportMeshWithFractures(const Mesh3D & mesh, const std::string filename);

	virtual void exportFractureJunctures(const Rigid_Mesh & mesh, const std::string filename);

	virtual void exportSolution(const Rigid_Mesh & mesh, const std::string filename, const Eigen::VectorXd & sol);

	virtual void exportSolutionOnFractures(const Rigid_Mesh & mesh, const std::string filename, const Eigen::VectorXd & sol);

	virtual void exportWithProperties(const Rigid_Mesh & mesh, const std::string filename, const Flag8bit propertiesType, const std::vector<Real> * property = NULL );

	virtual void exportWithProperties(const Mesh3D & mesh, const std::string filename, const Flag8bit propertiesType, const std::vector<Real> * property = NULL);

	virtual ~ExporterVTK() {};
};

#endif /* EXPORT_HPP_ */
