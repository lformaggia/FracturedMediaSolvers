/*!
 *  @file Import.hpp
 *  @brief Classes for loading files.
 */

#ifndef IMPORT_HPP_
#define IMPORT_HPP_

#include <FVCode3D/core/TypeDefinition.hpp>

namespace FVCode3D
{

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
     * @param mesh reference to a Mesh3D
     * @param properties reference to a PropertiesMap
     */
    Importer(const std::string filename, Mesh3D & mesh, PropertiesMap & properties):
        M_filename(filename), M_mesh(mesh), M_properties(properties) {}

    //! Import from a grid file
    /*!
     * @param fracturesOn if true, imports the fractures, else the fractures are disabled
     */
    virtual void import(bool fracturesOn) throw() = 0;

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
    virtual void addFractures();

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
    const Mesh3D & getMesh() const
        { return M_mesh; }

    //! Get properties
    /*!
     * @return reference to the properties
     */
    const PropertiesMap & getProperties() const
        { return M_properties; }

    //! Destructor
    virtual ~Importer() = default;

protected:

    //! Filename
    std::string M_filename;
    //! Reference to a Mesh3D
    Mesh3D & M_mesh;
    //! Reference to a PropertiesMap
    PropertiesMap & M_properties;

private:

    //! No default constructor
    Importer() = delete;

    //! No copy-constructor
    Importer(const Importer &) = delete;

    //! No assignment operator
    Importer & operator=(const Importer &) = delete;

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
     * @param mesh reference to a Mesh3D
     * @param properties reference to a PropertiesMap
     */
    ImporterMedit(const std::string filename, Mesh3D & mesh, PropertiesMap & properties):
        Importer(filename, mesh, properties) {}

    //! Import from a .mesh file
    /*!
     * Read points, polygons, polyhedra
     * @param fracturesOn if true, imports the fractures, else the fractures are disabled
     */
    virtual void import(bool fracturesOn = true) throw();

    //! Destructor
    virtual ~ImporterMedit() = default;

private:

    //! No default constructor
    ImporterMedit() = delete;

    //! No copy-constructor
    ImporterMedit(const ImporterMedit &) = delete;

    //! No assignment operator
    ImporterMedit & operator=(const ImporterMedit &) = delete;
};

//! Class used to read TetGen format files (.node, .face, .ele).
/*!
 * @class ImporterTetGen
 * This class allows to read TetGen format files (.node, .face, .ele).
 */
class ImporterTetGen : public ImporterMedit
{
public:

    //! Constructor
    /*!
     * Constructor from TetGen files
     * @param filename prefix filename of the .node, .face, .ele files
     * @param mesh reference to a Mesh3D
     * @param properties reference to a PropertiesMap
     */
    ImporterTetGen(const std::string filename, Mesh3D & mesh, PropertiesMap & properties):
        ImporterMedit(filename, mesh, properties) {}

    //! Import from .node, .face, .ele files
    /*!
     * Read points, faces, tetrahedra
     * @param fracturesOn if true, imports the fractures, else the fractures are disabled
     */
    virtual void import(bool fracturesOn = true) throw();

    //! Destructor
    virtual ~ImporterTetGen() = default;

private:

    //! No default constructor
    ImporterTetGen() = delete;

    //! No copy-constructor
    ImporterTetGen(const ImporterTetGen &) = delete;

    //! No assignment operator
    ImporterTetGen & operator=(const ImporterTetGen &) = delete;
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
     * @param mesh reference to a Mesh3D
     * @param properties reference to a PropertiesMap
     */
    ImporterTPFA(const std::string filename, Mesh3D & mesh, PropertiesMap & properties):
        Importer(filename, mesh, properties) {}

    //! Import from a standard grid file
    /*!
     * Read points, polygons, polyhedra, zone properties
     * @param fracturesOn if true, imports the fractures, else the fractures are disabled
     */
    virtual void import(bool fracturesOn = true) throw();

    //! Destructor
    virtual ~ImporterTPFA() = default;

private:

    //! No default constructor
    ImporterTPFA() = delete;

    //! No copy-constructor
    ImporterTPFA(const ImporterTPFA &) = delete;

    //! No assignment operator
    ImporterTPFA & operator=(const ImporterTPFA &) = delete;
};

//! Class used to read the OpenFOAM format file.
/*!
 * @class ImporterOpenFOAM
 * This class allows to read the OpenFOAM format file.
 */
class ImporterOpenFOAM : public Importer
{
public:

    //! Constructor
    /*!
     * Constructor from the OpenFOAM files
     * @param filename filename of the OpenFOAM files
     * @param mesh reference to a Mesh3D
     * @param properties reference to a PropertiesMap
     */
    ImporterOpenFOAM(const std::string filename, Mesh3D & mesh, PropertiesMap & properties):
        Importer(filename, mesh, properties) {}

    //! Import from the OpenFOAM files
    /*!
     * Read points, faces, owner, neighbour, boundary
     * @param fracturesOn if true, imports the fractures, else the fractures are disabled
     */
    virtual void import(bool fracturesOn = true) throw();

    //! Destructor
    virtual ~ImporterOpenFOAM() = default;

private:

    //! No default constructor
    ImporterOpenFOAM() = delete;

    //! No copy-constructor
    ImporterOpenFOAM(const ImporterOpenFOAM &) = delete;

    //! No assignment operator
    ImporterOpenFOAM & operator=(const ImporterOpenFOAM &) = delete;
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
     * @param mesh reference to a Mesh3D
     * @param properties reference to a PropertiesMap
     */
    ImporterForSolver(const std::string filename, Mesh3D & mesh, PropertiesMap & properties):
        Importer(filename, mesh, properties) {}

    //! Import from a file with boundary conditions and fracture network
    /*!
     * Read points, polygons, polyhedra, zone properties, BC ids, fracture network
     * @param fracturesOn if true, imports the fractures, else the fractures are disabled
     */
    virtual void import(bool fracturesOn = true) throw();

    //! Destructor
    virtual ~ImporterForSolver() = default;

private:
    //! No default constructor
    ImporterForSolver() = delete;

    //! No copy-constructor
    ImporterForSolver(const ImporterForSolver &) = delete;

    //! No assignment operator
    ImporterForSolver & operator=(const ImporterForSolver &) = delete;

    //! Hidden
    void addFractures(){};
    //! Hidden
    using Importer::addBCAndFractures;
};

} // namespace FVCode3D
#endif /* IMPORT_HPP_ */
