/*!
 *  @file ProxyRigidMesh.hpp
 *  @brief Class for unstructured mesh. Permits to edit the RigidMesh class.
 */

#ifndef PROXY_RIGID_MESH_HPP_
#define PROXY_RIGID_MESH_HPP_

#include <FVCode3D/core/TypeDefinition.hpp>
#include <FVCode3D/mesh/RigidMesh.hpp>

namespace FVCode3D
{

//! Class the allows to edit a Rigid_Mesh::Edge object
/*!
 * @class ProxyEdge
 * This class provides some methods that could be used to edit a Ridig_Mesh::Edge.
 */
class ProxyEdge
{
public:

    //! Typedef for Rigid_Mesh::Edge
    /*!
        @typedef Edge
        This type definition permits to treat Rigid_Mesh::Edge as a Edge.
    */
    typedef Rigid_Mesh::Edge Edge;

    //! No default constructor
    ProxyEdge() = delete;

    //! Constructor
    /*!
     * @param edge reference to a Rigid_Mesh::Border_Edge
     */
    ProxyEdge(Edge & edge): M_edge(edge) {}

    //! Default destructor
    ~ProxyEdge() = default;

    //! Set if the edge is a border edge (for the domain)
    /*!
     * @param isBorder True if the edge is a border edge (i.e. it separates a border facet)
     */
    void setBorderEdge(bool isBorder)
        { M_edge.M_isBorderEdge = isBorder; }

    //! Set if the edge is a tip, i.e. a border edge for the fracture
    /*!
     * @param True if the edge is a tip (i.e. it separates one of the fractures that it separates only once, i.e. it is a tip for at least a fracture)
     */
    void setTip(bool isTip)
        { M_edge.M_isTip = isTip; }

private:

    //! Reference to a Rigid_Mesh::Edge
    Edge & M_edge;
};

//! Class the allows to edit a Rigid_Mesh::Border_Edge object
/*!
 * @class ProxyBorderEdge
 * This class provides some methods that could be used to edit a Ridig_Mesh::Border_Edge.
 */
class ProxyBorderEdge
{
public:

    //! Typedef for Rigid_Mesh::Border_Edge
    /*!
        @typedef Border_Edge
        This type definition permits to treat Rigid_Mesh::Border_Edge as a Border_Edge.
    */
    typedef Rigid_Mesh::Border_Edge Border_Edge;

    //! No default constructor
    ProxyBorderEdge() = delete;

    //! Constructor
    /*!
     * @param edge reference to a Rigid_Mesh::Border_Edge
     */
    ProxyBorderEdge(Border_Edge & edge): M_edge(edge) {}

    //! Default destructor
    ~ProxyBorderEdge() = default;

    //! Get border ids
    /*!
     * @return the ids of the represented border facets (inherited from the facets)
     */
    std::set<UInt> & getBorderIds()
        {return M_edge.M_borderIds;}

private:

    //! Reference to a Rigid_Mesh::Border_Edge
    Border_Edge & M_edge;
};

//! Class the allows to edit a Rigid_Mesh object
/*!
 * @class ProxyRigidMesh
 * This class provides some methods that could be used to edit a Ridig_Mesh.
 */
class ProxyRigidMesh
{
public:

    //! Typedef for Rigid_Mesh::Edge
    /*!
        @typedef Edge
        This type definition permits to treat Rigid_Mesh::Edge as a Edge.
    */
    typedef Rigid_Mesh::Edge Edge;

    //! Typedef for Rigid_Mesh::Border_Tip_Edge
    /*!
        @typedef Border_Tip_Edge
        This type definition permits to treat Rigid_Mesh::Border_Tip_Edge as a Border_Tip_Edge.
    */
    typedef Rigid_Mesh::Border_Tip_Edge Border_Tip_Edge;

    //! @name Constructor & Destructor
    //@{

    //! No default constructor
    ProxyRigidMesh() = delete;

    //! Constructor
    /*!
     * @param mesh reference to a Rigid_Mesh
     */
    ProxyRigidMesh(Rigid_Mesh & mesh):M_mesh(mesh){}

    //! Default destructor
    ~ProxyRigidMesh() = default;

    //@}

    //! @name Methods
    //@{

    //! Get edges vector
    /*!
     * @return a reference to the vector that contains the edges of the mesh
     */
    std::vector<Edge> & getEdgesVector ()
        { return M_mesh.M_edges; }

    //! Get border tip edges ids vector
    /*!
     * @return a reference to the vector of the M_borderTipEdges
     */
    std::vector<Border_Tip_Edge> & getBorderTipEdgesIdsVector ()
        { return M_mesh.M_borderTipEdges; }

    //! Builds the vector of RigidMesh::Edge.
    /*!
     * @pre This method is used to build the Edge objects inside the RigidMesh
     * if they are not built yet.
     */
    void edgesBuilder()
        { M_mesh.edgesBuilder(); }

    //! Builds the vector of RigidMesh::Edge_ID.
    /*!
     * @pre This method is used to build the Edge_ID objects inside the RigidMesh
     * if they are not built yet.
     */
    void edgesIdBuilder()
        { M_mesh.edgesIdBuilder(); };

    //@}

private:

    //! Reference to a Rigid_Mesh
    Rigid_Mesh & M_mesh;
};

} // namespace FVCode3D

#endif // PROXY_RIGID_MESH_HPP_
