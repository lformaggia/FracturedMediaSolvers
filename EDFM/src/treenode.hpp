/**
 *  \file treenode.hpp
 *  Modified version of the code produced by Davide Pigoli and
 *  Daniele Prada.
 *  \author Luca Formaggia
 */

#ifndef TREENODE_HPP_
#define TREENODE_HPP_

#include <vector>
#include "bbox.hpp"
namespace ADT
{
/** \class TreeNode
 *  \brief Class defining a tree node.
 */
class TreeNode
{
public:
    typedef BBox<3> Shape;
protected:
    /** Position of the father node.
     *
     *    It's used by the algorithm for deleting a tree node.
     */
    int _father;
    /// Positions of left and right children.
    int _children[2];
    /// Geometric object to be stored in the node.
    Shape _shape;
    /** Additional informations to be stored in the node.
     *
     *    Maybe a more sophisticated object than a std::vector<int> can be more useful.
     */
    std::vector<int> _keys;
public:
    /**   Default constructor.
     *
     *    It's fundamental in creating a vector of TreeNode objects.
     */
    TreeNode() : _father (0), _shape()
    {
        _children[0] = 0;
        _children[1] = 0;
    }
    /**   Another constructor.
     *
     *    It uses the automatic copy constructor for _shape member.
     */
    TreeNode (Shape const& shape,
              std::vector<int> const& keys = std::vector<int>() ) :
        _father (0),
        _shape (shape),
        _keys (keys)
    {
        _children[0] = 0;
        _children[1] = 0;
    }
    /// Sets the father.
    inline void setfather (int const& ifth)
    {
        _father = ifth;
    }
    /**   \brief Sets a child.
     *
     *    \param[in] flag Index of the child to be set. \n
     *                    If:
     *                    <ul>
     *                    <li> flag = 0, set the left child
     *                    <li> flag = 1, set the right child
     *                    </ul>
     *    \param[in] child The child node.
     */
    inline void setchild (short int const& flag, int const& child)
    {
        _children[flag] = child;
    }
    /// Returns the father.
    inline int getfather() const
    {
        return _father;
    }
    /**   \brief Returns the id of a  child.
     *
     *    \param[in] flag Index of the child to be returned. \n
     *                    If:
     *                    <ul>
     *                    <li> flag = 0, return the left child
     *                    <li> flag = 1, return the right child
     *                    </ul>
     */
    inline int getchild (short int const& flag) const
    {
        return _children[flag];
    }
    //! Returns the bounding box (read only)
    Shape const& boundingBox() const
    {
        return _shape;
    };
    //! Set boundaing box
    void setBoundingBox (Shape const& s)
    {
        _shape = s;
    }
    /// Sets the additional informations to be stored in the node.
    inline void setkeys (std::vector<int> const& k)
    {
        _keys = k;
    }
    /// Gets the additional informations stored in the node.
    inline std::vector<int> getkeys() const
    {
        std::vector<int> info (_keys);
        return info;
    }
    /// Returns a reference to _shape member.
    inline Shape const& getshape() const
    {
        return _shape;
    }
    /// Returns a reference to _shape member.
    void setshape (Shape const& s)
    {
        _shape = s;
    }
};
}// end namespace ADT

#endif /* TREENODE_HPP_ */
