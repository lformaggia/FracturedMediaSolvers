/**
 *  \file adtree.hpp
 *  A modified version of the code produced by Pigoli Davide and Prada Daniele
 *  \author Luca Formaggia 2013
 */

#ifndef ADTREE_HPP_
#define ADTREE_HPP_

#include <iostream>
#include <vector>
#include <cmath>
#include "domain.hpp"
#include "header.hpp"
#include "treenode.hpp"

// Forward declaration to sort out some
// circular dependencies. Indeed we should take out
// explicit dependencies on CPGrid here and do a cleaner programming.
// ... but I am lazy.

namespace Geometry
{
  class CPgrid;
}

namespace ADT
{
  /** \class ADTree
   *  \brief Alternating binary range searching tree.
   */
  class ADTree
  {
  protected:
    //! a helper structure
    struct Helper
    {
      std::vector<double> xl;
      int ipoi;
      int lev;
    };
    /** The header.
     *
     *  It contains general information about the tree.
     */
    Header<Geometry::CPgrid> _header;
    /// Vector of tree nodes.
    std::vector<TreeNode> _data;
    /** \brief Adds a node to the tree.
        Remember that the node must be properly formed before
        adding it to
        It throws: <ul>
        <li> a
        TreeDomainError exception if the point is out of domain;
        <li> a TreeAlloc exception if there is no more space in the
        tree to add the node; <li> a LevRuntimeError if you exceed
        the limit set for the tree levels due to the inclusion of
        the node.  </ul>
    */
    int adtrb (TreeNode const& node);
    /** Searches dimension associated to a given level.
     *
     *    \param[in] lev The given level.
     *    \param[in] dim The number of dimensions used for the search.
     */
    inline int searchdim (int const& lev, int const& dim) const
    {
      return (lev % dim);
    }
    /** Finds delta associated to division at a given level.
     *
     *    \param[in] lev The given level.
     *    \param[in] dim The number of dimensions used for the search.
     */
    /*    inline double delta (int const& lev, int const& dim) const
      {
          return std::pow (0.5, int (lev / dim) + 1);
    }*/
    inline double delta (int const& lev, int const& dim) const
    {
      double res (0.5);
      for (int i = 0; i < int (lev / dim); ++i) res *= 0.5;
      return res;
    }
  public:
    /*!  Constructor froma corner point grid.
     */
    ADTree (Geometry::CPgrid const& grid);
    /// Returns a reference to the tree header.
    inline Header<Geometry::CPgrid> const& gettreeheader() const
    {
      return _header;
    }
    /** Adds a node to the tree.
     *    It calls the handlers of the exceptions that can be thrown by adtrb().
     *
     *    \param[in] coords Coordinates of the point.
     *    \param[in] keys Additional informations associated to the point.
     *
     *    The location of the current node in the tree is returned.
     */
    int addtreenode (TreeNode const&);
    /** Gets out the info related to a node at a given location
     *
     *    \param[in] loc Location of the searched node.
     *    \param[out] info The info stored in the node
     *    \return A constant reference to the node
     */
    TreeNode const&   getNode (int const& loc) const
    {
      return _data[loc];
    }
    /** Finds all points or (bounding) boxes that intersect a given
     * box.
     *
     *    \param[in]  region. Bounding Box
     *    \param[out] found Indices of the found
     *    elements.
     *    \return true if it has completed successfully, false otherwise.
     */
    bool search (BBox<3> region, std::vector<int>& found) const;
    /// Deletes a specified location in the tree.
    void deltreenode (int const& index);
    /// Outputs informations contained in the tree header.
    friend std::ostream& operator<< (std::ostream& ostr, ADTree const& myadt);
  };
}// end namespace
#endif /* ADTREE_HPP_ */
