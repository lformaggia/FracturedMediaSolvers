/**
 *  \file header.hpp
 *  A modified version of the code produced by Pigoli Davide and
 *  Prada Daniele
 *  \author Luca Formaggia 2013
 */

#ifndef HEADER_HPP_
#define HEADER_HPP_

#include <iostream>
#include <vector>
#include "domain.hpp"
#include "treenode.hpp"
#include "exception_handling.hpp"

namespace ADT
{
/** \class Header
 *  \brief It contains general information about the tree.
 *  \param  T The template parameter expresses the criterion for building the tree depending on whether you have to store points or boxes. \n
 *              So, you have:
 *              <ul>
 *              <li> T = Point<NDIMP>, in the first case
 *              <li> T = Box<NDIMP>, in the second case
 *              </ul>
 *              where NDIMP is the number of physical space dimensions.
 *
 *  Default copy constructor and default destructor work fine.
 */
template<class T>
class Header
{
protected:
    /// Tree memory locations.
    int tree_loc;
    /// Tree levels.
    int tree_lev;
    /// Number of physical space dimensions (typically 2 or 3).
    int ndimp;
    /// Number of pieces of information carried by the tree. The size of each node is nkey+2 words.
    int nkey;
    /// Number of dimensions used for the search (either ndimp or 2*ndimp).
    int ndimt;
    /// Number of logical locations currently used in the tree. Initialized to 0.
    int nele;
    /** @name Tree indices
     *  The use of iava and iend avoids the necessity of initializing the stack of available
     *  locations. In fact, the stack contains only locations that have been previously deleted from the tree.
     */
    //@{
    /** Next available location in the stack. Initialized to 1.
     *
     *  The stack of available nodes is a linked list with all nodes that have been previously deleted. If iava=0 the stack is empty.
     */
    int iava;
    /** Next available location in the yet not allocated part of the tree ("tree free store"). Initialized to 1.
     *
     *  The "tree free store" is the yet unassigned portion of the vector storing the tree.
     *  It effectively acts as a free storage area. If iend=1, the free storage is at his maximum size. iend never decreases.
     *  It may increase during insertions and will always remain unaltered after deletions.
     *  It isn't equal to nele.
     */
    int iend;
    //@}
    /// Tree's domain.
    Domain<T> tree_domain;
    /** A protected constructor.
     *
     *  \param[in] ntree Tree dimension needed.
     *  \param[in] nk Number of pieces of information carried by the tree.
     *  \param[in] d Tree's domain.
     *
     *  Public function createtreeheader must be used  to create an
     *  object of Header class.  This avoids the creation of a tree
     *  with more memory locations than a fixed limit.
     */
    Header (int const& ntree, int const& nk, Domain<T> const& d);
    /// Tries to set the number of tree memory locations (throws a
    /// LocLengthError exception if nt is out of range).
    void stml (int const& nt);
public:
    /** Default constructor.
     *
     *  It's fundamental in creating an ADTree object from a MeshFile::ff2dmesh or a MeshFile::ff3dmesh object.
     */
    Header() {}
    /// Gets the number of tree memory locations.
    inline int gettreeloc() const
    {
        return tree_loc;
    }
    /// Sets the number of tree memory locations (handles a LocLengthError exception).
    void settreeloc (int const& nt);
    /// Gets the number of tree levels.
    inline int gettreelev() const
    {
        return tree_lev;
    }
    /// Sets the number of tree levels.
    inline void settreelev (int const& nl)
    {
        tree_lev = nl;
    }
    /// Gets the number of physical space dimension.
    inline int getndimp() const
    {
        return ndimp;
    }
    /// Gets the number of pieces of information carried by the tree.
    inline int getnkey() const
    {
        return nkey;
    }
    /// Gets the number of dimensions used for the search.
    inline int getndimt() const
    {
        return ndimt;
    }
    /// Gets the number of logical locations currently used in the tree.
    inline int getnele() const
    {
        return nele;
    }
    /// Sets the number of logical locations currently used in the tree.
    inline void setnele (int const& ne)
    {
        nele = ne;
    }
    /// Gets the next available location in the stack.
    inline int getiava() const
    {
        return iava;
    }
    /// Sets the next available location in the stack.
    inline void setiava (int const& ia)
    {
        iava = ia;
    }
    /// Gets the next available location in the tree free store.
    inline int getiend() const
    {
        return iend;
    }
    /// Sets the next available location in the tree free store.
    inline void setiend (int const& ie)
    {
        iend = ie;
    }
    /// Gets the i-th coordinate of the origin of the tree's bounding box.
  std::vector<double> const& domainorig () const
    {
        return tree_domain.orig ();
    }
    /// Gets the i-th scaling factor of the tree's bounding box.
    std::vector<double> const& domainscal () const
    {
        return tree_domain.scal ();
    }
    /** Output operator.
     *
     *  It prints out:
     *  - number of tree memory locations;
     *  - number of tree levels;
     *  - number of physical space dimensions;
     *  - number of pieces of information carried by the tree;
     *  - number of dimensions used for the search;
     *  - number of logical locations currently used in the tree;
     *  - tree domain.
     */
    template<class S>
    friend std::ostream& operator<< (std::ostream&, Header<S> const&);
    /** \fn Header<S> createtreeheader(int const & nt, int const & nk, Domain<S> const & d)
     *  \brief Creates a tree header handling the exception condition in which tree maximum dimension is exceeded.
     *  \param[in] nt Tree dimension needed.
     *  \param[in] nk Number of pieces of information carried by the tree.
     *  \param[in] d Tree's domain.
     *
     *  The template parameter expresses the criterion for building the tree.
     */
    template<class S>
    friend Header<S> createtreeheader (int const& nt, int const& nk, Domain<S> const& d);
};

template<class T>
Header<T>::Header (int const& ntree, int const& nk, Domain<T> const& d) :
    tree_loc (ntree), tree_lev (0), ndimp (3), nkey (nk), ndimt (6), nele (0), iava (1), iend (1), tree_domain (d)
{

    std::vector<TreeNode> foo;
    if (foo.max_size() < unsigned (tree_loc + 1) )
        /* If there is no enough space to store the requested nodes and
         * the tree head, a LocLengthError exception is thrown.
         */
    {
        throw (LocLengthError<T> (foo.max_size(), tree_loc) );
    }
}

template<class T>
void Header<T>::stml (int const& nt)
{
    std::vector<TreeNode> foo;
    if (foo.max_size() < unsigned (nt + 1) )
        /* If there is no enough space to store the requested nodes and
         * the tree head, a LocLengthError exception is thrown.
         */
    {
        throw (LocLengthError<T> (foo.max_size(), nt) );
    }
    tree_loc = nt;
}

template<class T>
std::ostream& operator<< (std::ostream& ostr, Header<T> const& head)
{
    ostr << std::endl << std::endl;
    ostr << "General informations about the tree" << std::endl;
    ostr << "----------------------------------" << std::endl << std::endl;
    ostr << "Tree memory locations: " << head.tree_loc << std::endl;
    ostr << "Number of tree levels: " << head.tree_lev << std::endl;
    ostr << "Number of physical space dimension: " << head.ndimp << std::endl;
    ostr << "Number of pieces of information carried by the tree: " << head.nkey << std::endl;
    ostr << "Number of dimensions used for the search: " << head.ndimt << std::endl;
    ostr << "Number of logical locations currently used in the tree: " << head.nele;
    ostr << head.tree_domain << std::endl;

    return ostr;
}

template<class T>
Header<T> createtreeheader (int const& nt, int const& nk, Domain<T> const& d)
{
    try
    {
        Header<T> hd (nt, nk, d);
        return hd;
    }
    catch (LocLengthError<T> lo)
    {
        std::cout << std::endl << std::endl;
        std::cout << "warning!	createtreeheader : max dimension exceeded" << std::endl;
        std::cout << "the limit is " << lo.getmaxtreeloc() - 1
                  << " while needed at least " << lo.gettreeloc() << std::endl;
        std::exit (EXIT_FAILURE);
    }
}

template<class T>
void Header<T>::settreeloc (int const& nt)
{
    try
    {
        stml (nt);
    }
    catch (LocLengthError<T> lo)
    {
        std::cout << std::endl << std::endl;
        std::cout << "warning!	settreeloc : max dimension exceeded" << std::endl;
        std::cout << "the limit is " << lo.getmaxtreeloc() - 1
                  << " while requested " << lo.gettreeloc() << std::endl;
        std::cout << "increasing tree memory locations up to the limit" << std::endl;
        stml (lo.getmaxtreeloc() - 1);
    }
}

}// end namspace ADT

#endif /* HEADER_HPP_ */
