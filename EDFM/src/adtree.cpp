#include <stack>
#include <exception>
#include "adtree.hpp"
namespace ADT
{
  ADTree::ADTree (Geometry::CPgrid const& grid)
{
  Domain<Geometry::CPgrid> mydom (grid);
    int numCells = grid.Nx() * grid.Ny() * grid.Nz();
    Header<Geometry::CPgrid> myhead =
      createtreeheader<Geometry::CPgrid> (numCells, 3, mydom);
    _header = myhead;
    _data.reserve (_header.gettreeloc() + 1);
    std::vector<int> cellKeys (3, 0);
    _data.push_back (TreeNode (BBox<3>(), cellKeys) );
    std::vector<double> const& origin = _header.domainorig();
    std::vector<double> const& scale = _header.domainscal();
    double const* origin_ptr = & (origin[0]);
    double const* scale_ptr  = & (scale[0]);
    int numVertices (0);
    double vmin[3];
    double vmax[3];
    for (unsigned int i = 1; i <= grid.Nx(); ++i)
    {
        for (unsigned int j = 1; j <= grid.Ny(); ++j)
        {
            for (unsigned int k= 1; k <= grid.Nz(); ++k)
            {
                cellKeys[0] = static_cast<int> (i);
                cellKeys[1] = static_cast<int> (j);
                cellKeys[2] = static_cast<int> (k);
                Geometry::CPcell cell =
                    grid.cell (i, j, k);
                const std::vector<Geometry::Point3D>&
                cellVertices = cell.getVerticesVector();
                numVertices = cellVertices.size();
		std::vector<double> xcoor[3];
                xcoor[0].reserve (numVertices);
                xcoor[1].reserve (numVertices);
                xcoor[2].reserve (numVertices);
                for (std::vector<Geometry::Point3D>::const_iterator it = cellVertices.begin(); it != cellVertices.end(); ++it)
                {
                    xcoor[0].push_back (it->x);
                    xcoor[1].push_back (it->y);
                    xcoor[2].push_back (it->z);
                }
                for (unsigned int l = 0; l < 3; ++l)
                {
                    vmin[l] = *std::min_element (xcoor[l].begin(), xcoor[l].end() );
                    vmax[l] = *std::max_element (xcoor[l].begin(), xcoor[l].end() );
                }
                BBox<3> box (vmin, vmax);
                box.transform (origin_ptr, scale_ptr);
                this->addtreenode (TreeNode (box, cellKeys) );
            }
        }
    }
}


int ADTree::adtrb (TreeNode const& node)
{
    /*
     * We will trasverse the tree in preorder fashion.
     * ipoi and ifth will store the location of the current node
     * and that of the father node.
     *
     * A value 0 for the left_link or the right_link indicates
     * a void sibling address.
     * Therefore, if left_link = right_link = 0 the node is a leaf.
     *
     * In the algorithm, when ipoi = 0 we have
     * reached a void sibling, where we can add a new node.
     */
    int nele = _header.getnele();
    int iava = _header.getiava();
    int iend = _header.getiend();
    int dimt = _header.getndimt();

    // At the start ipoi is set to the address of the root and ifth to 0.
    int ipoi = _data[0].getchild (0);
    int ifth = 0;

    /* We scale the dimension of the "bounding box" of the T object
     * with coordinate values given by coords.
     */
    typedef TreeNode::Shape Shape;
    Shape  x (node.boundingBox() );
    for (int id = 0; id < dimt; ++id)
        if (x[id] > 1.0 || x[id] < 0.0)
        {
            throw std::runtime_error ("Node is not within domain bounds");
        }
    /* x now stores the scaled coordinated of the object's bounding
     * box: x_i \in (0,1)
     */

    /*
     * The variable currentlev will contain the current level.
     * We start from level 0.
     */
    int currentlev = 0;

    short int edge = 0;
    while (ipoi != 0)
    {
        // Get the current dimension.
        int id = searchdim (currentlev, dimt);

        /*
         * We take advantage of the fact that we are only descending
         * the tree. Then we recursively multiply by 2 the coordinate.
         */
        x[id] *= 2.;
        ifth = ipoi;

        if (x[id] < 1.)
            // Go to the left.
        {
            edge = 0;
        }
        else
        {
            // Go to the right.
            edge = 1;
            --x[id];
        }
        // Next level.
        ++currentlev;
        ipoi = _data[ifth].getchild (edge);
    }

    if ( iava == iend )
    {
        /*
         * The list of available node is empty, so we have to put the new node
         * in the yet unassigned portion of the vector storing the tree.
         */
        _data.push_back (node);
    }
    else
    {
        _data[iava] = node;
    }

    // Add the node in the next available location.
    ipoi = iava;
    // Set the father link to the new node.
    _data[ifth].setchild (edge, ipoi);
    // iava is the next available location.
    iava = _data[ipoi].getchild (0);
    ++nele;
    if (iava == 0)
    {
        if ( iend > _header.gettreeloc() )
        {
            // Not enough space.
            throw TreeAlloc<Shape>();
        }
        ++iend;
        iava = iend;
    }

    /*
     * Set left_link and right link of ipoi equal to 0: this node
     * is a leaf. This is already done by the constructor of TreeNode class
     * when we put the node in the "free store".
     */
    _data[ipoi].setchild (0, 0);
    _data[ipoi].setchild (1, 0);
    _data[ipoi].setfather (ifth);

    // Store back header informations.
    _header.setiend (iend);
    _header.setiava (iava);
    _header.setnele (nele);
    if (currentlev > _header.gettreelev() )
    {
        _header.settreelev (currentlev);
    }

    return ipoi;
}

int ADTree::addtreenode (TreeNode const& node)
{
    try
    {
        int iloc = adtrb (node);
        return iloc;
    }
    catch (TreeAlloc<TreeNode::Shape>)
    {
        // Handle a TreeAlloc exception.
        std::cout << "warning! not enough space" << std::endl;
        std::cout << "increasing tree memory locations up to 1.5 * number of current locations..." << std::endl;
        int locv = _header.gettreeloc();
        int delta = int (locv / 2);
        _header.settreeloc (locv + delta);
        _data.resize (_header.gettreeloc() + 1);
        int iloc = adtrb (node);
        return iloc;
    }
}

  bool ADTree::search (BBox<3> region, std::vector<int>& found) const
{
    // Start preorder traversal at level 0 (root).
    int ipoi = _data[0].getchild (0);
    int ipoiNext = 0;
    int dimp = _header.getndimp();
    int dimt = _header.getndimt();

    std::vector<double> const& origin (_header.domainorig() );
    std::vector<double> const& scale (_header.domainscal() );
    // scale region
    region.transform (& (origin[0]), & (scale[0]) );

    // xl is the origin point for searching.
    std::vector<double> xl (dimt, 0.0);


    std::stack<ADTree::Helper> _stack;
    ADTree::Helper helper;
    helper.xl.resize (dimt);

    found.clear();

    int lev = 0;
    
    // Check if the region is external to the domain
    for (int i = 0; i < dimp; ++i)
      {
	if(region[i]>1.0 || region[i+dimp]<0.0)
	  {
            return false;
	  }
      }
    
    /*
     * To traverse a non-empty binary tree in preorder,
     * perform the following operations recursively at each node, starting with the root node:
     * > visit the root;
     * > traverse the left subtree;
     * > traverse the right subtree.
     */
    while (ipoi != 0)
    {
        do
        {
            BBox<3> const& xel = _data[ipoi].boundingBox();
            bool flag = region.intersect (xel);
            if (flag)
            {
                found.push_back (ipoi);
            }

            // Traverse left subtree.
            ipoiNext = _data[ipoi].getchild (0);

            // Check if subtree intersects box.
            if (ipoiNext != 0)
            {
                int id = searchdim (lev, dimt);
                double amov = delta (lev, dimt);

                if (id < dimp)
                {
                    if (xl[id] > region[id + dimp])
                    {
                        ipoiNext = 0;
                    }
                }
                else
                {
                    if (xl[id] + amov < region[id - dimp])
                    {
                        ipoiNext = 0;
                    }
                }
            }

            /*
             * Left subtree is neither null nor 'external'.
             * Push info onto the stack.
             */
            if (ipoiNext != 0)
            {
                helper.ipoi = ipoi;
                helper.xl = xl;
                helper.lev = lev;
                _stack.push (helper);
                ipoi = ipoiNext;
                ++lev;
            }
        }
        while (ipoiNext != 0);

        // Traverse right subtree.
        ++lev;
        ipoi = _data[ipoi].getchild (1);
        do
        {
            while (ipoi == 0)
            {
                // If right_link is null we have to get the point from the stack.
                if (_stack.empty() )
                {
                    return !found.empty();
                }
                helper = _stack.top();
                _stack.pop();
                ipoi = _data[helper.ipoi].getchild (1);
                xl = helper.xl;
                lev = helper.lev+1;
            }

            /*
             * Check if the subtree intersects box. Otherwise set right_link = 0,
             * pop new node from the stack and adjourn level.
             *
             * lev-1 is the level of the parent node, which directs the search.
             */
            int id = searchdim (lev - 1, dimt);
            double amov = delta (lev - 1, dimt);
	    
            // Reset xl (we have just traversed the tree to the right).
            xl[id] += amov;

            // Check if the subtree intersects box.
            double x1 = xl[id];
            double x2 = xl[id] + amov;
            if (id < dimp)
            {
                if (x1 > region[id + dimp])
                {
                    ipoi = 0;
                }
            }
            else
            {
                if (x2 < region[id - dimp])
                {
                    ipoi = 0;
                }
            }
        }
        while (ipoi == 0);

    }

    return !found.empty();
}


void ADTree::deltreenode (int const& index)
{
    if (index < _header.getiend() )
    {
        int nele = _header.getnele();
        int ipoi = index;
        int iava = _header.getiava();
        int ifth = _data[ipoi].getfather();

        --nele;
        int lchild = _data[ipoi].getchild (0);
        int rchild = _data[ipoi].getchild (1);
        int whichchild;

        if (_data[ifth].getchild (0) == ipoi)
        {
            whichchild = 0;
        }
        else
        {
            whichchild = 1;
        }

        if ( (lchild == 0) && (rchild == 0) )
        {
            // If the location is already a leaf, we can just release it.
            _data[ifth].setchild (whichchild, 0);
        }
        else
        {
            int ipoiNext = ipoi;
            int flag = 0;
            int il = _data[ipoiNext].getchild (0);
            int ir = _data[ipoiNext].getchild (1);
            // Traverse the tree to find an empty leaf.
            while ( (il != 0) || (ir != 0) )
            {
                il = _data[ipoiNext].getchild (0);
                ir = _data[ipoiNext].getchild (1);
                if (il != 0)
                {
                    flag = 0;
                    ipoiNext = il;
                }
                else if (ir != 0)
                {
                    flag = 1;
                    ipoiNext = ir;
                }
            }
            _data[ifth].setchild (whichchild, ipoiNext);
            int foo_f = _data[ipoiNext].getfather();
            _data[foo_f].setchild (flag, 0);
            if (ipoiNext != lchild)
            {
                _data[ipoiNext].setchild (0, lchild);
            }
            if (ipoiNext != rchild)
            {
                _data[ipoiNext].setchild (1, rchild);
            }
        }

        // Add the erased location to iava
        _data[ipoi].setchild (0, iava);
        _data[ipoi].setchild (1, 0);

        // Store back header informations.
        _header.setiava (ipoi);
        _header.setnele (nele);
    }
}

std::ostream& operator<< (std::ostream& ostr, ADTree const& myadt)
{
    ostr << myadt._header;

    return ostr;
}
}// end namespace ADT
