#include "domain.hpp"
#include <iostream>
#include <sstream>
#include <string>
#include "adtree.hpp"
int main()
{
    using namespace std;
    using namespace ADT;
    using namespace Geometry;
    string gridfile ("../../data/EP_3D_grid.GRDECL");
    CPgrid grid (gridfile);
    ADTree tree (grid);
    cout << " Tree has been created" << std::endl << std::flush;
    cout << tree << std::endl;
    //    double xmin[3]={9.598e+06,5.143e+06,4006};
    //    double xmax[3]={9.600e+06,5.148e+06,4200};
    //    BBox<3> box(xmin,xmax);
    int numVertices (0);
    double vmin[3];
    double vmax[3];
    Geometry::CPcell cell =
      grid.cell (3, 2, 1);
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
    std::vector<int> listFound;
    bool yes=tree.search(box,listFound);
    vector<int> keys(3);
    if(yes)
      {
	cout<<"Found "<<listFound.size()<<" possible intersections with box"<<endl;
	cout<<box<<endl;
	cout<<"Mesh cells potentially intersecting box"<<endl;
	for(auto i : listFound)
	  {
	    TreeNode const & node=tree.getNode(i);
	    keys=node.getkeys();
	    for (auto j : keys){
	      cout<<j<<" ";
	    }
	    cout<<endl;
	  }
      }
    else
      {
	cout<<"No intersecting cells"<<endl;
      }
}


