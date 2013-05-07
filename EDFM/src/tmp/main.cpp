#include "domain.hpp"
#include <iostream>
#include <sstream>
#include <string>
#include "adtree.hpp"
#include "implHypParab.hpp"
#include "geomCPgrid.hpp"
#include <stdexcept>
int main()
{
    using namespace std;
    using namespace ADT;
    using namespace Geometry;
    string gridfile ("../../data/EP_3D_grid.GRDECL");
    CPgrid grid (gridfile);
    /*
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
      grid.cell (1, 1, 1);
    const std::vector<Geometry::Point3D>
      cellVerticesa = cell.getVerticesVector();
    cell=grid.cell (2, 2, 2);
    std::vector<Geometry::Point3D>
      cellVertices = cell.getVerticesVector();
    cellVertices.insert(cellVertices.end(),cellVerticesa.begin(),cellVerticesa.end());
    numVertices = cellVertices.size();
    std::cout<<"Num vertices="<<numVertices<<std::endl;
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
    std::cout<<box<<std::endl;
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
    
    Point3D P00(0,0,0);
    Point3D P01(0,1,0);
    Point3D P11(1,1,0);
    Point3D P10(1,0,0);
    bilinerarPatchPoints patch={P00,P10,P11,P10};
    */
    double x(0),y(0),z(0);
    std::vector<UInt> sol(3);
    grid.cell(3,2,1).showMe(std::cout);
    std::vector<Point3D> points;
    points.push_back(grid.cell(3,3,3).getVertex(2));
    points.push_back(grid.cell(3,3,3).getVertex(1));
    points.push_back(grid.cell(4,2,1).getVertex(5));
    points.push_back(grid.cell(3,3,6).getVertex(2));
    points.push_back(grid.cell(3,3,3).getVertex(8));
    //cout<<pippo.x<<" "<<pippo.y<<" "<<pippo.z<<endl;

    //    while(x!=-1. || y!=-1. || z!= -1.)
    for(auto i : points)
      {
	//	cout<<" Give Me x y and z (-1 -1 -1 to end)"<<endl;
	//cin>>x>>y>>z;
	x=i.x;
	y=i.y;
	z=i.z;
	cout <<"target"<<x<<" "<<y<<" "<<z<<endl;
	try{
	  grid.whereIs(Point3D(x,y,z),sol);
	  }
	catch (std::runtime_error const & e)
	  {
	    cout<<"POINT NOT IN GRID"<<std::endl;
	  }
	cout<<sol[0]<<" "<<sol[1]<<" "<<sol[2]<<endl;
      }
    
}

