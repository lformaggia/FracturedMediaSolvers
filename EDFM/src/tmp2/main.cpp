#include "trilinearElement.hpp"
#include <iostream>
int main()
{
  using namespace std;
  using namespace Geometry;
  double x,y,z;
  
  std::vector<Point3D> points(8);
  points[4]=Point3D(0.0,0.0,0.0);
  points[5]=Point3D(2.0,0.0,0.0);
  points[6]=Point3D(0.5,2.3,0.2);
  points[7]=Point3D(2.0,2.0,0.0);
  points[0]=Point3D(0.0,0.0,1.0);
  points[1]=Point3D(1.0,0.0,1.0);
  points[2]=Point3D(0.0,1.0,1.0);
  points[3]=Point3D(2.0,2.0,1.0);
  CPcell cell(points);
  TrilinearElement element(cell);
  InverseMapping invmap(element);
  Point3D pippo(grid.cell(3,3,3).getVertex(2));
  cout<<pippo.x<<" "<<pippo.y<<" "<<pippo.z<<endl;
  while(x!=-1. || y!=-1. || z!= -1.)
    {
      cout<<" Give Me x y and z (-1 -1 -1 to end)"<<endl;
      cin>>x>>y>>z;
      InvMapResult image=invmap(x,y,z);
      cout<<image<<std::endl;

    }
}


