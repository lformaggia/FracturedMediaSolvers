#include "delaunay.hpp"
#include "geomPoint3D.hpp"
#include "geomBilinearSurface.hpp"
#include "geomSegment.hpp"
#include <iostream>

int main()
{
  using namespace std;
  using namespace Geometry;
  Point2D a(0.,0.);
  Point2D b(100., 0);
  Point2D c(100.,100.);
  Point2D d(0.,100.);
  Point2D e(50.,50.);
  std::vector<Point2D> points;
  points.push_back(a);
  points.push_back(b);
  points.push_back(c);
  points.push_back(d);
  points.push_back(e);
  std::vector<IdTriplet> elements;
  delaun(points,elements);
  for (unsigned int i=0;i<elements.size();++i)
    cout<<"Element n."<<i<<"Points: "<<elements[i][0]<<" "<<
      elements[i][1]<<" "<<elements[i][2]<<endl;
  std::vector<unsigned int> bPoints=
    orderedBoundaryPoints(elements);
  for (unsigned int i=0;i<bPoints.size();++i)
    {
      cout<<bPoints[i]<<" ";
    }
  cout<<endl;
  Aligned align;
  cout<<"Aligned? "<<align(a,b,c)<<endl;
  cout<<"Aligned? "<<align(a,b,a+b)<<endl;
  vector<Point2D> pippo;
  pippo.reserve(bPoints.size());
  for (unsigned int i=0;i<bPoints.size();++i)
    {
      pippo.push_back(points[bPoints[i]]);
    }
  pippo.push_back(pippo[3]);
  pippo[3]=0.5*pippo[4]+0.5*pippo[2];
  for (unsigned int i=0;i<pippo.size();++i)
    {
      cout<<pippo[i]<<endl;
    }
  pippo=align.decimate(pippo);
  cout<<"decimated"<<endl;
  for (unsigned int i=0;i<pippo.size();++i)
    {
      cout<<pippo[i]<<endl;
    }
  cout<<"Testing bilinear surfaces"<<endl;
  Point3D a3(0., 0., -10.);
  Point3D b3(10., 0., 0.);
  Point3D c3(10., 10., 0.);
  Point3D d3(0., 10., 10.);
  BilinearSurface surf(a3,b3,c3,d3);
  Point3D where=surf.param(0.3,0.1);
  vector<Real> uv;
  Point3D result=surf.projection(where,uv);
  cout<<where<<" "<<result<<endl;
  cout<<uv[0]<<" " <<uv[1]<<endl;
  cout<<" *************** TESTING SEGMENTS *******************"<<std::endl;
  Segment2D A(Point2D(0,0),Point2D(1,1));
  Segment2D B(Point2D(1.001,0),Point2D(1.001,10));
  Point2D Result;
  bool intersect=A.intersectTheSegment(B,Result);
  std::cout<<intersect<<std::endl;
  if(intersect)Result.showMe();
  cout<<A.isIn(Result)<<endl;
  cout<<A.isIn(Point2D(5.0,5.0))<<endl;
  

}
