#include "delaunay.hpp"
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
}
