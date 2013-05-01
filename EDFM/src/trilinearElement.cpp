#include "trilinearElement.hpp"
#include <iostream>
# define DERIVATIVEX(DIR,DIRN,J)		\
  my_J(DIRN,J)=						\
      M_points[0].DIR*trilinearShapeDer<0,J>(xx,yy,zz)+	\
      M_points[1].DIR*trilinearShapeDer<1,J>(xx,yy,zz)+	\
      M_points[2].DIR*trilinearShapeDer<2,J>(xx,yy,zz)+	\
      M_points[3].DIR*trilinearShapeDer<3,J>(xx,yy,zz)+	\
      M_points[4].DIR*trilinearShapeDer<4,J>(xx,yy,zz)+	\
      M_points[5].DIR*trilinearShapeDer<5,J>(xx,yy,zz)+	\
      M_points[6].DIR*trilinearShapeDer<6,J>(xx,yy,zz)+	\
      M_points[7].DIR*trilinearShapeDer<7,J>(xx,yy,zz)


namespace Geometry
{
  //! Default values for inverse map algorithm
  const InvMapOption defaultMapOption=
    {1.e-6, 1.e-6,20,0.02};
  
  InvMapOption InverseMapping::M_options=defaultMapOption;
  
  std::ostream & operator <<(std::ostream & out, InvMapResult const & res)
  {
    out<<"Results of the inverse map iterations"<<std::endl;
    out<<std::boolalpha<<"Converged:"<<res.converged<<std::noboolalpha<<std::endl;
    out<<"In iterations: "<<res.numIter;
    out<<" Point found"<<res.finalPoint<<" is inside?: "<<std::boolalpha<<res.inside<<std::noboolalpha<<std::endl;
    return out;
  }

  TrilinearElement::TrilinearElement(const CPcell & cell): my_J(3,3)
  {
    M_points[0]=cell.getVertex(5);
    M_points[1]=cell.getVertex(6);
    M_points[2]=cell.getVertex(7);
    M_points[3]=cell.getVertex(8);
    M_points[4]=cell.getVertex(1);
    M_points[5]=cell.getVertex(2);
    M_points[6]=cell.getVertex(3);
    M_points[7]=cell.getVertex(4);
  }
   
  Point3D TrilinearElement::map(double const & xx, double const & yy, double const & zz)const
  {
    double x,y,z;
    x=
      M_points[0].x*trilinearShape<0>(xx,yy,zz)+
      M_points[1].x*trilinearShape<1>(xx,yy,zz)+
      M_points[2].x*trilinearShape<2>(xx,yy,zz)+
      M_points[3].x*trilinearShape<3>(xx,yy,zz)+
      M_points[4].x*trilinearShape<4>(xx,yy,zz)+
      M_points[5].x*trilinearShape<5>(xx,yy,zz)+
      M_points[6].x*trilinearShape<6>(xx,yy,zz)+
      M_points[7].x*trilinearShape<7>(xx,yy,zz);

    y=
      M_points[0].y*trilinearShape<0>(xx,yy,zz)+
      M_points[1].y*trilinearShape<1>(xx,yy,zz)+
      M_points[2].y*trilinearShape<2>(xx,yy,zz)+
      M_points[3].y*trilinearShape<3>(xx,yy,zz)+
      M_points[4].y*trilinearShape<4>(xx,yy,zz)+
      M_points[5].y*trilinearShape<5>(xx,yy,zz)+
      M_points[6].y*trilinearShape<6>(xx,yy,zz)+
      M_points[7].y*trilinearShape<7>(xx,yy,zz);

    z=
      M_points[0].z*trilinearShape<0>(xx,yy,zz)+
      M_points[1].z*trilinearShape<1>(xx,yy,zz)+
      M_points[2].z*trilinearShape<2>(xx,yy,zz)+
      M_points[3].z*trilinearShape<3>(xx,yy,zz)+
      M_points[4].z*trilinearShape<4>(xx,yy,zz)+
      M_points[5].z*trilinearShape<5>(xx,yy,zz)+
      M_points[6].z*trilinearShape<6>(xx,yy,zz)+
      M_points[7].z*trilinearShape<7>(xx,yy,zz);
    
    return Point3D(x,y,z);
  }

  gmm::dense_matrix<double> const &
  TrilinearElement::J(double const & xx, double const & yy, double const & zz)const
  {
    
    // X derivatives: \partial x/\partial xx_{j=0,2}
    DERIVATIVEX(x,0,0) ;
    DERIVATIVEX(x,0,1) ;	
    DERIVATIVEX(x,0,2) ;		
    
    // Y derivatives: \partial y/\partial xx_{j=0,2}
    DERIVATIVEX(y,1,0) ;
    DERIVATIVEX(y,1,1) ;	
    DERIVATIVEX(y,1,2) ;		
    
    // Z derivatives: \partial z/\partial xx_{j=0,2}
    DERIVATIVEX(z,2,0) ;
    DERIVATIVEX(z,2,1) ;	
    DERIVATIVEX(z,2,2) ;		
    
    return my_J;
  }
  
  InverseMapping::InverseMapping(TrilinearElement const & element):M_element(element)
  {
    M_BBmin[0]=M_BBmax[0]=M_element.M_points[0].x;
    M_BBmin[1]=M_BBmax[1]=M_element.M_points[0].y;
    M_BBmin[2]=M_BBmax[2]=M_element.M_points[0].z;
    for (unsigned int i=1;i<8;++i)
      {
	M_BBmin[0]=std::min(M_BBmin[0],M_element.M_points[i].x);
	M_BBmin[1]=std::min(M_BBmin[1],M_element.M_points[i].y);
	M_BBmin[2]=std::min(M_BBmin[2],M_element.M_points[i].z);
	M_BBmax[0]=std::max(M_BBmax[0],M_element.M_points[i].x);
	M_BBmax[1]=std::max(M_BBmax[1],M_element.M_points[i].y);
	M_BBmax[2]=std::max(M_BBmax[2],M_element.M_points[i].z);
      }
    // Expand bounding box by tolerance
    for (unsigned int i=0;i<3;++i)
      {
	double lenght=(M_BBmax[i]-M_BBmin[i])*InverseMapping::M_options.bboxTolerance;
	M_BBmin[i]-=lenght;
	M_BBmax[i]+=lenght;
      }
  }
  
  void  InverseMapping::F(const Point3D & in, const Point3D & target,Point3D & Fx) const
  {
    Fx=target-M_element.map(in.x,in.y,in.z);
  }
  
  InvMapResult InverseMapping::operator () (double const & x, double const & y, double const & z)
  const
  {
    InvMapResult res;
    // Eliminate trivial case: target outside bounding box.
    bool outside(false);
    double in[3]={x,y,z};
    for (unsigned int i=0;i<3;++i)
      {
	if(in[i] <M_BBmin[i])
	  {
	    outside=true;
	    res.direction[i]=-1;
	  }
	else if(in[i] >M_BBmax[i])
	  {
	    outside=true;
	    res.direction[i]=1;
	  }
	else res.direction[i]=0;
      }
    if (outside)
      {
	res.converged=false;
	res.numIter=0;
	res.inside=false;
	res.finalPoint=Point3D();
	return res;
      }
    unsigned int iteration=0;
    Point3D const target(x,y,z);
    Point3D lastValue(0.5,0.5,0.5);
    Point3D result(0.,0.,0.);
    Point3D Fx(0.,0.,0.);
    std::vector<double> fvalue(3);
    std::vector<double> delta(3);
    double res_iter=1.e+30;
    double res_resi=1.e+30;
    while(iteration++ < InverseMapping::M_options.maxIter
	  && res_iter > InverseMapping::M_options.iteratesTolerance
	  && res_resi > InverseMapping::M_options.residualTolerance
	  )
      {
	// COmpute residual (x-map(xx))
	this->F(lastValue,target,Fx);
	// Put into a vector to use gmm utilities
	fvalue[0]=Fx.x;
	fvalue[1]=Fx.y;
	fvalue[2]=Fx.z;
	res_resi=gmm::vect_norm2(fvalue);
	// Let's spare us a useless linear solution
	if(res_iter <= InverseMapping::M_options.iteratesTolerance)break;
	// COmpute Jacobian (infact is minus jacobian for Newton iterates!)
	gmm::lu_solve(
		      this->M_element.J(lastValue.x,lastValue.y,lastValue.z),
		      delta,
		      fvalue);
	res_iter=gmm::vect_norm2(delta);
	lastValue.x=lastValue.x+delta[0];
	lastValue.y=lastValue.y+delta[1];
	lastValue.z=lastValue.z+delta[2];
      }
    
    res.converged=
      res_resi <= InverseMapping::M_options.residualTolerance ||
      res_iter <= InverseMapping::M_options.iteratesTolerance;
    res.finalPoint=lastValue;
    res.numIter=iteration;
    
    in[0]=lastValue.x;
    in[1]=lastValue.y;
    in[2]=lastValue.z;
    // No tolernace, change if needed.
    // I want ot avoid a too large tolerance that would lead 
    // to false positive
    const double tolerance=1.e-8;
    for (unsigned int i=0;i<3;++i)
      {
	if(in[i] <-tolerance)
	  {
	    outside=true;
	    res.direction[i]=-1;
	  }
	else if(in[i] >1.+tolerance)
	  {
	    outside=true;
	    res.direction[i]=1;
	  }
	else res.direction[i]=0;
      }
    res.inside=!outside;
    return res;
  }
}
