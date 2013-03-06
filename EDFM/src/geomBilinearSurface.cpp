 /*!
 *	@file geomBilinearSurface.cpp
 *	@brief Class for Bilinear Surface in 3D space (definition).
 *
 *	@author Luca Turconi <lturconi@gmail.com>
 *  @date 22-09-2012
 *
 */

#include<cmath>
#include<limits>
#include<fstream>
#include<algorithm>
#include <iomanip>

#include<omp.h>

#include "geomBilinearSurface.hpp"

namespace Geometry
{
	
	// --------------------   Class BilinearSurface   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================
	BilinearSurface::BilinearSurface() : M_pA(), M_pB(), M_pC(), M_Lmax(),
		M_Triangulated(0), M_TriangulationToll(), M_TriangulatedSurface() {}
	
	BilinearSurface::BilinearSurface(const Point3D & a, const Point3D & b,
									 const Point3D & c, const Point3D & d ) :
		M_pA(a), M_pB(b), M_pC(c), M_pD(d),
		M_Triangulated(0), M_TriangulationToll(), M_TriangulatedSurface()
	{ this->setLmax(); }
		
	BilinearSurface::BilinearSurface(const BilinearSurface & b) :
		M_pA(b.A()), M_pB(b.B()), M_pC(b.C()), M_pD(b.D()),
		M_Triangulated(0), M_TriangulationToll(), M_TriangulatedSurface()
	{ this->setLmax(); }
	
	BilinearSurface::~BilinearSurface() {}

// ==================================================
// Set Methods
// ==================================================
	void BilinearSurface::setLmax()
	{
		std::vector<Real> L;
		L.push_back( (M_pA-M_pB).norm() );
		L.push_back( (M_pA-M_pC).norm() );
		L.push_back( (M_pA-M_pD).norm() );
		L.push_back( (M_pB-M_pC).norm() );
		L.push_back( (M_pB-M_pD).norm() );
		L.push_back( (M_pC-M_pD).norm() );
		
		M_Lmax = *std::max_element(L.begin(), L.end());
	}
	
// ==================================================
// Methods
// ==================================================

	
	void BilinearSurface::rescaleZ(Real scale, Real shift){
		M_pA.z=scale*(M_pA.z-shift)+shift;
		M_pB.z=scale*(M_pB.z-shift)+shift;
		M_pC.z=scale*(M_pC.z-shift)+shift;
		M_pD.z=scale*(M_pD.z-shift)+shift;

	}

	std::vector<Real> BilinearSurface::inv_param(const Point3D p, const Point3D normal){
		std::vector<Real> uv(3,0);
		uv[0]=0.5;
		uv[1]=0.5;
		Segment s(p+normal,p-normal);
					
		this->newtonIntersectionWith_uv(s,uv);
		
		return uv;
	}

	void BilinearSurface::makePlanar(){
		Point3D c, v1,v2,v3,v4,n;	
		c=0.25*(M_pA+M_pB+M_pC+M_pD);
		v1=M_pA-c;
		v2=M_pB-c;
		v3=M_pC-c;
		v4=M_pD-c;
		n=0.5*v1.cross(v2) + 0.5*v3.cross(v4);
		n=n/n.norm();
		M_pA=M_pA - (n.dot(M_pA-c))*n;
		M_pB=M_pB - (n.dot(M_pB-c))*n;
		M_pC=M_pC - (n.dot(M_pC-c))*n;
		M_pD=M_pD - (n.dot(M_pD-c))*n;

	}

	Real BilinearSurface::areaFault(){
		Point3D lato1(M_pA),lato2(M_pB), lato3(M_pC), lato4(M_pD), NN1,NN2;
		lato1=lato1-M_pB;
		lato2=lato2-M_pC;
		lato3=lato3-M_pD;
		lato4=lato4-M_pA;
		NN1=lato1.cross(lato2);
		NN2=lato3.cross(lato4);
		Real Area;
		Area=0.5*NN1.norm()+0.5*NN2.norm();
		return Area;
		
	}

	Real BilinearSurface::maxCurvature(const Real & u, const Real & v) const
	{
		// Coefficients of I fundamental form
		Real E( this->S_u(u,v).dot(this->S_u(u,v)) );
		Real F( this->S_u(u,v).dot(this->S_v(u,v)) );
		Real G( this->S_v(u,v).dot(this->S_v(u,v)) );
		 
		// Coefficients of II fundamental form
// 		Real L(0);
		Real M( this->S_uv().dot(this->normal(u,v)) );
// 		Real N(0);
		
		// Gaussian curvature
// 		Real K( (L*N - M*M)/(E*G - F*F) );
		Real K( -M*M/(E*G - F*F) );
		
		// Mean Curvature
// 		Real H( (E*N + G*L + 2*F*M)/(2*(E*G - F*F)) );
		Real H( 2*F*M/(2*(E*G - F*F)) );
		
		// Principal curvature
		Real k1( H + std::sqrt(H*H - K) );
		Real k2( H - std::sqrt(H*H - K) );

		if(std::fabs(k1)>std::fabs(k2))
		{
			return std::fabs(k1);
		}
		
		return std::fabs(k2);
	}
	
	Real BilinearSurface::maxCurvature() const
	{
		Real k;
		Real kMax(0);
		
		for(Real u=0; u<=1; u+=0.01)
		{
			for(Real v=0; v<=1; v+=0.01)
			{
				k = this->maxCurvature(u,v);
				if(k>kMax)
					kMax = k;
			}
		}
		
		return kMax;
	}
	
	void BilinearSurface::buildTriangulation(const Real & toll)
	{
		if( toll==M_TriangulationToll && M_Triangulated )
			return;
		
		M_TriangulatedSurface.clear();
		std::cout << "maxCurv = " << this->maxCurvature() << std::endl;
		Real h(std::sqrt(2*toll*M_Lmax/this->maxCurvature()));
		
		std::cout << " buildTriangulation: START" << std::endl;
		std::cout << "    h = " << h << std::endl;
		
		
		std::cout << "    Lmax = " << M_Lmax << std::endl;
		
		UInt nElemPerEdge = M_Lmax/h + 1;
		std::cout << "    nElemPerEdge = " << nElemPerEdge << std::endl;

		M_TriangulatedSurface.resize(nElemPerEdge*nElemPerEdge*2);
		
		#pragma omp parallel shared(nElemPerEdge)
		{
			Triangle t;
			Real u(0), v(0);
			UInt id(0);
			#pragma omp for schedule(dynamic,1) private(u,v,id,t)
			for(UInt i=0; i<nElemPerEdge; ++i)
			{
				u = i * 1./nElemPerEdge;
				v = 0.;
				for(UInt j=0; j<nElemPerEdge; ++j)
				{
				//std::cout << " u = " << u << "  v = " << v << std::endl;
					t.setA( this->param(u,v) );
					t.setB( this->param(u,v+1./nElemPerEdge) );
					t.setC( this->param(u+1./nElemPerEdge,v) );
					t.setLmax();
					id = 2*i*nElemPerEdge + 2*j;
					#pragma omp critical
					{ M_TriangulatedSurface[id] = t; }
					
					t.setA( this->param(u,v+1./nElemPerEdge) );
					t.setB( this->param(u+1./nElemPerEdge,v+1./nElemPerEdge) );
					t.setC( this->param(u+1./nElemPerEdge,v) );
					t.setLmax();
					id += 1;
					#pragma omp critical
					{ M_TriangulatedSurface[id] = t; }
					
					v += 1./nElemPerEdge;
				}
			}	
				
		}
		std::cout << " buildTriangulation: DONE" << std::endl;
		M_Triangulated = 1;
		M_TriangulationToll = toll;
	}
	
	
	bool BilinearSurface::approxIsIntersectedBy(const Segment & s, const bool & stdDivision) const
	{
		// split bilinear surface in two triangles (two different ways)
		Triangle t1(M_pA,M_pB,M_pC);
		Triangle t2(M_pA,M_pC,M_pD);
		
		if(!stdDivision)
		{
			Triangle t(M_pA,M_pB,M_pD);
			t1=t;
			Triangle tt(M_pB,M_pC,M_pD);
			t2=tt;
		}
		
		// Test
		if( t1.isIntersectedBy(s) || t2.isIntersectedBy(s) )
			return 1;
		
		return 0;
	}

	Point3D BilinearSurface::newtonIntersectionWith(const Segment & s, const Real & toll,
			const UInt & maxIter) const
	{
// 		std::cout << "NEWTON: " << std::endl;
// 		std::cout << " ------------------------------ " << std::endl;
// 		std::cout << " segment: " << s << std::endl;
		
		Point3D inter;
		
		
		// Newton Method1: global problem, without constraints imposition
		// ----------------------------------------------------------------
		std::vector<Real> res(3);
		res[0]=1;
		std::vector<Real> dx(3,0), xNew(3,0), x(3, 0.5);
		// x = (u,v,l) : CI x = (0.5,0.5,0.5)
		std::vector<Real> b(3);
		
		// columns matrix A = Jacobian matrix
		std::vector<Real> a1(3); // first column
		std::vector<Real> a2(3);
		std::vector<Real> a3(3);
		Real normRes(toll+1);
		
		Real detA, detA1, detA2, detA3;
		
		UInt it=0;
		
// 		std::cout << " Newton START " << std::endl;
		// Newton Method
		while( normRes>toll && it<maxIter)
		{
				// Build Matrix
			a1[0] = (1-x[1])*(M_pB.x-M_pA.x) + x[1]*(M_pC.x-M_pD.x);
			a2[0] = (1-x[0])*(M_pD.x-M_pA.x) + x[0]*(M_pC.x-M_pB.x);
			a3[0] = s.A().x-s.B().x;

			a1[1] = (1-x[1])*(M_pB.y-M_pA.y) + x[1]*(M_pC.y-M_pD.y);
			a2[1] = (1-x[0])*(M_pD.y-M_pA.y) + x[0]*(M_pC.y-M_pB.y);
			a3[1] = s.A().y-s.B().y;

			a1[2] = (1-x[1])*(M_pB.z-M_pA.z) + x[1]*(M_pC.z-M_pD.z);
			a2[2] = (1-x[0])*(M_pD.z-M_pA.z) + x[0]*(M_pC.z-M_pB.z);
			a3[2] = s.A().z-s.B().z;

				// Build b = -F
			b[0] = s.param(x[2]).x - this->param(x[0],x[1]).x;
			b[1] = s.param(x[2]).y - this->param(x[0],x[1]).y;
			b[2] = s.param(x[2]).z - this->param(x[0],x[1]).z;

			// Solving linear sistem (Cramer's rule)
			
			detA = a1[0]*a2[1]*a3[2] + a1[2]*a2[0]*a3[1] + a1[1]*a2[2]*a3[0]
					- a1[2]*a2[1]*a3[0] - a1[0]*a2[2]*a3[1] - a1[1]*a2[0]*a3[2];
			
			detA1 = b[0]*a2[1]*a3[2] + b[2]*a2[0]*a3[1] + b[1]*a2[2]*a3[0]
					- b[2]*a2[1]*a3[0] - b[0]*a2[2]*a3[1] - b[1]*a2[0]*a3[2];
			
			detA2 = a1[0]*b[1]*a3[2] + a1[2]*b[0]*a3[1] + a1[1]*b[2]*a3[0]
					- a1[2]*b[1]*a3[0] - a1[0]*b[2]*a3[1] - a1[1]*b[0]*a3[2];
			
			detA3 = a1[0]*a2[1]*b[2] + a1[2]*a2[0]*b[1] + a1[1]*a2[2]*b[0]
					- a1[2]*a2[1]*b[0] - a1[0]*a2[2]*b[1] - a1[1]*a2[0]*b[2];
		
			dx[0] = detA1 / detA;
			dx[1] = detA2 / detA;
			dx[2] = detA3 / detA;
			
			// x + dx --> xNew
			xNew[0] = x[0] + dx[0];
			xNew[1] = x[1] + dx[1];
			xNew[2] = x[2] + dx[2];
			
			// xNew - x --> res
			res[0] = xNew[0] - x[0];
			res[1] = xNew[1] - x[1];
			res[2] = xNew[2] - x[2];
			
			normRes = std::sqrt( res[0]*res[0] + res[1]*res[1] + res[2]*res[2] );
			
			x = xNew;
			
			it++;
		}
		
// 		std::cout << " x[0]: " << x[0] << std::endl;
// 		std::cout << " x[1]: " << x[1] << std::endl;
// 		std::cout << " x[2]: " << x[2] << std::endl;
		
		if(it>=maxIter)
		{
// 			std::cout << " Newton not convergent: no intersection found!" << std::endl;
			inter.x = std::numeric_limits<Real>::quiet_NaN();
			inter.y = std::numeric_limits<Real>::quiet_NaN();
			inter.z = std::numeric_limits<Real>::quiet_NaN();
			return inter;
		}
	//commentato da anna per cercare anche le intersezioni virtuali
	
	/*	if( x[0]<(-toll) || x[0]>(1+toll) || x[1]<(-toll) || x[1]>(1+toll) || x[2]<(-toll) || x[2]>(1+toll) )
		{
// 			std::cout << " Newton: no intersection found!" << std::endl;
// 			inter = this->param(x[0],x[1]);
// 			std::cout << " Intersection: " << inter << std::endl;
			inter.x = std::numeric_limits<Real>::quiet_NaN();
			inter.y = std::numeric_limits<Real>::quiet_NaN();
			inter.z = std::numeric_limits<Real>::quiet_NaN();
			return inter;
		}*/
		
		// ----------------------------------------------------------------
		
// 		std::cout << " Newton END " << std::endl;
// 		std::cout << " -- Newton Iterations: " << it << std::endl;
		
		inter = this->param(x[0],x[1]);
// 		std::cout << " Intersection: " << inter << std::endl;
		
		return inter;
	}
	
	void BilinearSurface::newtonIntersectionWith_uv(const Segment & s,std::vector<Real> & uv, const Real & toll,
			const UInt & maxIter) const
	{
//		std::cout << "NEWTON: " << std::endl;
// 		std::cout << " ------------------------------ " << std::endl;
// 		std::cout << " segment: " << s << std::endl;
		
		Point3D inter;
	
		
		// Newton Method1: global problem, without constraints imposition
		// ----------------------------------------------------------------
		std::vector<Real> res(3);
		res[0]=1;
		std::vector<Real> dx(3,0), xNew(3,0), x(3, 0.5);
		// x = (u,v,l) : CI x = (0.5,0.5,0.5)
		std::vector<Real> b(3);
		
		// columns matrix A = Jacobian matrix
		std::vector<Real> a1(3); // first column
		std::vector<Real> a2(3);
		std::vector<Real> a3(3);
		Real normRes(toll+1);
		
		Real detA, detA1, detA2, detA3;
		
		UInt it=0;
		
// 		std::cout << " Newton START " << std::endl;
		// Newton Method
		while( normRes>toll && it<maxIter)
		{
				// Build Matrix
			a1[0] = (1-x[1])*(M_pB.x-M_pA.x) + x[1]*(M_pC.x-M_pD.x);
			a2[0] = (1-x[0])*(M_pD.x-M_pA.x) + x[0]*(M_pC.x-M_pB.x);
			a3[0] = s.A().x-s.B().x;

			a1[1] = (1-x[1])*(M_pB.y-M_pA.y) + x[1]*(M_pC.y-M_pD.y);
			a2[1] = (1-x[0])*(M_pD.y-M_pA.y) + x[0]*(M_pC.y-M_pB.y);
			a3[1] = s.A().y-s.B().y;

			a1[2] = (1-x[1])*(M_pB.z-M_pA.z) + x[1]*(M_pC.z-M_pD.z);
			a2[2] = (1-x[0])*(M_pD.z-M_pA.z) + x[0]*(M_pC.z-M_pB.z);
			a3[2] = s.A().z-s.B().z;

				// Build b = -F
			b[0] = s.param(x[2]).x - this->param(x[0],x[1]).x;
			b[1] = s.param(x[2]).y - this->param(x[0],x[1]).y;
			b[2] = s.param(x[2]).z - this->param(x[0],x[1]).z;

			// Solving linear sistem (Cramer's rule)
			
			detA = a1[0]*a2[1]*a3[2] + a1[2]*a2[0]*a3[1] + a1[1]*a2[2]*a3[0]
					- a1[2]*a2[1]*a3[0] - a1[0]*a2[2]*a3[1] - a1[1]*a2[0]*a3[2];
			
			detA1 = b[0]*a2[1]*a3[2] + b[2]*a2[0]*a3[1] + b[1]*a2[2]*a3[0]
					- b[2]*a2[1]*a3[0] - b[0]*a2[2]*a3[1] - b[1]*a2[0]*a3[2];
			
			detA2 = a1[0]*b[1]*a3[2] + a1[2]*b[0]*a3[1] + a1[1]*b[2]*a3[0]
					- a1[2]*b[1]*a3[0] - a1[0]*b[2]*a3[1] - a1[1]*b[0]*a3[2];
			
			detA3 = a1[0]*a2[1]*b[2] + a1[2]*a2[0]*b[1] + a1[1]*a2[2]*b[0]
					- a1[2]*a2[1]*b[0] - a1[0]*a2[2]*b[1] - a1[1]*a2[0]*b[2];
		
			dx[0] = detA1 / detA;
			dx[1] = detA2 / detA;
			dx[2] = detA3 / detA;
			
			// x + dx --> xNew
			xNew[0] = x[0] + dx[0];
			xNew[1] = x[1] + dx[1];
			xNew[2] = x[2] + dx[2];
			
			// xNew - x --> res
			res[0] = xNew[0] - x[0];
			res[1] = xNew[1] - x[1];
			res[2] = xNew[2] - x[2];
			
			normRes = std::sqrt( res[0]*res[0] + res[1]*res[1] + res[2]*res[2] );
			
			x = xNew;
			
			it++;
		}
		
// 		std::cout << " x[0]: " << x[0] << std::endl;
// 		std::cout << " x[1]: " << x[1] << std::endl;
// 		std::cout << " x[2]: " << x[2] << std::endl;
	inter = this->param(x[0],x[1]);
		if(it>=maxIter || fabs(x[0])>100 || fabs(x[1])>100 || fabs(x[2])>100 )
		{
 		//	std::cout << " Newton not convergent: no intersection found!" << std::endl;
			inter.x = std::numeric_limits<Real>::quiet_NaN();
			inter.y = std::numeric_limits<Real>::quiet_NaN();
			inter.z = std::numeric_limits<Real>::quiet_NaN();
			uv[2]=-1.;
			
		}
		//commentato da anna per cercare anche le intersezioni virtuali
	
		if( x[0]<(-toll) || x[0]>(1+toll) || x[1]<(-toll) || x[1]>(1+toll) )
		{
// 			std::cout << " Newton: no intersection found!" << std::endl;
// 			inter = this->param(x[0],x[1]);
// 			std::cout << " Intersection: " << inter << std::endl;
		//	inter.x = std::numeric_limits<Real>::quiet_NaN();
		//	inter.y = std::numeric_limits<Real>::quiet_NaN();
		//	inter.z = std::numeric_limits<Real>::quiet_NaN();
	 		uv[0]=x[0];
			uv[1]=x[1];
			uv[2]=1.;
			
		}
		
		// ----------------------------------------------------------------
		
// 		std::cout << " Newton END " << std::endl;
// 		std::cout << " -- Newton Iterations: " << it << std::endl;
		//commentato da anna per cercare anche le intersezioni virtuali
	
		if( x[0]>(-toll) && x[0]<(1+toll) && x[1]>(-toll) && x[1]<(1+toll) )
		{
		
// 		std::cout << " Intersection: " << inter << std::endl;
		uv[0]=x[0];
		uv[1]=x[1];
		uv[2]=0.0;
		
		}
		if (!s.isIn(inter)) {uv[2]=-1;}
		
	}
	
	Point3D BilinearSurface::approxIntersectionWith(const Segment & s, const bool & stdDivision) const
	{
		// split bilinear surface in two triangles (two different ways)
		Triangle t1(M_pA,M_pB,M_pC);
		Triangle t2(M_pA,M_pC,M_pD);
		
		if(!stdDivision)
		{
			Triangle t(M_pA,M_pB,M_pD);
			t1=t;
			Triangle tt(M_pB,M_pC,M_pD);
			t2=tt;
		}
		
		Point3D inter( t1.intersectionWith(s) );
		if( inter.x==inter.x && inter.y==inter.y && inter.z==inter.z )
			return inter;
		
		inter = t2.intersectionWith(s);
		return inter;
	}
	
	Point3D BilinearSurface::approxNewtonIntersectionWith(const Segment & s,
		 const Real & toll, const UInt & maxIter) const
	{
		Point3D p( std::numeric_limits<Real>::quiet_NaN(),
					   std::numeric_limits<Real>::quiet_NaN(),
					   std::numeric_limits<Real>::quiet_NaN() );
		
		if( !this->approxIsIntersectedBy(s,0) && !this->approxIsIntersectedBy(s,1) )
			return p;
		
		return this->newtonIntersectionWith(s,toll,maxIter);
	}
	
	Point3D BilinearSurface::approxRefinedIntersectionWith(const Segment & s, const Real & toll)
	{
		if(!M_Triangulated || toll != M_TriangulationToll)
			this->buildTriangulation(toll);
		
		Point3D p( std::numeric_limits<Real>::quiet_NaN(),
				   std::numeric_limits<Real>::quiet_NaN(),
				   std::numeric_limits<Real>::quiet_NaN() );
		
		std::vector<Triangle>::const_iterator it;
		for(it=M_TriangulatedSurface.begin();
			it!=M_TriangulatedSurface.end() && !( p.x==p.x && p.y==p.y && p.z==p.z );
			++it)
		{
			p = (*it).intersectionWith(s);
		}
		return p;
	}

	bool BilinearSurface::isIntersectedBy(BilinearSurface altra, std::vector<Point3D>& punti) const
	{	
		gmm::size_type cont(0);
		Segment s1(altra.A(),altra.B());
		std::vector<Real> uv1(3,0.);
		this->newtonIntersectionWith_uv(s1,uv1);

		if (uv1[0]>=0 && uv1[0]<=1 && uv1[1]>=0 && uv1[1]<=1 && uv1[2]>=0&& s1.isIn(this->param(uv1[0],uv1[1]))) {cont+=1; 
		punti.push_back(this->param(uv1[0],uv1[1]));}
		Segment s2(altra.B(),altra.C());
		std::vector<Real> uv2(3,0.); 
		this->newtonIntersectionWith_uv(s2,uv2);	
		if (uv2[0]>=0 && uv2[0]<=1 && uv2[1]>=0 && uv2[1]<=1&& uv2[2]>=0 && s2.isIn(this->param(uv2[0],uv2[1]))) {cont+=1;
		punti.push_back(this->param(uv2[0],uv2[1]));}
		Segment s3(altra.C(),altra.D());
		std::vector<Real> uv3(3,0.);
		this->newtonIntersectionWith_uv(s3,uv3);	
		if (uv3[0]>=0 && uv3[0]<=1 && uv3[1]>=0 && uv3[1]<=1&& uv3[2]>=0 && s3.isIn(this->param(uv3[0],uv3[1]))) {cont+=1;
		punti.push_back(this->param(uv3[0],uv3[1]));}
		Segment s4(altra.D(),altra.A());
		std::vector<Real> uv4(3,0.);
		this->newtonIntersectionWith_uv(s4,uv4);	
		if (uv4[0]>=0 && uv4[0]<=1 && uv4[1]>=0 && uv4[1]<=1&& uv4[2]>=0 && s4.isIn(this->param(uv4[0],uv4[1]))) {cont+=1;
		punti.push_back(this->param(uv4[0],uv4[1]));}				
		Segment s5(altra.D(),altra.B());
		std::vector<Real> uv5(3,0.);
		this->newtonIntersectionWith_uv(s5,uv5);	
		if (uv5[0]>=0 && uv5[0]<=1 && uv5[1]>=0 && uv5[1]<=1 && uv5[2]>=0&& s5.isIn(this->param(uv5[0],uv5[1]))) {cont+=1;
		punti.push_back(this->param(uv5[0],uv5[1]));}
		Segment s6(altra.C(),altra.A());
		std::vector<Real> uv6(3,0.);
		this->newtonIntersectionWith_uv(s6,uv6);	
		if (uv6[0]>=0 && uv6[0]<=1 && uv6[1]>=0 && uv6[1]<=1&& uv6[2]>=0 && s6.isIn(this->param(uv6[0],uv6[1]))) {cont+=1;
		punti.push_back(this->param(uv6[0],uv6[1]));}
		if (cont>0) {return true;}
		else {return false;}
	
	}
	
	bool BilinearSurface::exportVtk(const std::string & fileName) const
	{
		std::fstream filestr;
		
		filestr.open (fileName.c_str(), std::ios_base::out);
	
		if (filestr.is_open())
		{
			std::cout << std::endl << " File: " << fileName << ", successfully opened";
		}
		else
		{
			std::cerr << std::endl << " *** Error: file not opened *** " << std::endl << std::endl;
			return  0;
		}
		
		std::cout << std::endl << " Exporting bilinear surface in Vtk format... " << std::endl;
		
		UInt nCells = 1;
		UInt nPoints = 4;
		UInt CellType = 9; // for VTK_HEXAHEDRON
		
		// Header
		filestr << "# vtk DataFile Version 3.1" << std::endl;
		filestr << "this is a file created for Paraview" << std::endl;
		filestr << "ASCII" << std::endl;
		filestr << "DATASET UNSTRUCTURED_GRID" << std::endl;
		filestr << std::endl;	// The fifth line is empty.
		
		filestr << std::scientific << std::setprecision(10);
		
		// Pointdata
		filestr << "POINTS " << nPoints << " double" << std::endl;
		filestr << M_pA.x << " " << M_pA.y
				<< " " << M_pA.z << std::endl;
		filestr << M_pB.x << " " << M_pB.y
				<< " " << M_pB.z << std::endl;
		filestr << M_pC.x << " " << M_pC.y
				<< " " << M_pC.z << std::endl;
		filestr << M_pD.x << " " << M_pD.y
				<< " " << M_pD.z << std::endl;
		filestr << std::endl;
		
		// Celldata
		filestr << "CELLS " << nCells << " " << 5 << std::endl;
		filestr << "4 0 1 2 3" << std::endl;
		filestr << std::endl;
		
		filestr << "CELL_TYPES " << nCells << std::endl;
		filestr << CellType << std::endl;
		filestr << std::endl;
		
		filestr.close();
		
		return 1;
	}
	
	bool BilinearSurface::exportTriangulationVtk(const std::string & fileName) const
	{
		std::fstream filestr;
		
		filestr.open (fileName.c_str(), std::ios_base::out);
	
		if (filestr.is_open())
		{
			std::cout << std::endl << " File: " << fileName << ", successfully opened";
		}
		else
		{
			std::cerr << std::endl << " *** Error: file not opened *** " << std::endl << std::endl;
			return  0;
		}
		
		std::cout << std::endl << " Exporting bilinear surface triangulation in Vtk format... " << std::endl;
		
		UInt nCells = M_TriangulatedSurface.size();
		UInt nPoints = 3*nCells;
		UInt CellType = 5; // for VTK_TRIANGLE
		
		// Header
		filestr << "# vtk DataFile Version 3.1" << std::endl;
		filestr << "this is a file created for Paraview" << std::endl;
		filestr << "ASCII" << std::endl;
		filestr << "DATASET UNSTRUCTURED_GRID" << std::endl;
		filestr << std::endl;	// The fifth line is empty.
		
		filestr << std::scientific << std::setprecision(10);
		
		// Pointdata
		filestr << "POINTS " << nPoints << " double" << std::endl;
		for(std::vector<Triangle>::const_iterator it=M_TriangulatedSurface.begin();
			it != M_TriangulatedSurface.end(); ++it)
		{
			filestr << (*it).A().x << " " << (*it).A().y
				<< " " << (*it).A().z << std::endl;
			filestr << (*it).B().x << " " << (*it).B().y
				<< " " << (*it).B().z << std::endl;
			filestr << (*it).C().x << " " << (*it).C().y
				<< " " << (*it).C().z << std::endl;
		}
		filestr << std::endl;
		
		// Celldata
		filestr << "CELLS " << nCells << " " << 4*nCells << std::endl;
		for(UInt i=0; i < M_TriangulatedSurface.size(); ++i)
		{
			filestr << "3 " << 3*i << " " << 3*i+1 << " " << 3*i+2 << std::endl;
		}
		filestr << std::endl;
		
		filestr << "CELL_TYPES " << nCells << std::endl;
		for(UInt i=0; i < M_TriangulatedSurface.size(); ++i)
		{
			filestr << CellType << std::endl;
		}
		filestr << std::endl;
		
		filestr.close();
		
		return 1;
	}
	
	void BilinearSurface::showMe(std::ostream  & out) const
	{
		out << "Type = BilinearSurface : " << std::endl;
		out << " M_pA : ";
		M_pA.showMe();
		out << " M_pB : ";
		M_pB.showMe();
		out << " M_pC : ";
		M_pC.showMe();
		out << " M_pD : ";
		M_pD.showMe();
	}
	
// ==================================================
// Operators
// ==================================================
	std::ostream& operator<<(std::ostream & ostr, const BilinearSurface & b)
	{
		ostr << " [ " << b.A() << " , "
				<< b.B() << " , "
				<< b.C() << " , "
				<< b.D() << " ]" << std::flush;
		return ostr;
	}
	
} // namespace Geometry

