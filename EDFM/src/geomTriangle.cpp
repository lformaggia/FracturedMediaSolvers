 /*!
 *	@file geomTriangle.cpp
 *	@brief Class for Triangle in 3D space (definition). 
 *
 *	@author Luca Turconi <lturconi@gmail.com>
 *  @date 22-09-2012
 *
 */

#include<cmath>
#include<limits>
#include<fstream>
#include<iomanip>
#include<vector>
#include<algorithm>
 
#include "geomTriangle.hpp"

namespace Geometry
{
	
	// --------------------   Class Triangle   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================
	Triangle::Triangle() : M_pA(), M_pB(), M_pC(), M_Lmax() {}
	
	Triangle::Triangle(const Point3D & a, const Point3D & b, const Point3D & c) :
		M_pA(a), M_pB(b), M_pC(c)
	{ this->setLmax(); }
	
	Triangle::Triangle(const Triangle & t) : M_pA(t.A()), M_pB(t.B()), M_pC(t.C())
	{ this->setLmax(); }
	
	Triangle::~Triangle() {}

// ==================================================
// Set Methods
// ==================================================	
	void Triangle::setLmax()
	{
		std::vector<Real> L;
		L.push_back( (M_pA-M_pB).norm() );
		L.push_back( (M_pA-M_pC).norm() );
		L.push_back( (M_pC-M_pB).norm() );
		
		M_Lmax = *std::max_element(L.begin(), L.end());
	}
	
// ==================================================
// Methods
// ==================================================	
	bool Triangle::containPoint(const Point3D & p) const
	{
		if( !this->coplanarWithPoint(p) )
		{ return 0; }
		
		Real s(this->area());
		
		Triangle tri1(p,M_pA,M_pB);
		Real s1(tri1.area());

		Triangle tri2(p,M_pA,M_pC);
		Real s2(tri2.area());

		Triangle tri3(p,M_pB,M_pC);
		Real s3(tri3.area());
		
		return ( (s1+s2+s3)/s <= 1+eps );
	}
	
	bool Triangle::isIntersectedBy(const Segment & s) const
	{
		if( !s.intersectThePlaneOf(*this) )
			return 0;
		
		if( !this->containPoint( s.intersectionWithThePlaneOf(*this) ) )
			return 0;
		
		return 1;
	}
	
	Point3D Triangle::intersectionWith(const Segment & s) const
	{
		Point3D p(std::numeric_limits<Real>::quiet_NaN(),
					  std::numeric_limits<Real>::quiet_NaN(),
					  std::numeric_limits<Real>::quiet_NaN() );
		
		if( !s.intersectThePlaneOf(*this) )
			return p;
		
		Point3D inter;
		inter = s.intersectionWithThePlaneOf(*this);
		
		if( !this->containPoint(inter) )
			return p;

		return inter;
	}

	std::vector<Point3D> Triangle::getGaussNodes(UInt deg){
		Point3D p;
		std::vector<Point3D> nodigauss;
		if (deg==2){
		p=2./3.*M_pA + 1./6.*M_pB + 1./6.*M_pC; 
		nodigauss.push_back(p);
		p=1./6.*M_pA + 2./3.*M_pB + 1./6.*M_pC;
		nodigauss.push_back(p);
		p=1./6.*M_pA +1./6.*M_pB +2./3.*M_pC;
		nodigauss.push_back(p);
		}
		if (deg==4){
		p=0.816847572*M_pA + 0.091576213*M_pB + 0.091576213*M_pC; 
		nodigauss.push_back(p);
		p= 0.091576213*M_pA + 0.816847572*M_pB +  0.091576213*M_pC;
		nodigauss.push_back(p);
		p= 0.091576213*M_pA + 0.091576213*M_pB +0.816847572*M_pC;
		nodigauss.push_back(p);
         
		p=0.1081030181*M_pA + 0.44594849*M_pB +0.44594849*M_pC; 
		nodigauss.push_back(p);
		p= 0.44594849*M_pA + 0.1081030181*M_pB + 0.44594849*M_pC;
		nodigauss.push_back(p);
		p= 0.44594849*M_pA +0.44594849*M_pB +0.1081030181*M_pC;
		nodigauss.push_back(p);
		}
	
		return nodigauss;
	}

	std::vector<Real> Triangle::getGaussWeights(UInt deg){	
		std::vector<Real>  pesigauss;
		if (deg==2){
		pesigauss.push_back(1./3.);
		pesigauss.push_back(1./3.);
		pesigauss.push_back(1./3.);
		}
		if (deg==4){
		pesigauss.push_back(0.10995174);
		pesigauss.push_back(0.10995174);
		pesigauss.push_back(0.10995174);
		pesigauss.push_back(0.2233815896);
		pesigauss.push_back(0.2233815896);
		pesigauss.push_back(0.2233815896);
		}
		return pesigauss;
		
	}
	
	bool Triangle::exportVtk(const std::string & fileName) const
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
		
		std::cout << std::endl << " Exporting triangle in Vtk format... " << std::endl;
		
		UInt nCells = 1;
		UInt nPoints = 3;
		UInt CellType = 5; // for VTK_HEXAHEDRON
		
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
		filestr << std::endl;
		
		// Celldata
		filestr << "CELLS " << nCells << " " << 4 << std::endl;
		filestr << "3 0 1 2" << std::endl;
		filestr << std::endl;
		
		filestr << "CELL_TYPES " << nCells << std::endl;
		filestr << CellType << std::endl;
		filestr << std::endl;
		
		filestr.close();
		
		return 1;
	}
	
	void Triangle::showMe(std::ostream  & out) const
	{
		out << "Type = Triangle : " << std::endl;
		out << " M_pA : ";
		M_pA.showMe();
		out << " M_pB : ";
		M_pB.showMe();
		out << " M_pC : ";
		M_pC.showMe();
		out << " Area = " << this->area() << std::endl;
		out << " Perimeter = " << this->perimeter() << std::endl;
		out << " Normal = " << this->normal() << std::endl;
	}
	
// ==================================================
// Operators
// ==================================================
	std::ostream& operator<<(std::ostream & ostr, const Triangle & t)
	{
		ostr << " [ " << t.A() << " , "
				<< t.B() << " , "
				<< t.C() << " ]" << std::flush;
		return ostr;
	}


} // namespace Geometry
