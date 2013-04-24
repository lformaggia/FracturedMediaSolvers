 /*!
 *	@file geomCPgridElements.cpp
 *	@brief Classes for Corner Point Grid elements: pillar and cell (definition).
 *
 *	@author Luca Turconi <lturconi@gmail.com>
 *  @date 22-09-2012
 *
 */  

#include <algorithm>
#include <iomanip>
 
#include "geomCPgridElements.hpp"
#include "geomCPgrid.hpp"


namespace Geometry
{

	// --------------------   Class CPpillar   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================
	CPpillar::CPpillar() : Line(), M_i(0), M_j(0) {}

	CPpillar::CPpillar(const Real & xA, const Real & yA, const Real & zA,
					   const Real & xB, const Real & yB, const Real & zB,
					   const UInt & i, const UInt & j) :
		Line(xA,yA,zA,xB,yB,zB), M_i(i), M_j(j) {}
	
	CPpillar::CPpillar(const Point3D & a, const Point3D & b, const UInt & i, const UInt & j) :
		Line(a,b), M_i(i), M_j(j) {}
	
	CPpillar::~CPpillar() {}
	
// ==================================================
// Methods
// ==================================================
	Segment CPpillar::getLimitedPillar(const CPgrid & grid) const
	{
		Real zmax, zmin;
		std::vector<Real> zpil;
		
		UInt idx(0), idx1(0), idx2(0);
		UInt i(M_i);
		UInt j(M_j);
		
		for(UInt k=1; k<=grid.Nz(); ++k)
		{
			// Vertex 1 of cell (i,j.k)
			if( i<grid.Nx()+1 && j<grid.Ny()+1 )
			{
				idx1 = (i-1)*2 + (j-1)*2*(2*grid.Nx()) + (k-1)*2*grid.Ny()*(4*grid.Nx());
				zpil.push_back(grid.zcorn()[idx1]);
			}
			// Vertex 2 of cell (i-1,j,k)
			if( i>1 )
			{
				idx =  idx1 - 1;
				zpil.push_back(grid.zcorn()[idx]);
			}
			// Vertex 3 of cell (i,j-1.k)
			if( j>1 )
			{
				idx = idx1 - 2*grid.Nx();
				zpil.push_back(grid.zcorn()[idx]);
			}
			// Vertex 4 of cell (i-1,j-1,k)
			if( i>1 && j>1 )
			{
				idx -= 1;
				zpil.push_back(grid.zcorn()[idx]);
			}
			// Vertex 5 of cell (i,j,k)
			if( i<grid.Nx()+1 && j<grid.Ny()+1 )
			{
				idx2 = idx1 + (4*grid.Nx())*grid.Ny();
				zpil.push_back(grid.zcorn()[idx2]);
			}
			// Vertex 6 of cell (i-1,j,k)
			if( i>1 )
			{
				idx = idx2 - 1;
				zpil.push_back(grid.zcorn()[idx]);
			}
			// Vertex 7 of cell (i,j-1,k)
			if( j>1 )
			{
				idx = idx2 - 2*grid.Nx();
				zpil.push_back(grid.zcorn()[idx]);
			}
			// Vertex 8 of cell (i-1,j-1,k)
			if( i>1 && j>1 )
			{
				idx -= 1;
				zpil.push_back(grid.zcorn()[idx]);
			}
		}
		
		// @todo: Estrarre max e min dal vector!
		zmin = *std::min_element(zpil.begin(),zpil.end());
		zmax = *std::max_element(zpil.begin(),zpil.end());
			
		Point3D top( this->getPointAtZ(zmax) );
		Point3D bottom( this->getPointAtZ(zmin) );
		
		Segment pil(top, bottom);
		return pil;
	}

	void CPpillar::showMe(std::ostream  & out) const
	{
		out << "Type = CPpillar : " << std::endl;
		out << " M_pA : ";
		M_pA.showMe();
		out << " M_pB : ";
		M_pB.showMe();
		out << " (i,j) : Type = (UInt,UInt) : ( "
					<< M_i << " , " << M_j << " )" << std::endl;
	}
	
	
	// --------------------   Class CPcell   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================
	CPcell::CPcell() : M_vertex(), M_actnum(), M_i(0), M_j(0), M_k(0) {}
	
	CPcell::CPcell(const std::vector<Point3D> & v,
				   const UInt & i, const UInt & j, const UInt & k, const UInt & actnum) :
		M_vertex(v), M_actnum(actnum), M_i(i), M_j(j), M_k(k) {}
	
	CPcell::CPcell(const CPcell & c) : M_vertex(c.getVerticesVector()), M_actnum(c.getActnum()),
		M_i(c.i()), M_j(c.j()), M_k(c.k()) {}
	
	CPcell::~CPcell() {}

// ==================================================
// Get Methods
// ==================================================
	Point3D CPcell::getVertex(const UInt & i) const
	{
		if(i==0 || i>8)
		{
			std::cerr << " *** Error: invalid vertex number ***" << std::endl;
			Point3D p;
			return p;
		}
			
		return M_vertex[i-1];
	}
	
	Segment CPcell::getEdge(const UInt & i) const
	{
		Segment edge;
		
		if(i==0 || i>12)
		{
			std::cerr << " *** Error: invalid edge number ***" << std::endl;
			return edge;
		}
		
		switch (i)
		{
			case 1: edge = Segment(this->getVertex(1),this->getVertex(2));
				break;
			case 2: edge = Segment(this->getVertex(2),this->getVertex(4));
				break;
			case 3: edge = Segment(this->getVertex(4),this->getVertex(3));
				break;
			case 4: edge = Segment(this->getVertex(3),this->getVertex(1));
				break;
			case 5: edge = Segment(this->getVertex(1),this->getVertex(5));
				break;
			case 6: edge = Segment(this->getVertex(2),this->getVertex(6));
				break;
			case 7: edge = Segment(this->getVertex(4),this->getVertex(8));
				break;
			case 8: edge = Segment(this->getVertex(3),this->getVertex(7));
				break;
			case 9: edge = Segment(this->getVertex(5),this->getVertex(6));
				break;
			case 10: edge = Segment(this->getVertex(6),this->getVertex(8));
				break;
			case 11: edge = Segment(this->getVertex(8),this->getVertex(7));
				break;
			case 12: edge = Segment(this->getVertex(7),this->getVertex(5));
				break;
		}
		return edge;
	}
	
// ==================================================
// Methods
// ==================================================
	bool CPcell::exportVtk(const std::string & fileName) const
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
		
		std::cout << std::endl << " Exporting cell in Vtk format... " << std::endl;
		
		UInt nCells = 1;
		UInt nPoints = 8;
		UInt CellType = 12; // for VTK_HEXAHEDRON
		
		// Header
		filestr << "# vtk DataFile Version 3.1" << std::endl;
		filestr << "this is a file created for Paraview" << std::endl;
		filestr << "ASCII" << std::endl;
		filestr << "DATASET UNSTRUCTURED_GRID" << std::endl;
		filestr << std::endl;	// The fifth line is empty.
		
		filestr << std::scientific << std::setprecision(10);
		
		// Pointdata
		filestr << "POINTS " << nPoints << " double" << std::endl;
		for(UInt i=1; i<=8; ++i)
		{
			filestr << this->getVertex(i).x << " "
					<< this->getVertex(i).y << " "
					<< this->getVertex(i).z << std::endl;
		}
		filestr << std::endl;
		
		// Celldata
		filestr << "CELLS " << nCells << " " << 9 << std::endl;
		filestr << "8 4 5 7 6 0 1 3 2" << std::endl;
		filestr << std::endl;
		
		filestr << "CELL_TYPES " << nCells << std::endl;
		filestr << CellType << std::endl;
		filestr << std::endl;
		
		filestr.close();
		
		return 1;
	}


    void CPcell::showMe(std::ostream  & out) const
	{
		out << "Type = CPcell : " << std::endl;
		out << "            Vertex 1 : "; M_vertex[0].showMe();
		out << "   3----4   Vertex 2 : "; M_vertex[1].showMe();
		out << "  /|   /|   Vertex 3 : "; M_vertex[2].showMe();
		out << " 1----2 |   Vertex 4 : "; M_vertex[3].showMe();
		out << " | 7--|-8   Vertex 5 : "; M_vertex[4].showMe();
		out << " |/   |/    Vertex 6 : "; M_vertex[5].showMe();
		out << " 5----6     Vertex 7 : "; M_vertex[6].showMe();
		out << "            Vertex 8 : "; M_vertex[7].showMe();
		out << " M_actnum : Type = UInt : " << M_actnum << std::endl;
		out << " (i,j,k) : Type = (UInt,UInt,UInt) : ( "
			<< M_i << " , " << M_j << " , " << M_k <<" )" << std::endl;

	}

double CPcell::set_volume() const{
double v1,v2,v3,v4,v5;
Point3D lato1, lato2,lato3, provv;
std::vector<int> indici1, indici2, indici3, indici4, indici5;
indici1.push_back(0);indici1.push_back(1);indici1.push_back(2);indici1.push_back(4);
indici2.push_back(1);indici2.push_back(4);indici2.push_back(5);indici2.push_back(7);
indici3.push_back(6);indici3.push_back(7);indici3.push_back(4);indici3.push_back(2);
indici4.push_back(1);indici4.push_back(2);indici4.push_back(3);indici4.push_back(7);
indici5.push_back(4);indici5.push_back(7);indici5.push_back(1);indici5.push_back(2);

lato1=M_vertex[indici1[0]]-M_vertex[indici1[1]];
lato2=M_vertex[indici1[0]]-M_vertex[indici1[2]];
lato3=M_vertex[indici1[0]]-M_vertex[indici1[3]];
provv=lato1.cross(lato2);
v1=1./6.*lato3.dot(provv);

lato1=M_vertex[indici2[0]]-M_vertex[indici2[1]];
lato2=M_vertex[indici2[0]]-M_vertex[indici2[2]];
lato3=M_vertex[indici2[0]]-M_vertex[indici2[3]];
provv=lato1.cross(lato2);
v2=1./6.*lato3.dot(provv);

lato1=M_vertex[indici3[0]]-M_vertex[indici3[1]];
lato2=M_vertex[indici3[0]]-M_vertex[indici3[2]];
lato3=M_vertex[indici3[0]]-M_vertex[indici3[3]];
provv=lato1.cross(lato2);
v3=1./6.*lato3.dot(provv);

lato1=M_vertex[indici4[0]]-M_vertex[indici4[1]];
lato2=M_vertex[indici4[0]]-M_vertex[indici4[2]];
lato3=M_vertex[indici4[0]]-M_vertex[indici4[3]];
provv=lato1.cross(lato2);
v4=1./6.*lato3.dot(provv);

lato1=M_vertex[indici5[0]]-M_vertex[indici5[1]];
lato2=M_vertex[indici5[0]]-M_vertex[indici5[2]];
lato3=M_vertex[indici5[0]]-M_vertex[indici5[3]];
provv=lato1.cross(lato2);
v5=1./6.*lato3.dot(provv);
return fabs(v1)+fabs(v2)+fabs(v3)+fabs(v4)+fabs(v5);
}

bool CPcell::isIn(Point3D & punto1, Real toll ) const
{//ciclo sulle facce
int f1[4]={0,1,2,3};
int f2[4]={4,5,6,7};
int f3[4]={0,1,4,5};
int f4[4]={1,3,5,7};
int f5[4]={2,3,6,7};
int f6[4]={0,2,4,6};
double vol2(0);
double vol1=set_volume();
int* facce[6]={&f1[0],&f2[0],&f3[0],&f4[0],&f5[0],&f6[0] };
for (int f=0;f<6;++f){
Point3D lato1,lato2,lato3,provv;

lato1=punto1-M_vertex[facce[f][0]];
lato2=punto1-M_vertex[facce[f][1]];
lato3=punto1-M_vertex[facce[f][2]];
provv=lato1.cross(lato2);
vol2+=1./6.*fabs(provv.dot(lato3));
lato1=punto1-M_vertex[facce[f][1]];
lato2=punto1-M_vertex[facce[f][2]];
lato3=punto1-M_vertex[facce[f][3]];
provv=lato1.cross(lato2);
vol2+=1./6.*fabs(provv.dot(lato3));

}

if (fabs(vol1-vol2)<toll*vol1){
return true;}
else{
return false;}
}

bool CPcell::intersectTheFace(Segment S,  int quale,   Point3D & risultato) const
	{
		int f1[4]={0,1,3,2};
		int f2[4]={4,5,7,6};
		int f3[4]={0,1,5,4};
		int f4[4]={1,3,7,5};
		int f5[4]={2,3,7,6};
		int f6[4]={0,2,6,4};
		int* facce[6]={&f1[0],&f2[0],&f3[0],&f4[0],&f5[0],&f6[0] };
		
                BilinearSurface faccia( M_vertex[facce[quale][0]], M_vertex[facce[quale][1]], M_vertex[facce[quale][2]], M_vertex[facce[quale][3]] ); 
std::vector<Real> uv_point(3,0.); 

	//	risultato=faccia.approxIntersectionWith(S);
faccia.newtonIntersectionWith_uv(S,uv_point,1e-10,100);
			risultato=faccia.param(uv_point[0], uv_point[1]);
	       if (uv_point[2]==0) { return true;}
		return false;
 
	}

bool CPcell::segmentIntersectCell(Segment S, std::vector<Point3D> &vettpunti) const
{
	
	Point3D A(S.A()), B(S.B());
	for (UInt i=0; i<6;++i){
	Point3D ris;
	if (this->intersectTheFace(S, i, ris)) {vettpunti.push_back(ris); }
	}
	if (vettpunti.size()<2)
	{
		if (this->isIn(A)) vettpunti.push_back(A);
		if (this->isIn(B)) vettpunti.push_back(B);
	}
	if (vettpunti.size()>1){
	return true;
	}
	else
	{return false;}
}
	
} // namespace Geometry


