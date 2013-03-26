 /*!
 *	@file geomFault.cpp
 *	@brief Class for Fault (definition).
 *
 *	@author Luca Turconi <lturconi@gmail.com>
 *  @date 22-09-2012
 *
 */ 

#include<omp.h>
//#include "/usr/lib/gcc/x86_64-unknown-linux-gnu/4.7.1/include/omp.h"
#include<limits>
#include<fstream>
#include <iomanip>
#include "fracture.hpp"
#include "cutCellProperties.hpp"

namespace Geometry
{
	// --------------------   Class Fault   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================
	Fracture::Fracture() : BilinearSurface() {}
	
	Fracture::Fracture(const Point3D & a, const Point3D & b, const Point3D & c, const Point3D & d) :
			BilinearSurface(a,b,c,d)
		{ this->setLmax();
CDARCY=(M_isMetric==true)?0.008527:0.001127;
 }
	
	Fracture::Fracture(const BilinearSurface & b) : BilinearSurface(b)
		{ this->setLmax(); 
CDARCY=(M_isMetric==true)?0.008527:0.001127;}		


	Fracture::~Fracture() {}
	
	void Fracture::setProperties(Real perm, Real compr, Real aperture)
	{	
		M_perm=perm;
		M_aperture=aperture;
		M_compr=compr;
	}

	void Fracture::setGeoProp(CProp & geoprop)
	{
		M_gridpointer=geoprop.getgridpointer();
		M_Ne=0;
		for (gmm::size_type i=0; i<geoprop.getI().size();++i){
			gmm::size_type NN= geoprop.getI()[i]-1+(geoprop.getJ()[i]-1)*geoprop.getNx()+(geoprop.getK()[i]-1)*geoprop.getNx()*geoprop.getNy()+1;
			Real aa(geoprop.getAreas()[NN]);
			if (aa>0){
			M_Ne=M_Ne+1;
			M_areas.push_back(aa);

			M_Dmedio.push_back(geoprop.getDmedio()[NN]);
			M_CG.push_back(geoprop.getCG()[NN]);
			M_S1x.push_back(geoprop.getSx(1)[i]);
			M_S1y.push_back(geoprop.getSy(1)[i]);
			M_S1z.push_back(geoprop.getSz(1)[i]);
			M_S2x.push_back(geoprop.getSx(2)[i]);
			M_S2y.push_back(geoprop.getSy(2)[i]);
			M_S2z.push_back(geoprop.getSz(2)[i]);
			M_ipos.push_back(geoprop.getI()[i]);
			M_jpos.push_back(geoprop.getJ()[i]);
			M_kpos.push_back(geoprop.getK()[i]);

			CellaF cella;

			for (gmm::size_type kk=0; kk<geoprop.getPoints(i).size();++kk){
			cella.puntiAree.push_back(geoprop.getPoints(i)[kk]);
			
			cella.raggi.push_back(geoprop.getPoints(i)[kk]-geoprop.getCG()[NN]);

			}
			cella.theta.push_back(0);
			for (gmm::size_type kk=1; kk<geoprop.getPoints(i).size();++kk){
			Real tt,rr;
			Point3D p1(cella.raggi[kk]),p2(cella.raggi[0]),p3(cella.raggi[0]);
			p3=p3.cross(this->normal(0,0));

			if (p1.norm()>0) {p1=p1/p1.norm();} 
			if (p2.norm()>0) {p2=p2/p2.norm();}
			if (p3.norm()>0) {p3=p3/p3.norm();}

			tt=p1.dot(p2);	
			rr=p1.dot(p3);	

			if (tt>=-1 && tt<=1){
			if (tt>=0 && rr>=0){
			cella.theta.push_back(acos(tt));}
			if (tt<0 && rr>=0){
			cella.theta.push_back(acos(tt));}
			if (tt>=0 && rr<0){
			cella.theta.push_back(2*3.14-acos(tt));}
			if (tt<0 && rr<0){
			cella.theta.push_back(2*3.14-acos(tt));}
			}else{ if (tt<-1) {cella.theta.push_back(3.14);}
if (tt>1) {cella.theta.push_back(0);}}
			}
		
	
			M_celle.push_back(cella);
			

			}
		}
		this->setNormaltoEdges();
		this->setHalfLength();
		this->setTransmFF();
		this->setTransmFM();
		if (M_isintby.size()>0){
		for (gmm::size_type j=0;j<M_isintby.size();++j){
			for (gmm::size_type k=0; k<geoprop.getMdmedioInt(j).size();++k){
			M_inter[j].M_dmedio.push_back(geoprop.getMdmedioInt(j)[k]); 
		}
		}		
		}
	}

	void Fracture::setNormaltoEdges()
	{
		for (gmm::size_type i=0; i<M_S1x.size();++i)
		{
			Point3D ll;
			Point3D NN;
			ll=M_S1x[i].A()-M_S1x[i].B();
			NN=ll.cross(this->normal(0.5,0.5));
			if (NN.norm()>0){M_N1x.push_back(NN/NN.norm());}
			else {M_N1x.push_back(NN);}
			ll=M_S1y[i].A()-M_S1y[i].B();
			NN=ll.cross(this->normal(0.5,0.5));
			if (NN.norm()>0){M_N1y.push_back(NN/NN.norm());}
			else {M_N1y.push_back(NN);}
			ll=M_S1z[i].A()-M_S1z[i].B();
			NN=ll.cross(this->normal(0.5,0.5));
			if (NN.norm()>0){M_N1z.push_back(NN/NN.norm());}
			else {M_N1z.push_back(NN);}
			ll=M_S2x[i].A()-M_S2x[i].B();
			NN=ll.cross(this->normal(0.5,0.5));
			if (NN.norm()>0){M_N2x.push_back(NN/NN.norm());}
			else {M_N2x.push_back(NN);}
			ll=M_S2y[i].A()-M_S2y[i].B();
			NN=ll.cross(this->normal(0.5,0.5));
			if (NN.norm()>0){M_N2y.push_back(NN/NN.norm());}
			else {M_N2y.push_back(NN);}
			ll=M_S2z[i].A()-M_S2z[i].B();
			NN=ll.cross(this->normal(0.5,0.5));
			if (NN.norm()>0){M_N2z.push_back(NN/NN.norm());}
			else {M_N2z.push_back(NN);}

		}
	}



	void Fracture::setHalfLength()
	{
		for (gmm::size_type i=0; i<M_S1x.size();++i)
		{
			Point3D ll;
			ll=2*(0.5*M_S1x[i].A()+0.5*M_S1x[i].B()-M_CG[i]);
			M_L1x.push_back(ll);
			ll=2*(0.5*M_S1y[i].A()+0.5*M_S1y[i].B()-M_CG[i]);
			M_L1y.push_back(ll);
			ll=2*(0.5*M_S1z[i].A()+0.5*M_S1z[i].B()-M_CG[i]);
			M_L1z.push_back(ll);
			ll=2*(0.5*M_S2x[i].A()+0.5*M_S2x[i].B()-M_CG[i]);
			M_L2x.push_back(ll);
			ll=2*(0.5*M_S2y[i].A()+0.5*M_S2y[i].B()-M_CG[i]);
			M_L2y.push_back(ll);
			ll=2*(0.5*M_S2z[i].A()+0.5*M_S2z[i].B()-M_CG[i]);
			M_L2z.push_back(ll);


		}
	}

	void Fracture::setHalfTransmFF()
	{
		for (gmm::size_type i=0; i<M_S1x.size();++i)
		{
			Real TT;
			TT=2*M_perm*CDARCY*M_S1x[i].length()*M_aperture*fabs(M_L1x[i].dot(M_N1x[i]))/(M_L1x[i].norm()*M_L1x[i].norm());
			M_TH1X.push_back(TT);
			TT=2*M_perm*CDARCY*M_S1y[i].length()*M_aperture*fabs(M_L1y[i].dot(M_N1y[i]))/(M_L1y[i].norm()*M_L1y[i].norm());
			M_TH1Y.push_back(TT);
			TT=2*M_perm*CDARCY*M_S1z[i].length()*M_aperture*fabs(M_L1z[i].dot(M_N1z[i]))/(M_L1z[i].norm()*M_L1z[i].norm());
			M_TH1Z.push_back(TT);
			TT=2*M_perm*CDARCY*M_S2x[i].length()*M_aperture*fabs(M_L2x[i].dot(M_N2x[i]))/(M_L2x[i].norm()*M_L2x[i].norm());
			M_TH2X.push_back(TT);
			TT=2*M_perm*CDARCY*M_S2y[i].length()*M_aperture*fabs(M_L2y[i].dot(M_N2y[i]))/(M_L2y[i].norm()*M_L2y[i].norm());
			M_TH2Y.push_back(TT);
			TT=2*M_perm*CDARCY*M_S2z[i].length()*M_aperture*fabs(M_L2z[i].dot(M_N2z[i]))/(M_L2z[i].norm()*M_L2z[i].norm());
			M_TH2Z.push_back(TT);
		}
	}

	void Fracture::setTransmFF()
	{	
		this->setHalfTransmFF();
		for (gmm::size_type i=0; i<M_S1x.size();++i)
		{ 
			Real TT(0);
			for (gmm::size_type j=0; j<M_S1x.size();++j && j!=i)
			{
				if (M_ipos[j]==M_ipos[i]+1 && M_jpos[i]==M_jpos[j] && M_kpos[i]==M_kpos[j]) {
				TT=1./(1./M_TH2X[i]+1./M_TH1X[j]);
				}
			}
		        M_TX.push_back(TT);
			TT=0;
			for (gmm::size_type j=0; j<M_S1x.size();++j && j!=i)
			{
				if (M_jpos[j]==M_jpos[i]+1 && M_ipos[i]==M_ipos[j] && M_kpos[i]==M_kpos[j]) {
				TT=1./(1./M_TH2Y[i]+1./M_TH1Y[j]);
				}
			}
		        M_TY.push_back(TT);
			TT=0;
			for (gmm::size_type j=0; j<M_S1x.size();++j && j!=i)
			{
				if (M_kpos[j]==M_kpos[i]+1 && M_jpos[i]==M_jpos[j] && M_ipos[i]==M_ipos[j]) {
				TT=1./(1./M_TH2Z[i]+1./M_TH1Z[j]);
				}
			}
		        M_TZ.push_back(TT);
		}
	}

	void Fracture::setTransmFM()
	{	
		for (gmm::size_type i=0; i<M_S1x.size();++i)
		{ 
			Real k_mX(1),k_mY(1), k_mZ(1);
			Real TT(0);
			Point3D NN(this->normal(0,0));
			if (M_Dmedio[i]>0){
			TT=CDARCY*M_areas[i]*(k_mX*NN.x*NN.x + k_mY*NN.y*NN.y + k_mZ*NN.z*NN.z)/M_Dmedio[i];
			
			}
		M_TM.push_back(TT);
		}
		
	}
	
	void Fracture::setInt(gmm::size_type i, Fracture &altra, std::vector<Point3D> punti, bool completo){
	UInt dove;
	if(completo==false){
		M_isintby.push_back(i);

}
	else{

	std::vector<gmm::size_type>::iterator it;

  // iterator to vector element:
  it = find (M_isintby.begin(), M_isintby.end(), i);
	dove=it-M_isintby.begin();

	}

		IntFrac nuova;
		
	for (gmm::size_type kk=0; kk<punti.size();++kk)
	{
		nuova.puntiInt.push_back(punti[kk]);
		
	}

	Segment SS(maxSegment(punti));
//	Point3D AA(0,0,0);
//	Segment SS(AA,AA);
	nuova.SMax=SS;
	Point3D provv(SS.A());
	provv=provv-SS.B();
	nuova.Normale=provv.cross(this->normal(0.5,0.5));

	for (gmm::size_type ii=0; ii<M_ipos.size();++ii)
		{
			for (gmm::size_type jj=0; jj<altra.M_ii().size();++jj )
			{
			   if (M_ipos[ii]==altra.M_ii()[jj] && M_jpos[ii]==altra.M_jj()[jj] &&M_kpos[ii]==altra.M_kk()[jj])
				{	
				
					CPcell CC(M_gridpointer->cell(M_ipos[ii],M_jpos[ii],M_kpos[ii]));
					std::vector<Point3D> pp;
					if (CC.segmentIntersectCell(SS,pp)){
					
					nuova.M_i.push_back(M_ipos[ii]);
					nuova.M_j.push_back(M_jpos[ii]);
					nuova.M_k.push_back(M_kpos[ii]);
					nuova.Which1.push_back(ii);
					nuova.Which2.push_back(jj);
					
					
					
					nuova.M_ax.push_back(fabs(pp[0].x-pp[1].x));
					nuova.M_ay.push_back(fabs(pp[0].y-pp[1].y));
					nuova.M_az.push_back(fabs(pp[0].z-pp[1].z));
					Real TF1F2;
					
					nuova.M_dmedio.push_back(M_inter[dove].M_dmedio[ii]);
					Point3D diff(pp[0]);
					diff=diff-pp[1];
					TF1F2=CDARCY*M_perm*diff.norm()*M_aperture/nuova.M_dmedio[ii];
					nuova.M_Tf1f2.push_back(TF1F2);

				}
	
				}
			}
		}

		if (completo==false){
		M_inter.push_back(nuova);}
		else{M_inter[dove]=nuova;

		}

	}

	bool Fracture::exportFracture(std::ofstream & myfile, gmm::size_type i)
	{
		myfile<< "FRACTURE "<<i+1 <<std::endl;
		if (!isAreaOK()) {myfile<<"ATTENZIONE: errore aree superiore al 2%"<<std::endl;}
		myfile<<"Ne  "<<M_Ne<<std::endl;
		myfile<<"N\tI\tJ\tK\tPV\t\tAREA\t\tDMEAN\t\tTMF"<<std::endl;

		for (gmm::size_type n=0;n<M_Ne;++n)
		{
			myfile<<n+1<<"\t"<<M_ipos[n]<<"\t"<<M_jpos[n]<<"\t"<<M_kpos[n]<<"\t"<<M_areas[n]*this->aperture()<<"\t\t"<<M_areas[n]<<"\t\t"
			<<M_Dmedio[n]<<"\t\t"<<M_TM[n]<<std::endl;	
		}
		
		myfile<<"Transmissibility FF"<<std::endl;
		myfile<<"N\tI\tJ\tK\tThalfX1\tThalfY1\tThalfZ1\tThalfX2\tThalfY2\tThalfZ2\tTX\tTY\tTZ"<<std::endl;

		for (gmm::size_type n=0;n<M_Ne;++n)
		{
			myfile<<n+1<<"\t"<<M_ipos[n]<<"\t"<<M_jpos[n]<<"\t"<<M_kpos[n]<<"\t"<<M_TH1X[n]<<"\t"<<M_TH1Y[n]<<"\t"
			<<M_TH1Z[n]<<"\t"<<M_TH2X[n]<<"\t"<<M_TH2Y[n]<<"\t"
			<<M_TH2Z[n]<<"\t"<<M_TX[n]<<"\t"<<M_TY[n]<<"\t"<<M_TZ[n]<<"\t"<<M_S1z[n].length()<<std::endl;	
		}
		
		myfile<< "Intersections "<<std::endl;	
		
		if (M_isintby.size()>0){
			for (gmm::size_type ff=0;ff<M_isintby.size();++ff)
			{
				myfile<<M_isintby[ff]+1<<"   "<<std::endl;
				myfile<<"MY_EL\t"<<" OTHER_EL\t"<<"i\t"<<"j\t"<<"k\t"<<  "dmedio\t"<<"Tf1f2"<<std::endl;
				for (gmm::size_type kk=0; kk<M_inter[ff].M_i.size();++kk)
				{
					myfile<<M_inter[ff].Which1[kk]+1<<"\t"<<M_inter[ff].Which2[kk]+1<<"\t"<<M_inter[ff].M_i[kk]<<
					"  "<<M_inter[ff].M_j[kk]<<"\t"<<M_inter[ff].M_k[kk]<<"\t"<<M_inter[ff].M_dmedio[kk]<<"\t"<<M_inter[ff].M_Tf1f2[kk]<<std::endl;
				}
			}
		}
		else {
			myfile<< "NONE"<<std::endl;
		}


			
	
	}

Segment maxSegment(std::vector<Point3D> & punti)
	{
		Real dist;
		Point3D diff;
		UInt i, j;
		for (gmm::size_type ii=0;ii<punti.size();++ii)
		{
			for (gmm::size_type jj=0;jj<punti.size();++jj)
			{
				diff=punti[ii]-punti[jj];
				if (dist<diff.norm()) { dist=diff.norm();
				i=ii;
				j=jj;	
				}
				
			}	

		}
		// LF: Sicuri che i e j vengano sempre assegnati???
	Segment S(punti[i],punti[j]);
	return S;	
	}

	bool Fracture::exportVtk(const std::string & fileName) const
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
		
		std::cout << std::endl << " Exporting fracture mesh in Vtk format... " << std::endl;
		
	
		// Header
		filestr << "# vtk DataFile Version 3.1" << std::endl;
		filestr << "this is a file created for Paraview" << std::endl;
		filestr << "ASCII" << std::endl;
		filestr << "DATASET POLYDATA" << std::endl;
		filestr << std::endl;	// The fifth line is empty.
		
		filestr << std::scientific << std::setprecision(10);
		gmm::size_type npoints(0);

		for (gmm::size_type kk=0; kk<M_celle.size();++kk){
		npoints=npoints+M_celle[kk].puntiAree.size();
		}
		// Pointdata
		filestr << "POINTS " << npoints << " double" << std::endl;
		for (gmm::size_type kk=0; kk<M_celle.size();++kk){

for (gmm::size_type ii=0; ii<M_celle[kk].puntiAree.size();++ii){
}
			
			std::vector<Real> tt(M_celle[kk].theta);
			for (gmm::size_type ii=0; ii<M_celle[kk].puntiAree.size();++ii){
			
			Point3D punto;
//--ciclo per ordinare i punti in senso orario
		        std::vector<Real>::iterator it;
 			it=min_element(tt.begin(),tt.end());
			punto=M_celle[kk].puntiAree[it-tt.begin()];
		
			tt[it-tt.begin()]=100.;


//-----------
		filestr <<punto.x<<"  "<< punto.y<<"  "<< punto.z<<std::endl;
		}
		}
	
		// Pointdata
		filestr << "POLYGONS " << M_celle.size() << "  "<<M_celle.size()+npoints << std::endl;
		gmm::size_type contatore(0);
		for (gmm::size_type kk=0; kk<M_celle.size();++kk){
			
		filestr << M_celle[kk].puntiAree.size();
		for (gmm::size_type ii=0; ii<M_celle[kk].puntiAree.size();++ii){
		filestr <<"  "<< contatore;
		contatore=contatore+1;
		}
		filestr<<std::endl;
		
		}
		
		filestr.close();
		
		return 1;
	}

bool Fracture::isAreaOK()   
{
	Real Areatot(this->areaFault());
	Real Areapezzi(0);
	for (gmm::size_type kk=0; kk<M_areas.size();++kk){
		Areapezzi=Areapezzi+M_areas[kk];
	}
	if (fabs(Areatot-Areapezzi)/Areatot*100>2){
	return false;	
	}
	else {return true;}
}

} // namespace Geometry
