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
			M_Sx.push_back(geoprop.getSx()[i]);
			M_Sy.push_back(geoprop.getSy()[i]);
			M_Sz.push_back(geoprop.getSz()[i]);
			M_ipos.push_back(geoprop.getI()[i]);
			M_jpos.push_back(geoprop.getJ()[i]);
			M_kpos.push_back(geoprop.getK()[i]);
			}
		}
		this->setNormaltoEdges();
		this->setHalfLength();
		this->setTransmFF();
		this->setTransmFM();
		/*if (M_isintby.size()>0){
		for (gmm::size_type j=0;j<M_isintby.size();++j){
			for (gmm::size_type k=0; k<geoprop.getMdmedioInt(j).size();++k){
			M_inter[j].M_dmedio.push_back(geoprop.getMdmedioInt(j)[k]); 
		}
		}		
		}*/
	}

	void Fracture::setNormaltoEdges()
	{
		for (gmm::size_type i=0; i<M_Sx.size();++i)
		{
			Point3D ll;
			Point3D NN;
			ll=M_Sx[i].A()-M_Sx[i].B();
			NN=ll.cross(this->normal(0.5,0.5));
			if (NN.norm()>0){M_Nx.push_back(NN/NN.norm());}
			else {M_Nx.push_back(NN);}
			ll=M_Sy[i].A()-M_Sy[i].B();
			NN=ll.cross(this->normal(0.5,0.5));
			if (NN.norm()>0){M_Ny.push_back(NN/NN.norm());}
			else {M_Ny.push_back(NN);}
			ll=M_Sz[i].A()-M_Sz[i].B();
			NN=ll.cross(this->normal(0.5,0.5));
			if (NN.norm()>0){M_Nz.push_back(NN/NN.norm());}
			else {M_Nz.push_back(NN);}

		}
	}



	void Fracture::setHalfLength()
	{
		for (gmm::size_type i=0; i<M_Sx.size();++i)
		{
			Point3D ll;
			ll=2*(M_Sx[i].A()-M_CG[i]);
			M_Lx.push_back(ll);
			ll=2*(M_Sy[i].A()-M_CG[i]);
			M_Ly.push_back(ll);
			ll=2*(M_Sz[i].A()-M_CG[i]);
			M_Lz.push_back(ll);

		}
	}

	void Fracture::setHalfTransmFF()
	{
		for (gmm::size_type i=0; i<M_Sx.size();++i)
		{
			Real TT;
			TT=2*M_perm*CDARCY*M_Sx[i].length()*M_aperture*fabs(M_Lx[i].dot(M_Nx[i]))/(M_Lx[i].norm()*M_Lx[i].norm());
			M_THX.push_back(TT);
			TT=2*M_perm*CDARCY*M_Sy[i].length()*M_aperture*fabs(M_Ly[i].dot(M_Ny[i]))/(M_Ly[i].norm()*M_Ly[i].norm());
			M_THY.push_back(TT);
			TT=2*M_perm*CDARCY*M_Sz[i].length()*M_aperture*fabs(M_Lz[i].dot(M_Nz[i]))/(M_Lz[i].norm()*M_Lz[i].norm());
			M_THZ.push_back(TT);
		}
	}

	void Fracture::setTransmFF()
	{	
		this->setHalfTransmFF();
		for (gmm::size_type i=0; i<M_Sx.size();++i)
		{ 
			Real TT(0);
			for (gmm::size_type j=0; j<M_Sx.size();++j && j!=i)
			{
				if (M_ipos[j]==M_ipos[i]+1) {
				TT=1./(1./M_THX[i]+1./M_THX[j]);
				}
			}
		        M_TX.push_back(TT);
			TT=0;
			for (gmm::size_type j=0; j<M_Sx.size();++j && j!=i)
			{
				if (M_jpos[j]==M_jpos[i]+1) {
				TT=1./(1./M_THY[i]+1./M_THY[j]);
				}
			}
		        M_TY.push_back(TT);
			TT=0;
			for (gmm::size_type j=0; j<M_Sx.size();++j && j!=i)
			{
				if (M_kpos[j]==M_kpos[i]+1) {
				TT=1./(1./M_THZ[i]+1./M_THZ[j]);
				}
			}
		        M_TZ.push_back(TT);
		}
	}

	void Fracture::setTransmFM()
	{	
		for (gmm::size_type i=0; i<M_Sx.size();++i)
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

	//Segment SS(maxSegment(punti));
	Point3D AA(0,0,0);
	Segment SS(AA,AA);
	nuova.SMax=SS;

	nuova.Normale=altra.normal(0,0) - (this->normal(0,0)).dot(altra.normal(0,0))*this->normal(0,0);

	for (gmm::size_type ii=0; ii<M_ipos.size();++ii)
		{
			for (gmm::size_type jj=0; jj<altra.M_ii().size();++jj )
			{
			   if (M_ipos[ii]==altra.M_ii()[jj] && M_jpos[ii]==altra.M_jj()[jj] &&M_kpos[ii]==altra.M_kk()[jj])
				{	
					nuova.M_i.push_back(M_ipos[ii]);
					nuova.M_j.push_back(M_jpos[ii]);
					nuova.M_k.push_back(M_kpos[ii]);
					nuova.Which1.push_back(ii);
					nuova.Which2.push_back(jj);
					
				/*	CPcell CC(M_gridpointer->cell(M_ipos[ii],M_jpos[ii],M_kpos[ii]));
					std::vector<Point3D> pp(CC.segmentIntersectCell(SS));
					
					nuova.M_ax.push_back(fabs(pp[0].x-pp[1].x));
					nuova.M_ay.push_back(fabs(pp[0].y-pp[1].y));
					nuova.M_az.push_back(fabs(pp[0].z-pp[1].z));
					Real THX,THY,THZ;
					nuova.M_dmedio=M_inter[dove].M_dmedio;

					THX=CDARCY*M_perm*fabs(pp[0].x-pp[1].x)*M_aperture/nuova.M_dmedio[ii];*/

				
	
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
		myfile<<"Ne  "<<M_Ne<<std::endl;
		myfile<<"N\tI\tJ\tK\tPV\tAREA\tDMEAN\tTMF"<<std::endl;

		for (gmm::size_type n=0;n<M_Ne;++n)
		{
			myfile<<n+1<<"\t"<<M_ipos[n]<<"\t"<<M_jpos[n]<<"\t"<<M_kpos[n]<<"\t"<<M_areas[n]*this->aperture()<<"\t"<<M_areas[n]<<"\t"
			<<M_Dmedio[n]<<"\t"<<M_TM[n]<<std::endl;	
		}
		
		myfile<<"Transmissibility FF"<<std::endl;
		myfile<<"N\tI\tJ\tK\tThalfX\tThalfY\tThalfZ\tTX\tTY\tTZ"<<std::endl;

		for (gmm::size_type n=0;n<M_Ne;++n)
		{
			myfile<<n+1<<"\t"<<M_ipos[n]<<"\t"<<M_jpos[n]<<"\t"<<M_kpos[n]<<"\t"<<M_THX[n]<<"\t"<<M_THY[n]<<"\t"
			<<M_THZ[n]<<"\t"<<M_TX[n]<<"\t"<<M_TY[n]<<"\t"<<M_TZ[n]<<std::endl;	
		}
		
		myfile<< "Intersections "<<std::endl;	
		
		if (M_isintby.size()>0){
			for (gmm::size_type ff=0;ff<M_isintby.size();++ff)
			{
				myfile<<M_isintby[ff]+1<<"   "<<std::endl;
				myfile<<"MY_EL  "<<" OTHER_EL "<<" i "<<" j "<<" k "<<std::endl;
				for (gmm::size_type kk=0; kk<M_inter[ff].M_i.size();++kk)
				{
					myfile<<M_inter[ff].Which1[kk]+1<<"       "<<M_inter[ff].Which2[kk]+1<<"         "<<M_inter[ff].M_i[kk]<<
					"  "<<M_inter[ff].M_j[kk]<<"  "<<M_inter[ff].M_k[kk]<<std::endl;
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
	Segment S(punti[i],punti[j]);
	return S;	
	}


} // namespace Geometry
