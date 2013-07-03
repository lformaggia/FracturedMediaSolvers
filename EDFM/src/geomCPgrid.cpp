 /*!
 *	@file geomCPgrid.cpp
 *	@brief Class for Corner Point Grid management (definition).
 *
 *	@author Luca Turconi <lturconi@gmail.com>
 *  @date 22-09-2012
 *
 */

#include <iomanip>
#include <stdexcept>
#include <memory>
#include "geomCPgrid.hpp"
#include "EclipseFileTranslator.hpp"
#include "trilinearElement.hpp" 
 
 namespace Geometry
{
	
	// --------------------   Class CPgrid   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================
    CPgrid::CPgrid(const std::string & filename, const bool & scriptWithSpecialChar, std::string direzione, Real angle,  bool rotate_z) :
		M_permDefined(0), M_tranDefined(0)
	{
		std::vector<UInt> buf;
		UInt nval;
		
		std::cout << std::endl << "Reading section SPECGRID..." << std::endl;
		if( !EclipseFile::readSection(filename,"SPECGRID",buf,scriptWithSpecialChar,3) )
			std::cerr << " *** Error: in reading section SPECGRID! *** " << std::endl;
			
		M_Nx = buf[0];
		M_Ny = buf[1];
		M_Nz = buf[2];
		
		nval = (M_Nx+1)*(M_Ny+1)*6;
		
		std::cout << std::endl << "Reading section INCLUDE..." << std::endl;

		std::string nomefile1("data/");
		
		if( readTextLine(filename,"INCLUDE", nomefile1, scriptWithSpecialChar,0) ){
			std::cout << std::endl << "Reading section COORD..." << std::endl;	
			if( !EclipseFile::readSection(nomefile1,"COORD",M_coord,scriptWithSpecialChar,nval) ){
			std::cerr << " *** Error: in reading section COORD! *** " << std::endl;
			}
		}
		else{
			std::cout << std::endl << "Reading section COORD..." << std::endl;
			if( !EclipseFile::readSection(filename,"COORD",M_coord,scriptWithSpecialChar,nval) ){
			std::cerr << " *** Error: in reading section COORD! *** " << std::endl;
			}
		}

		

		
		nval = 8*M_Nx*M_Ny*M_Nz;
		std::string nomefile2("data/");
		if( readTextLine(filename,"INCLUDE", nomefile2, scriptWithSpecialChar,1) ){
			std::cout << std::endl << "Reading section ZCORN..." << std::endl;	
			if( !EclipseFile::readSection(nomefile2,"ZCORN",M_zcorn,scriptWithSpecialChar,nval) ){
			std::cerr << " *** Error: in reading section ZCORN! *** " << std::endl;
			}
		}
		else{
			std::cout << std::endl << "Reading section ZCORN..." << std::endl;
			if( !EclipseFile::readSection(filename,"ZCORN",M_zcorn,scriptWithSpecialChar,nval) ){
			std::cerr << " *** Error: in reading section ZCORN! *** " << std::endl;
			}
		}
		
		
		nval = M_Nx*M_Ny*M_Nz;

		std::string nomefile3("data/");
		if( readTextLine(filename,"INCLUDE", nomefile3, scriptWithSpecialChar,2) ){
			std::cout << std::endl << "Reading section ACTNUM..." << std::endl;	
			if( !EclipseFile::readSection(nomefile3,"ACTNUM",M_actnum,scriptWithSpecialChar,nval) ){
			std::cerr << " *** Error: in reading section ACTNUM! *** " << std::endl;
			}
		}
		else{
			std::cout << std::endl << "Reading section ACTNUM..." << std::endl;
			if( !EclipseFile::readSection(filename,"ACTNUM",M_actnum,scriptWithSpecialChar,nval) ){
			std::cerr << " *** Error: in reading section ACTNUM! *** " << std::endl;
			}
		}
		
		CPcell cc(this->cell(1,1,M_Nz));
		CPcell cc2(this->cell(M_Nx,M_Ny,1));
		M_gridBB.push_back(cc.getVertex(1).x);
		M_gridBB.push_back(cc2.getVertex(8).x);
		M_gridBB.push_back(cc.getVertex(1).y);
		M_gridBB.push_back(cc2.getVertex(8).y);
		M_gridBB.push_back(cc.getVertex(8).z);
		M_gridBB.push_back(cc2.getVertex(1).z);

	
		M_permeability.setDim(M_Nx,M_Ny,M_Nz);
		M_permeability.setIsotropy();
		M_transmissibility.setDim(M_Nx,M_Ny,M_Nz);

		if (rotate_z)
		{	
			CPcell cc3(this->cell(1,1,M_Nz));
			Point3D p1(cc3.getVertex(1)),p2(cc3.getVertex(2));
			Real theta(-atan((p2.y-p1.y)/(p2.x-p1.x)));
					
			for (gmm::size_type ii=0;ii<M_coord.size()/3;++ii){
				gmm::size_type i(3*ii);
				Point3D ppp(M_coord[i],M_coord[i+1],M_coord[i+2]);
				ppp=ppp-p1;
				M_coord[i+2]=ppp.z + p1.z;
				M_coord[i]=ppp.x*cos(theta)-ppp.y*sin(theta) + p1.x;
				M_coord[i+1]=sin(theta)*ppp.x + cos(theta)*ppp.y + p1.y;
				Point3D centro(0.5*M_gridBB[0] + 0.5*M_gridBB[1], 0.5*M_gridBB[2] + 0.5*M_gridBB[3],0.5*M_gridBB[4] + 0.5*M_gridBB[5]);
				Point3D p(M_coord[i],M_coord[i+1],M_coord[i+2]);
				Point3D pp(apply_shear(p-centro,angle, direzione)+centro);
				M_coord[i]=pp.x;
				M_coord[i+1]=pp.y;
				M_coord[i+2]=pp.z;
			}
		}
	}

	CPgrid::~CPgrid() {}
	
// ==================================================
// Get Methods
// ==================================================
    CPpillar CPgrid::pillar(const UInt & i, const UInt & j) const
	{
		UInt l = ((i-1)+(j-1)*(M_Nx+1))*6;
		CPpillar pil(M_coord[l],M_coord[l+1],M_coord[l+2],M_coord[l+3],M_coord[l+4],M_coord[l+5],i,j);
		return pil;
	}
	
    CPcell CPgrid::cell(const UInt & i, const UInt & j, const UInt & k) const
	{
		CPpillar pilA = this->pillar(i,j);
		CPpillar pilB = this->pillar(i+1,j);
		CPpillar pilC = this->pillar(i,j+1);
		CPpillar pilD = this->pillar(i+1,j+1);
		
		std::vector<Point3D> vertex;
		
		// Vertex 1
		UInt idx1 = (i-1)*2 + (j-1)*2*(2*M_Nx) + (k-1)*2*M_Ny*(4*M_Nx);
		Real z = M_zcorn[idx1];
		vertex.push_back(pilA.getPointAtZ(z));
		// Vertex 2
		UInt idx = idx1 + 1;
		z = M_zcorn[idx];
		vertex.push_back(pilB.getPointAtZ(z));
		// Vertex 3
		idx = idx1 + 2*M_Nx;
		z = M_zcorn[idx];
		vertex.push_back(pilC.getPointAtZ(z));
		// Vertex 4
		idx += 1;
		z = M_zcorn[idx];
		vertex.push_back(pilD.getPointAtZ(z));
		// Vertex 5
		UInt idx2 = idx1 + (4*M_Nx)*M_Ny;
		z = M_zcorn[idx2];
		vertex.push_back(pilA.getPointAtZ(z));
		// Vertex 6
		idx = idx2 + 1;
		z = M_zcorn[idx];
		vertex.push_back(pilB.getPointAtZ(z));
		// Vertex 7
		idx = idx2 + 2*M_Nx;
		z = M_zcorn[idx];
		vertex.push_back(pilC.getPointAtZ(z));
		// Vertex 8
		idx += 1;
		z = M_zcorn[idx];
		vertex.push_back(pilD.getPointAtZ(z));
		
		idx = (i-1) + (j-1)*M_Nx + (k-1)*M_Nx*M_Ny; 
		UInt actnum = M_actnum[idx]; 
		
		CPcell cell(vertex,i,j,k,actnum);
		return cell;
	}
      

  //-----------------------------------------------------------------
  /*     OLD VERSION USES BISECTION 
    void CPgrid::whereIs(Point3D const & p,  std::vector<UInt> & sol) const 
  {
  InvMapResult in;
   
  UInt iL(1),iR(M_Nx),jL(1),jR(M_Ny),kL(M_Nz),kR(1);
  UInt iM(1),jM(1),kM(1);
  UInt iMOld(0),jMOld(0),kMOld(0);
  in.inside=false;

    // First iteration is always carried out!
    iM=std::floor(0.5*(iL+iR)+0.5);
    jM=std::floor(0.5*(jL+jR)+0.5);
    kM=std::floor(0.5*(kL+kR)+0.5);
    
    iM=std::max(std::min(iM,1u),M_Nx);
    jM=std::max(std::min(jM,1u),M_Ny);
    kM=std::max(std::min(kM,1u),M_Nz);
    
    do{
      
      sol[0]=iM;
      sol[1]=jM;
      sol[2]=kM;
      
      CPcell cc(cell(iM,jM,kM));
      
      in=cc.isIn2(p);
      // if is inside it is OK!!! No waste of time
      if(in.inside)break;
      
      // NON CAPISCO!
	// if (!in.inside)
	// {
	// if (iM==M_Nx-1) {iM=M_Nx;}
	// if (jM==M_Ny-1) {jM=M_Ny;}
	// if (kM==M_Nz-1) {kM=M_Nz;}
	// }
      if (in.direction[0]==-1) {iR=iM;}
      if (in.direction[0]== 1) {iL=iM;}
      if (in.direction[1]==-1) {jR=jM;}
      if (in.direction[1]== 1) {jL=jM;}
      if (in.direction[2]==-1) {kR=kM;}
      if (in.direction[2]== 1) {kL=kM;}
      
      iMOld=iM;
      jMOld=jM;
      kMOld=kM;
      
      iM=std::floor(0.5*(iL+iR)+0.5);
      jM=std::floor(0.5*(jL+jR)+0.5);
      kM=std::floor(0.5*(kL+kR)+0.5);
      
      iM=std::max(std::min(iM,1u),M_Nx);
      jM=std::max(std::min(jM,1u),M_Ny);
      kM=std::max(std::min(kM,1u),M_Nz);
    }
    while (!in.inside && (iMOld!=iM || jMOld!=jM || kMOld!=kM) );
     
    if (!in.inside) throw std::runtime_error("Cannot locate cell containing point in isIn2()");
  }
*/
  //-----------------------------------------------------------------
  // New version used search tree
  void CPgrid::whereIs(Point3D const & p,  std::vector<UInt> & sol) const 
  {
    ADT::ADTree const * tree(this->searchTree());
    volatile bool found;
    double vmin[3];
    double vmax[3];
    vmin[0]=vmax[0]=p.x;
    vmin[1]=vmax[1]=p.y;
    vmin[2]=vmax[2]=p.z;
    ADT::BBox<3> box(vmin,vmax);
    std::vector<int> listFound; 
    tree->search(box,listFound);
#ifdef VERBOSE
    std::cout<<"Found  "<<listFound.size()<< 
      " possible containing cells  "<<std::endl;
#endif
#pragma omp parallel shared(tree,listFound,sol,found)
    {
      std::vector<int> keys(3);
#pragma omp for
      for(unsigned int it=0;
	  it<listFound.size();++it)
	{
#pragma omp flush (found)
	  if(!found)
	    {
	      keys=tree->getNode(listFound[it]).getkeys();
	      UInt i=static_cast<UInt>(keys[0]);
	      UInt j=static_cast<UInt>(keys[1]);
	      UInt k=static_cast<UInt>(keys[2]);
	      CPcell cc(cell(i,j,k));
	      //bool in=cc.isIn(p);
	      InvMapResult in=cc.isIn2(p);
	      // if is inside it is OK!!! 
	      if(in.inside)
		{
#pragma omp critical
		  {
		    found=true;
#ifdef VERBOSE
		    std::cout<<in<<std::endl;
#endif
		    sol[0]=i;
		    sol[1]=j;
		    sol[2]=k;
		  }
		}
	    }
	}
    }	  
    if (!found) throw std::runtime_error("Cannot locate cell containing point in isIn2()");
  }
  
  void CPgrid::buildBB(Point3D A, Point3D B, Point3D C, Point3D D, std::vector<UInt> & BB ) const
  {
		std::vector<UInt> II;
		std::vector<UInt> JJ;
		std::vector<UInt> KK;
		std::vector<UInt> solA(3,0);
		whereIs(A,solA);
		II.push_back(solA[0]);
		JJ.push_back(solA[1]);
		KK.push_back(solA[2]);
		std::vector<UInt> solB(3,0);
		whereIs(B,solB);
		II.push_back(solB[0]);
		JJ.push_back(solB[1]); 
		KK.push_back(solB[2]);
		std::vector<UInt> solC(3,0);
		whereIs(C,solC);
		II.push_back(solC[0]); 
		JJ.push_back(solC[1]); 
		KK.push_back(solC[2]);
		std::vector<UInt> solD(3,0);
		whereIs(D,solD);
	        II.push_back(solD[0]);
		JJ.push_back(solD[1]);
		KK.push_back(solD[2]);

		std::vector<UInt>::iterator it = std::min_element (II.begin(), II.end());
		BB.push_back(*it);
	        it = std::max_element (II.begin(), II.end());
		BB.push_back(*it);
    	        it = std::min_element (JJ.begin(), JJ.end());
		BB.push_back(*it);
	        it = std::max_element (JJ.begin(), JJ.end());
		BB.push_back(*it);  
  		it = std::min_element (KK.begin(), KK.end());
		BB.push_back(*it);
	        it = std::max_element (KK.begin(), KK.end());
		BB.push_back(*it);

	}


//-----------------------------------------------------------------

	const PermeabilityField & CPgrid::getPermeabilityField() const
	{
		if(!M_permDefined)
			std::cerr << " *** Error: permeability field not defined ***" << std::endl;
		
		return M_permeability;
	}
		
	const TransmissibilityField & CPgrid::getTransmissibilityField() const
	{
		if(!M_tranDefined)
			std::cerr << " *** Error: transmissibility field not defined ***" << std::endl;
		
		return M_transmissibility;
	}

// ==================================================
// Methods
// ==================================================

	void CPgrid::rescaleZ(Real scale,Real shift){
		for (unsigned int i=0; i<M_zcorn.size(); ++i){
			M_zcorn[i]=scale*(M_zcorn[i]-shift)+shift;
		}
		
	}
	
	bool CPgrid::importPermFromFile(const std::string & filename,  const bool & iso,
									const bool & scriptWithSpecialChar)
	{
		char ans='a';
		
		while(M_permDefined && ans!='y' && ans!='n' && ans!='Y' && ans!='N'){
			std::cout << "Do you want to overwrite previous permeability field? (Y/n) : ";
			std::cin >> ans;
		}
		if(ans=='n' || ans=='N')
		{
			std::cerr << std::endl << " *** Error: Permeability field not imported! *** ";
			return 0;
		}
		
		M_permDefined=1;
		
		return M_permeability.importFromFile(filename, iso, scriptWithSpecialChar);
	}

	bool CPgrid::exportPermToFile(const std::string & filename, const std::ios_base::openmode & mode)
	{
		if(M_permDefined==0){
			std::cerr << std::endl << " *** Error: Permeability field not defined! *** ";
			return 0;
		}
		
		return M_permeability.exportToFile(filename,mode);
	}
	
	bool CPgrid::importTranFromFile(const std::string & filename,
									const bool & scriptWithSpecialChar)
	{
		char ans='a';
		
		while(M_tranDefined && ans!='y' && ans!='n' && ans!='Y' && ans!='N'){
			std::cout << "Do you want to overwrite previous transmissibility field? (Y/n) : ";
			std::cin >> ans;
		}
		if(ans=='n' || ans=='N')
		{
			std::cerr << std::endl << " *** Error: Transmissibility field not imported! *** ";
			return 0;
		}
		
		M_tranDefined=1;
		
		return M_transmissibility.importFromFile(filename, scriptWithSpecialChar);
	}
	
	bool CPgrid::exportTranToFile(const std::string & filename, const std::ios_base::openmode & mode)
	{
		if(M_tranDefined==0){
			std::cerr << std::endl << " *** Error: Transmissibility field not defined! *** ";
			return 0;
		}
		
		return M_transmissibility.exportToFile(filename,mode);
	}
		
	bool CPgrid::exportToFile(const std::string & filename, const std::ios_base::openmode & mode)
	{
		if(M_coord.size()==0)
		{
			std::cerr << std::endl << " *** Error: CP Grid not defined! *** ";
			return 0;
		}
		
		std::cout << std::endl << "Writing section SPECGRID..." << std::endl;
		std::vector<UInt> specgrid(3);
		specgrid[0] = M_Nx; specgrid[1] = M_Ny; specgrid[2] = M_Nz;
		if( !EclipseFile::writeSection(specgrid, "SPECGRID", filename, mode) )
		{
			std::cerr << std::endl << " *** Error: in writing section SPECGRID! *** "
					  << std::endl << std::endl;
			return  0;
		}
		
		std::cout << std::endl << "Writing section COORD..." << std::endl;
		if( !EclipseFile::writeSection(M_coord, "COORD", filename) )
		{
			std::cerr << std::endl << " *** Error: in writing section COORD! *** "
					  << std::endl << std::endl;
			return  0;
		}
		
		std::cout << std::endl << "Writing section ZCORN..." << std::endl;
		if( !EclipseFile::writeSection(M_zcorn, "ZCORN", filename) )
		{
			std::cerr << std::endl << " *** Error: in writing section ZCORN! *** "
					  << std::endl << std::endl;
			return  0;
		}
		
		std::cout << std::endl << "Writing section ACTNUM..." << std::endl;
		if( !EclipseFile::writeSection(M_actnum, "ACTNUM", filename) )
		{
			std::cerr << std::endl << " *** Error: in writing section ACTNUM! *** "
					  << std::endl << std::endl;
			return  0;
		}
		
		return 1;
	}
	
	bool CPgrid::exportVtk(const std::string & fileName) const
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
		
		std::cout << std::endl << " Exporting grid in Vtk format... " << std::endl;
		
		UInt nCells = M_Nx*M_Ny*M_Nz;
		UInt nPoints = 8*nCells;
		UInt CellType = 12; // for VTK_HEXAHEDRON
		UInt id1;
		Point3D p;
		
		
		// Header
		filestr << "# vtk DataFile Version 3.1" << std::endl;
		filestr << "this is a file created for Paraview" << std::endl;
		filestr << "ASCII" << std::endl;
		filestr << "DATASET UNSTRUCTURED_GRID" << std::endl;
		filestr << std::endl;	// The fifth line is empty.
		
		filestr << std::scientific << std::setprecision(10);
		
		// Pointdata
		filestr << "POINTS " << nPoints << " double" << std::endl;
		for(UInt i=1; i<=M_Nx; ++i)
		{
			for(UInt j=1; j<=M_Ny; ++j)
			{
				for(UInt k=1; k<=M_Nz; ++k)
				{
					for(UInt e=1; e<=8; ++e)
					{
						//p = this->celli,j,k).getVertex(e);
						filestr << this->cell(i,j,k).getVertex(e).x << " "
								<< this->cell(i,j,k).getVertex(e).y << " "
								<< this->cell(i,j,k).getVertex(e).z << std::endl;
					}
				}
			}
		}
		filestr << std::endl;
		
		// Celldata
		filestr << "CELLS " << nCells << " " << 9*nCells << std::endl;
		for(UInt i=1; i<=M_Nx; ++i)
		{
			for(UInt j=1; j<=M_Ny; ++j)
			{
				for(UInt k=1; k<=M_Nz; ++k)
				{
					id1 = ( (i-1)*M_Ny*M_Nz + (j-1)*M_Nz + (k-1) )*8;
						// id1 = id vertex 1 of cell (i,j,k)
					filestr << "8 " << id1+4 << " "	// id5 = id vertex 5 of cell (i,j,k)
							<< id1+5 << " "			// id6 = id vertex 6 of cell (i,j,k)
							<< id1+7 << " "			// id8 = id vertex 8 of cell (i,j,k)
							<< id1+6 << " "			// id7 = id vertex 7 of cell (i,j,k)
							<< id1 << " "			// id1 = id vertex 1 of cell (i,j,k)
							<< id1+1 << " "			// id2 = id vertex 2 of cell (i,j,k)
							<< id1+3 << " "			// id4 = id vertex 4 of cell (i,j,k)
							<< id1+2 << std::endl;	// id3 = id vertex 3 of cell (i,j,k)
				}
			}
		}
		filestr << std::endl;
		
		filestr << "CELL_TYPES " << nCells << std::endl;
			
		for(UInt i=0; i<nCells; ++i)
			filestr << CellType<< std::endl;
		filestr << std::endl;
		
		filestr.close();
		
		return 1;
	}
	
	bool CPgrid::exportVtk(const std::string & fileName, std::vector<double> & campo,const std::string & tag, int mode) const
	{
		std::fstream filestr;
		
if (mode==0){
		filestr.open (fileName.c_str(), std::ios_base::out);}
	else{
		filestr.open (fileName.c_str(), std::ios_base::app | std::ios_base::out);}

		if (filestr.is_open())
		{
			std::cout << std::endl << " File: " << fileName << ", successfully opened";
		}
		else
		{
			std::cerr << std::endl << " *** Error: file not opened *** " << std::endl << std::endl;
			return  0;
		}
		
		std::cout << std::endl << " Exporting grid in Vtk format... " << std::endl;
		
		UInt nCells = M_Nx*M_Ny*M_Nz;
		UInt nPoints = 8*nCells;
		UInt CellType = 12; // for VTK_HEXAHEDRON
		UInt id1;
		Point3D p;
		
		if (mode==0){ //primo campo che scrivo
		// Header
		filestr << "# vtk DataFile Version 3.1" << std::endl;
		filestr << "this is a file created for Paraview" << std::endl;
		filestr << "ASCII" << std::endl;
		filestr << "DATASET UNSTRUCTURED_GRID" << std::endl;
		filestr << std::endl;	// The fifth line is empty.
		
		filestr << std::scientific << std::setprecision(10);
		
		// Pointdata
		filestr << "POINTS " << nPoints << " double" << std::endl;
		for(UInt i=1; i<=M_Nx; ++i)
		{
			for(UInt j=1; j<=M_Ny; ++j)
			{
				for(UInt k=1; k<=M_Nz; ++k)
				{
					for(UInt e=1; e<=8; ++e)
					{
						//p = this->cell(i,j,k).getVertex(e);
						filestr << this->cell(i,j,k).getVertex(e).x << " "
								<< this->cell(i,j,k).getVertex(e).y << " "
								<< this->cell(i,j,k).getVertex(e).z << std::endl;
					}
				}
			}
		}
		filestr << std::endl;
		
		// Celldata
		filestr << "CELLS " << nCells << " " << 9*nCells << std::endl;
		for(UInt i=1; i<=M_Nx; ++i)
		{
			for(UInt j=1; j<=M_Ny; ++j)
			{
				for(UInt k=1; k<=M_Nz; ++k)
				{
					id1 = ( (i-1)*M_Ny*M_Nz + (j-1)*M_Nz + (k-1) )*8;
						// id1 = id vertex 1 of cell (i,j,k)
					filestr << "8 " << id1+4 << " "	// id5 = id vertex 5 of cell (i,j,k)
							<< id1+5 << " "			// id6 = id vertex 6 of cell (i,j,k)
							<< id1+7 << " "			// id8 = id vertex 8 of cell (i,j,k)
							<< id1+6 << " "			// id7 = id vertex 7 of cell (i,j,k)
							<< id1 << " "			// id1 = id vertex 1 of cell (i,j,k)
							<< id1+1 << " "			// id2 = id vertex 2 of cell (i,j,k)
							<< id1+3 << " "			// id4 = id vertex 4 of cell (i,j,k)
							<< id1+2 << std::endl;	// id3 = id vertex 3 of cell (i,j,k)
				}
			}
		}
		filestr << std::endl;
		
		filestr << "CELL_TYPES " << nCells << std::endl;
			
		for(UInt i=0; i<nCells; ++i)
			filestr << CellType<< std::endl;
		filestr << std::endl;
		
 		filestr << "CELL_DATA " << nCells <<std::endl;
filestr << "SCALARS "<< tag<<" FLOAT" <<std::endl;
    		filestr << "LOOKUP_TABLE default"<<std::endl;

for(UInt i=1; i<=M_Nx; ++i)
		{
			for(UInt j=1; j<=M_Ny; ++j)
			{
				for(UInt k=1; k<=M_Nz; ++k)
				{

filestr << campo[ i-1+(j-1)*M_Nx+(k-1)*M_Nx*M_Ny+1]<<std::endl;
				}
			}
		}
		 filestr << std::endl;}
else
{

filestr << "CELL_DATA " << nCells <<std::endl;
		filestr << "SCALARS "<< tag<<" FLOAT" <<std::endl;
    		filestr << "LOOKUP_TABLE default"<<std::endl;

for(UInt i=1; i<=M_Nx; ++i)
		{
			for(UInt j=1; j<=M_Ny; ++j)
			{
				for(UInt k=1; k<=M_Nz; ++k)
				{

filestr << campo[ i-1+(j-1)*M_Nx+(k-1)*M_Nx*M_Ny+1]<<std::endl;
				}
			}
		}
}

		filestr.close();
		
		return 1;
	}

	

	void CPgrid::showMe(std::ostream  & out) const
	{
		out << "Type = CPgrid : " << std::endl;
		out << " Dim = (Nx,Ny,Nz) : Type = (UInt,UInt,UInt) : ( "
			<< M_Nx << " , " << M_Ny << " , " << M_Nz << " )" << std::endl;
		out << " Coord : Type = vector<Real> : size = " << M_coord.size() << std::endl;
		out << " Zcorn : Type = vector<Real> : size = " << M_zcorn.size() << std::endl;
		out << " Actnum : Type = vector<UInt> : size = " << M_actnum.size() << std::endl;
		out << " PermeabilityField defined : ";

		if(M_permDefined)
		{
			out << "TRUE" << std::endl;
			M_permeability.showMe();
		}else{
			out << "FALSE" << std::endl;
		}
		out << " TransmissibilityField defined : ";
		if(M_tranDefined)
		{
			out << "TRUE" << std::endl;
			M_transmissibility.showMe();
		}else{
			out << "FALSE" << std::endl;
		}
	}

bool CPgrid::readTextLine(const std::string fileName, const std::string sectionName,
						 std::string  &data, const bool specialChar, int skip)
	{		
		std::string curLine, keyword;
		std::fstream filestr;
		
		std::cout << " Reading... " << std::endl;
		std::cout << " fileName : " << fileName << std::endl;
		std::cout << " sectionName : " << sectionName << std::endl;
		
		filestr.open (fileName.c_str(), std::fstream::in);
	
		if (filestr.is_open())
		{
			std::cout << std::endl << " File successfully opened" << std::endl;
		}
		else
		{
			std::cerr << std::endl << " *** Error: file not opened *** " << std::endl << std::endl;
			return  0;
		}

		int pos;
		int cont(-1);	
		while(!filestr.eof() )
		{
			std::getline(filestr,curLine);
		
			// remove the last character ("^M") from the string
			// otherwise it creates problems
			// ^M = character used by windows at the end of each line in text files
			// remove comments
			pos = curLine.find("--");
			curLine = curLine.substr(0,pos);
			pos=curLine.find(sectionName);
				
				if( pos!=-1 )
				{
				std::getline(filestr,curLine);
				cont+=1;
				if(specialChar)
				{curLine.erase(curLine.end()-1);}
				int found=curLine.find_first_of("/");
				curLine.resize(found);
				found=curLine.find_last_of("'");
				curLine.resize(found);
				curLine.erase(curLine.begin());
				if (cont==skip){
				data.append(curLine);
				return 1;}
				}
	
			
		}
			return 0;
		
	}

	
} // namespace Geometry
