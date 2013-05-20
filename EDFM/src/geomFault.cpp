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
 
#include "geomFault.hpp"
#include "bbox.hpp"


namespace Geometry
{
	// --------------------   Class Fault   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================
	Fault::Fault() : BilinearSurface() {}
	
	Fault::Fault(const Point3D & a, const Point3D & b, const Point3D & c, const Point3D & d) :
			BilinearSurface(a,b,c,d)
		{ this->setLmax(); }
			
	Fault::Fault(const BilinearSurface & b) : BilinearSurface(b)
		{ this->setLmax(); }

	Fault::~Fault() {}
	
// ==================================================
// Methods
// ==================================================
	bool Fault::approxIsIntersectedByCell(const CPcell & c,const bool & stdDivision) const
	{
		return ( this->approxIsIntersectedBy(c.getEdge(1),stdDivision) ||
				 this->approxIsIntersectedBy(c.getEdge(2),stdDivision) ||
				 this->approxIsIntersectedBy(c.getEdge(3),stdDivision) ||
				 this->approxIsIntersectedBy(c.getEdge(4),stdDivision) ||
				 this->approxIsIntersectedBy(c.getEdge(5),stdDivision) ||
				 this->approxIsIntersectedBy(c.getEdge(6),stdDivision) ||
				 this->approxIsIntersectedBy(c.getEdge(7),stdDivision) ||
				 this->approxIsIntersectedBy(c.getEdge(8),stdDivision) ||
				 this->approxIsIntersectedBy(c.getEdge(9),stdDivision) ||
				 this->approxIsIntersectedBy(c.getEdge(10),stdDivision) ||
				 this->approxIsIntersectedBy(c.getEdge(11),stdDivision) ||
				 this->approxIsIntersectedBy(c.getEdge(12),stdDivision) );
	}
	
	void Fault::newtonIntersectionWithCell
		(const CPcell & c, Intersect::CellIntersections & cellinter,
		 const Real & toll, const UInt & maxIter) const
	{
		cellinter.i() = c.i();
		cellinter.j() = c.j();
		cellinter.k() = c.k();
		Intersect::Intersection inter;
		Point3D interPoint;
		std::vector<Real> uv_point(3,0.); 
		for(UInt i=1; i<=12; ++i)
		{
			this->newtonIntersectionWith_uv(c.getEdge(i),uv_point,toll,maxIter);
			interPoint=this->param(uv_point[0], uv_point[1]);

		if( uv_point[2]>=0)
			{
				if (uv_point[2]==0) {inter.first = i;  }
				else  {inter.first = 100+i;  }
				inter.second = interPoint;
				cellinter.insert(inter);
			}
		}
	}
	
	void Fault::approxIntersectionWithCell
		(const CPcell & c, Intersect::CellIntersections & cellinter,
		 const bool & stdDivision) const
	{
		cellinter.i() = c.i();
		cellinter.j() = c.j();
		cellinter.k() = c.k();
		Intersect::Intersection inter;
		Point3D interPoint;
		
		for(UInt i=1; i<=12; ++i)
		{
			interPoint = this->approxIntersectionWith(c.getEdge(i),stdDivision);
			
			if( interPoint.x==interPoint.x && interPoint.y==interPoint.y && interPoint.z==interPoint.z )
			{
				inter.first = i;
				inter.second = interPoint;
				cellinter.insert(inter);
			}
		}
	}
	
	void Fault::approxNewtonIntersectionWithCell
		(const CPcell & c, Intersect::CellIntersections & cellinter,
		 const Real & toll, const UInt & maxIter) const
	{
		cellinter.i() = c.i();
		cellinter.j() = c.j();
		cellinter.k() = c.k();
		Intersect::Intersection inter;
		Point3D interPoint;
		
		for(UInt i=1; i<=12; ++i)
		{
			interPoint = this->approxNewtonIntersectionWith(c.getEdge(i),toll,maxIter);
			
			if( interPoint.x==interPoint.x && interPoint.y==interPoint.y && interPoint.z==interPoint.z )
			{
				inter.first = i;
				inter.second = interPoint;
				cellinter.insert(inter);
			}
		}
	}
	
	void Fault::approxRefinedIntersectionWithCell
		(const CPcell & c, Intersect::CellIntersections & cellinter,
		 const Real & toll)
	{
		cellinter.i() = c.i();
		cellinter.j() = c.j();
		cellinter.k() = c.k();
		Intersect::Intersection inter;
		Point3D interPoint;
		
		for(UInt i=1; i<=12; ++i)
		{
			interPoint = this->approxRefinedIntersectionWith(c.getEdge(i),toll);
			
			if( interPoint.x==interPoint.x && interPoint.y==interPoint.y && interPoint.z==interPoint.z )
			{
				inter.first = i;
				inter.second = interPoint;
				cellinter.insert(inter);
			}
		}
	}

	void Fault::newtonIntersectionWithGrid_FOR3
		(const CPgrid & g, Intersect::GridIntersections & gridInter,
		 const Real & toll, const UInt & maxIter) const
	{
		#pragma omp parallel shared(g,gridInter,toll, maxIter)
		{
			Intersect::CellIntersections cellInter;
		
			#pragma omp for schedule(dynamic,1)
			for(UInt i=1; i<=g.Nx(); ++i)
			{
				for(UInt j=1; j<=g.Ny(); ++j)
				{
					for(UInt k=1; k<=g.Nz(); ++k)
					{
						if( g.cell(i,j,k).getActnum() )
						{
							this->newtonIntersectionWithCell(g.cell(i,j,k), cellInter, toll, maxIter);
							gmm::size_type counter(0);
							for (Intersect::CellIntersections_Iterator_Type cc=cellInter.begin();cc!=cellInter.end();++cc){
								if ((*cc).first<100) {counter=counter+1;}
							}
							if(cellInter.size()>0 )//&& counter>0)
							{
								#pragma omp critical
								{ gridInter.insert(cellInter); }
							}
							cellInter.clearIntersections();
						}
					}
				}
			}
		}
	}

  void Fault::newtonIntersectionWithGrid_FOR3OPT
  (const CPgrid & g, Intersect::GridIntersections & gridInter,
   const Real & toll, const UInt & maxIter) const
  {	
    ADT::ADTree const * tree(g.searchTree());
    double vmin[3];
    double vmax[3];
    // Find bounding box of the fault.
    vmin[0]=std::min(this->A().x,
		     std::min(this->B().x,
			      std::min(this->C().x,this->D().x)
			      )
		     );
    vmin[1]=std::min(this->A().y,
		     std::min(this->B().y,
			      std::min(this->C().y,this->D().y)
			      )
		     );
    vmin[2]=std::min(this->A().z,
		     std::min(this->B().z,
			      std::min(this->C().z,this->D().z)
			      )
		     );
    vmax[0]=std::max(this->A().x,
		     std::max(this->B().x,
			      std::max(this->C().x,this->D().x)
			      )
		     );
    vmax[1]=std::max(this->A().y,
		     std::max(this->B().y,
			      std::max(this->C().y,this->D().y)
			      )
		     );
    vmax[2]=std::max(this->A().z,
		     std::max(this->B().z,
			      std::max(this->C().z,this->D().z)
			      )
		     );
    
    ADT::BBox<3> box(vmin,vmax);
    std::vector<int> listFound; 
    // Use binary tree search
    tree->search(box,listFound);
#ifdef VERBOSE
    std::cout<<"Found  "<<listFound.size()<< 
      " possible intersections  "<<std::endl;
#endif
#pragma omp parallel shared(g,gridInter,toll, maxIter,tree,listFound)
    {
      Intersect::CellIntersections cellInter;
      std::vector<int> keys(3);// cambiare in unsigned
#pragma omp for schedule(dynamic,1)
      for(unsigned int it=0;
	  it<listFound.size();++it)
	{
	  keys=tree->getNode(listFound[it]).getkeys();
	  UInt i=keys[0];
	  UInt j=keys[1];
	  UInt k=keys[2];
	  if( g.cell(i,j,k).getActnum() )
	    {
	      this->newtonIntersectionWithCell(g.cell(i,j,k), cellInter, toll, maxIter);
	      if(cellInter.size()>0)
		{
#pragma omp critical
		  { gridInter.insert(cellInter); }
		}
	      cellInter.clearIntersections();
	    }
	}
    }
  }

	
	void Fault::approxIntersectionWithGrid_FOR3
		(const CPgrid & g, Intersect::GridIntersections & gridInter,
		 const bool & stdDivision) const
	{
		#pragma omp parallel shared(g,gridInter,stdDivision)
		{
			Intersect::CellIntersections cellInter;
		
			#pragma omp for schedule(dynamic,1)
			for(UInt i=1; i<=g.Nx(); ++i)
			{
				for(UInt j=1; j<=g.Ny(); ++j)
				{
					for(UInt k=1; k<=g.Nz(); ++k)
					{
						if( g.cell(i,j,k).getActnum() )
						{
							this->approxIntersectionWithCell(g.cell(i,j,k), cellInter, stdDivision);
							if(cellInter.size()>0)
							{
								#pragma omp critical
								{ gridInter.insert(cellInter); }
							}
							cellInter.clearIntersections();
						}
					}
				}
			}
		}
	}
	
	void Fault::approxNewtonIntersectionWithGrid_FOR3
		(const CPgrid & g, Intersect::GridIntersections & gridInter,
		 const Real & toll, const UInt & maxIter) const
	{
		#pragma omp parallel shared(g,gridInter, toll, maxIter)
		{
			Intersect::CellIntersections cellInter;
		
			#pragma omp for schedule(dynamic,1)
			for(UInt i=1; i<=g.Nx(); ++i)
			{
				for(UInt j=1; j<=g.Ny(); ++j)
				{
					for(UInt k=1; k<=g.Nz(); ++k)
					{
						if( g.cell(i,j,k).getActnum() )
						{
							this->approxNewtonIntersectionWithCell(g.cell(i,j,k), cellInter);
							if(cellInter.size()>0)
							{
								#pragma omp critical
								{ gridInter.insert(cellInter); }
							}
							cellInter.clearIntersections();
						}
					}
				}
			}
		}
	}
	
  /* New bisection version */
  void Fault::newtonIntersectionWithGrid_BISECTION(
						   const CPgrid & g, Intersect::GridIntersections & gridInter,
						   const Real & toll, const UInt & maxIter) const
  {	
    std::vector<UInt> BB;			
    g.buildBB(this->A(),this->B(),this->C(),this->D(),BB);
    UInt LX,RX,LY,RY,LZ,RZ;
    LX=(BB[0]>0)? BB[0]:1;
    LY=(BB[2]>0)? BB[2]:1;
    LZ=(BB[4]>0)? BB[4]:1;
    RX=(BB[1]<g.Nx()+1)? BB[1]:g.Nx();
    RY=(BB[3]<g.Ny()+1)? BB[3]:g.Ny();
    RZ=(BB[5]<g.Nz()+1)? BB[5]:g.Ny();
    std::cout << LX<<"  "<<RX<<"  "<<LY<<"  "<<RY<<"  "<<LZ <<"  "<<RZ<<std::endl;
#pragma omp parallel shared(g,gridInter,toll, maxIter)
    {
      Intersect::CellIntersections cellInter;
      
#pragma omp for schedule(dynamic,1)
      for(UInt i=LX; i<=RX; ++i)
	{
	  for(UInt j=LY; j<=RY; ++j)
	    {
	      for(UInt k=LZ; k<=RZ; ++k)
		{
		  if( g.cell(i,j,k).getActnum() )
		    {
		      this->newtonIntersectionWithCell(g.cell(i,j,k), cellInter, toll, maxIter);
		      if(cellInter.size()>0)
			{
#pragma omp critical
			  { gridInter.insert(cellInter); }
			}
		      cellInter.clearIntersections();
		    }
		}
	    }
	}
    }
  }
  
  /* end new bisection version */


	void Fault::approxRefinedIntersectionWithGrid_FOR3
		(const CPgrid & g, Intersect::GridIntersections & gridInter,
		 const Real & toll)
	{
		this->buildTriangulation(toll);
		
		#pragma omp parallel shared(g,gridInter,toll)
		{
			Intersect::CellIntersections cellInter;
		
			#pragma omp for schedule(dynamic,1)
			for(UInt i=1; i<=g.Nx(); ++i)
			{
				for(UInt j=1; j<=g.Ny(); ++j)
				{
					for(UInt k=1; k<=g.Nz(); ++k)
					{
						if( g.cell(i,j,k).getActnum() )
						{
							this->approxRefinedIntersectionWithCell(g.cell(i,j,k), cellInter, toll);
							if(cellInter.size()>0)
							{
								#pragma omp critical
								{ gridInter.insert(cellInter); }
							}
							cellInter.clearIntersections();
						}
					}
				}
			}
		}		
	}
	
	void Fault::newtonIntersectionWithGridEdge
		(const Intersect::GridEdgeMap & gridEdge,
		 Intersect::GridIntersectionMap & gridInter,
		 const Real & toll, const UInt & maxIter) const
	{
		
		//Parallel Version
		#pragma omp parallel shared(gridEdge,gridInter,toll,maxIter)
		{
			#ifdef _OPENMP
				int numThreads = omp_get_num_threads();
				int threadID = omp_get_thread_num(); // from 0 to numThreads-1
			#else
				int numThreads = 1;
				int threadID = 0;
			#endif
			
			Intersect::GridIntersectionMapElement elem;
		
			Intersect::GridEdgeMap_Const_Iterator_Type it = gridEdge.begin();
			it = Intersect::safe_advancer(it,gridEdge.end(),threadID);
			
			while( it != gridEdge.end())
			{
				elem.first = this->newtonIntersectionWith(it->first,toll,maxIter);
				if( elem.first.x==elem.first.x && elem.first.y==elem.first.y && elem.first.z==elem.first.z )
				{
					elem.second = it->second;
					#pragma omp critical
					{ gridInter.insert(elem); }
				}
				it = Intersect::safe_advancer(it,gridEdge.end(),numThreads);
			}
		}
		
	}
	
	void Fault::approxIntersectionWithGridEdge
		(const Intersect::GridEdgeMap & gridEdge,
		 Intersect::GridIntersectionMap & gridInter,
		 const bool & stdDivision) const
	{
		// Parallel Version
		#pragma omp parallel shared(gridEdge,gridInter,stdDivision)
		{
			#ifdef _OPENMP
				int numThreads = omp_get_num_threads();
				int threadID = omp_get_thread_num(); // from 0 to numThreads-1
			#else
				int numThreads = 1;
				int threadID = 0;
			#endif
			
			Intersect::GridIntersectionMapElement elem;
		
			Intersect::GridEdgeMap_Const_Iterator_Type it = gridEdge.begin();
			it = Intersect::safe_advancer(it,gridEdge.end(),threadID);
			
			while( it != gridEdge.end())
			{
				elem.first = this->approxIntersectionWith(it->first,stdDivision);
				if( elem.first.x==elem.first.x && elem.first.y==elem.first.y && elem.first.z==elem.first.z )
				{
					elem.second = it->second;
					#pragma omp critical
					{ gridInter.insert(elem); }
				}
				it = Intersect::safe_advancer(it,gridEdge.end(),numThreads);
			}
		}
		
	}
	
	void Fault::approxNewtonIntersectionWithGridEdge
		(const Intersect::GridEdgeMap & gridEdge,
		 Intersect::GridIntersectionMap & gridInter,
		 const Real & toll, const UInt & maxIter) const
	{
		// Parallel Version
		#pragma omp parallel shared(gridEdge,gridInter,toll,maxIter)
		{
			#ifdef _OPENMP
				int numThreads = omp_get_num_threads();
				int threadID = omp_get_thread_num(); // from 0 to numThreads-1
			#else
				int numThreads = 1;
				int threadID = 0;
			#endif
			
			Intersect::GridIntersectionMapElement elem;
		
			Intersect::GridEdgeMap_Const_Iterator_Type it = gridEdge.begin();
			it = Intersect::safe_advancer(it,gridEdge.end(),threadID);
			
			while( it != gridEdge.end())
			{
				elem.first = this->approxNewtonIntersectionWith(it->first,toll,maxIter);
				if( elem.first.x==elem.first.x && elem.first.y==elem.first.y && elem.first.z==elem.first.z )
				{
					elem.second = it->second;
					#pragma omp critical
					{ gridInter.insert(elem); }
				}
				it = Intersect::safe_advancer(it,gridEdge.end(),numThreads);
			}
		}
		
	}
	
	void Fault::approxRefinedIntersectionWithGridEdge
		(const Intersect::GridEdgeMap & gridEdge,
		 Intersect::GridIntersectionMap & gridInter,
		 const Real & toll)
	{
		this->buildTriangulation(toll);
		
		// Parallel Version
		#pragma omp parallel shared(gridEdge,gridInter,toll)
		{
			#ifdef _OPENMP
				int numThreads = omp_get_num_threads();
				int threadID = omp_get_thread_num(); // from 0 to numThreads-1
			#else
				int numThreads = 1;
				int threadID = 0;
			#endif
			
			Intersect::GridIntersectionMapElement elem;
		
			Intersect::GridEdgeMap_Const_Iterator_Type it = gridEdge.begin();
			it = Intersect::safe_advancer(it,gridEdge.end(),threadID);
			
			while( it != gridEdge.end())
			{
				elem.first = this->approxRefinedIntersectionWith(it->first,toll);
				if( elem.first.x==elem.first.x && elem.first.y==elem.first.y && elem.first.z==elem.first.z )
				{
					elem.second = it->second;
					#pragma omp critical
					{ gridInter.insert(elem); }
				}
				it = Intersect::safe_advancer(it,gridEdge.end(),numThreads);
			}
		}
		
	}
	
	void Fault::newtonIntersectionWithGrid_EDGEMAP
		(const CPgrid & g, Intersect::GridIntersections & gridInter,
		 const Real & toll, const UInt & maxIter) const
	{
		Intersect::GridEdgeMap gridEdge(g);
		Intersect::GridIntersectionMap gridInterMap(g);
		
		this->newtonIntersectionWithGridEdge(gridEdge,gridInterMap,toll,maxIter);
		
		gridInterMap.exportToGridIntersections(gridInter);
	}
	
	void Fault::approxIntersectionWithGrid_EDGEMAP
		(const CPgrid & g, Intersect::GridIntersections & gridInter,
		 const bool & stdDivision) const
	{
		Intersect::GridEdgeMap gridEdge(g);
		Intersect::GridIntersectionMap gridInterMap(g);
		
		this->approxIntersectionWithGridEdge(gridEdge,gridInterMap,stdDivision);
		
		gridInterMap.exportToGridIntersections(gridInter);
	}
	
	void Fault::approxNewtonIntersectionWithGrid_EDGEMAP
		(const CPgrid & g, Intersect::GridIntersections & gridInter,
		 const Real & toll, const UInt & maxIter) const
	{
		Intersect::GridEdgeMap gridEdge(g);
		Intersect::GridIntersectionMap gridInterMap(g);
		
		this->approxNewtonIntersectionWithGridEdge(gridEdge,gridInterMap,toll,maxIter);
		
		gridInterMap.exportToGridIntersections(gridInter);
	}
	
	void Fault::approxRefinedIntersectionWithGrid_EDGEMAP
		(const CPgrid & g, Intersect::GridIntersections & gridInter,
		 const Real & toll)
	{
		Intersect::GridEdgeMap gridEdge(g);
		Intersect::GridIntersectionMap gridInterMap(g);
		
		this->approxRefinedIntersectionWithGridEdge(gridEdge,gridInterMap,toll);
		
		gridInterMap.exportToGridIntersections(gridInter);
	}
	
} // namespace Geometryg.Nz()
