/*!
 *	@file Triangulation2D.cpp
 *	@brief Class for unstructured triangular mesh in 2D space (definition).
 *
 *	@author Luca Turconi <luca.turconi@moxoff.com>
 *  @date 29-01-2013
 *
 */ 


#include "Triangulation2D.hpp"

#include <boost/lexical_cast.hpp>

#include<cmath>
#include<fstream>
#include<functional>


namespace Geometry{

	// --------------------   Class Fracture2D   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================
	Triangulation2D::Triangulation2D() {}
	
	Triangulation2D::Triangulation2D(const Triangulation2D & t) :
		M_cdt(t.getTriangulation()), 
		M_fn(t.getFractureNetwork()),
		M_fnConstraints(t.getFractureNetworkConstraints()),
		M_approxFractureNetwork(t.getApproxFractureNetwork()) {}
	
	Triangulation2D::Triangulation2D(const std::vector<Geometry::Point2D> & borderNodes)
	{
		std::vector<Geometry::Point2D>::const_iterator it = borderNodes.begin();
		
		// insert first node
		Vertex_handle v0 = M_cdt.insert(Point(it->x(),it->y()));
		Vertex_handle vA = v0;
		Vertex_handle vB;
		++it;
		
		while( it != borderNodes.end() )
		{
			// insert node
			vB = M_cdt.insert(Point(it->x(),it->y()));
			// insert constraint (it-1,it)
			M_cdt.insert_constraint(vA,vB);
			
			++it;
			vA = vB;
		}
		
		// insert last constraint (begin(),end()) [closure]
		M_cdt.insert_constraint(vB,v0);
	}
	
	Triangulation2D::~Triangulation2D() {}
	
// ==================================================
// Methods
// ==================================================
	std::pair<UInt,Real> Triangulation2D::maxErrorOnFactureLength() const
	{
		std::pair<UInt,Real> maxError;
		maxError.second = 0;
		Real newError;
		for( UInt f = 0; f < M_approxFractureNetwork.size(); ++f )
		{
			newError = std::fabs( approxFractureLength(M_approxFractureNetwork[f])
								  - M_fn.getNetwork()[f].length() );
								  
//			std::cout << " Fract[" << f << "] : Err = " << newError << std::endl;
//			std::cout << "    -> length() = " << M_fn.getNetwork()[f].length() << std::endl;
			if(newError > maxError.second)
			{
				maxError.first = f;
				maxError.second = newError;
			}
		}
		return maxError;
	}

	std::pair<UInt,Real> Triangulation2D::maxRelativeErrorOnFactureLength() const
	{
		std::pair<UInt,Real> maxError;
		maxError.second = 0;
		Real newError;
		for( UInt f = 0; f < M_approxFractureNetwork.size(); ++f )
		{
			newError = std::fabs( approxFractureLength(M_approxFractureNetwork[f])
								  - M_fn.getNetwork()[f].length() )
					   / M_fn.getNetwork()[f].length();
					   
//			std::cout << " Fract[" << f << "] : RelErr = " << newError << std::endl;
			if(newError > maxError.second)
			{
				maxError.first = f;
				maxError.second = newError;
			}
		}
		return maxError;
	}

	Triangulation2D::Vertex_handle Triangulation2D::addConstraint
		( const Geometry::Point2D & p )
	{
		return M_cdt.insert( Point(p.x(),p.y()) );
	}
	
	Triangulation2D::FnConstraint Triangulation2D::addConstraint
		( const Geometry::Point2D & p1, const Geometry::Point2D & p2 )
	{
		std::pair<Vertex_handle,Vertex_handle> ret;
//		ret.first = M_cdt.insert( Point(p1.x(),p1.y()) );
//		ret.second = M_cdt.insert( Point(p2.x(),p2.y()) );
		ret.first = M_cdt.insert( p1 );
		ret.second = M_cdt.insert( p2 );

		M_cdt.insert_constraint(ret.first,ret.second);
		return ret;
	}
	
	Triangulation2D::FnConstraint Triangulation2D::addConstraint
		( const Geometry::Segment2D & s )
	{
		return this->addConstraint(s.source(),s.target());
	}

	void Triangulation2D::addFractureNetwork( const Geometry::FractureNetwork2D & fn )
	{
		M_fn = fn;
		M_approxFractureNetwork.resize(fn.size());
		M_fnConstraints.resize(fn.size());
		UInt i=0;
		
		for( std::vector<Geometry::Fracture2D>::const_iterator 
			it = fn.getNetwork().begin();
		     it != fn.getNetwork().end(); ++it )
		{
			M_fnConstraints[i] = this->addConstraint( it->segment() );
			++i;
		}
	}

	void Triangulation2D::buildTriangulation( const Real alphaMin,
						  const Real maxElementEdge )
	{
		Real a(std::sin( 4*std::atan(1) * alphaMin/180 ));
		// build mesher
		Mesher mesher(M_cdt);

		// set mesher criteria
		mesher.set_criteria(Criteria(a*a,maxElementEdge));
// 0.125 is the default shape bound. It corresponds to abound 20.6 degree.
// 0.5 is the upper bound on the length of the longuest edge.

		// refine mesh
		mesher.refine_mesh();
		
		std::cout << "Number of vertices: " << M_cdt.number_of_vertices() << std::endl;
	}
	
	void Triangulation2D::buildAdaptiveTriangulation( const Real alphaMin,
					 		  const Real hmin,
					 		  const Real hmax,
					 		  const Real transitionRegion,
					 		  const Real powerLawExponent )
	{
		Real a(std::sin( 4*std::atan(1) * alphaMin/180 ));
		// build mesher
		AdaptiveMesher mesher(M_cdt);

		// set mesher criteria
		mesher.set_criteria( 
			AdaptiveCriteria(&M_fn,a*a,hmin,hmax,transitionRegion,powerLawExponent) );
// AdaptiveCriteria(shapeBound, hmin, hmax, transitionRegion, powerLawExponent)

		// refine mesh
		mesher.refine_mesh();
		
		std::cout << "Number of vertices: " << M_cdt.number_of_vertices() << std::endl;
	}

	void Triangulation2D::findApproxFractureNetwork()
	{
		for( UInt f=0; f<M_fn.size(); ++f )
		{
			M_approxFractureNetwork[f].clear();
			
			for( Vertices_in_constraint_iterator it =
				M_cdt.vertices_in_constraint_begin( M_fnConstraints[f].first,
								    M_fnConstraints[f].second );
			     it != M_cdt.vertices_in_constraint_end( M_fnConstraints[f].first,
								     M_fnConstraints[f].second );
			     ++it )
			{
				M_approxFractureNetwork[f].push_back( Vertex_handle(*it) );
			}
		}
	}
	
	bool Triangulation2D::exportApproxFractureNetworkUniqueVtk(const std::string & fileName) const
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
		
		std::cout << std::endl <<
			" Exporting Triangulation2D::ApproxFractureNetwork in Vtk format... "
			<< std::endl;
		
		UInt nPoints(0);		

		for( std::vector<ApproxFracture>::const_iterator it = M_approxFractureNetwork.begin();
			 it != M_approxFractureNetwork.end(); ++it )
		{
			nPoints += it->size();
		}
		
		UInt nCells( M_approxFractureNetwork.size() );
		
		UInt CellType = 4; // for VTK_POLY_LINE
		
		// Header
		filestr << "# vtk DataFile Version 3.1" << std::endl;
		filestr << "this is a file created for Paraview" << std::endl;
		filestr << "ASCII" << std::endl;
		filestr << "DATASET UNSTRUCTURED_GRID" << std::endl;
		filestr << std::endl;	// The fifth line is empty.
		
		filestr << std::scientific << std::setprecision(10);
		
		// Pointdata
		filestr << "POINTS " << nPoints << " double" << std::endl;

		for( std::vector<ApproxFracture>::const_iterator it = M_approxFractureNetwork.begin();
			 it != M_approxFractureNetwork.end(); ++it )
		{
			for( ApproxFracture::const_iterator jt = it->begin(); jt != it->end(); ++jt )
			{
				filestr << (*jt)->point().x() << " " << (*jt)->point().y() << " 0" << std::endl;
				filestr << std::endl;
			}	
		}
		
		// Celldata
		UInt nVal(nPoints+nCells);
		
		filestr << "CELLS " << nCells << " " << nVal << std::endl;
		
		UInt id(0);
		
		for( std::vector<ApproxFracture>::const_iterator it = M_approxFractureNetwork.begin();
			 it != M_approxFractureNetwork.end(); ++it )
		{
			filestr << it->size();
			for( ApproxFracture::const_iterator jt = it->begin(); jt != it->end(); ++jt )
			{
				filestr << " " << id;
				++id;
			}
			filestr << std::endl;
		}
		
		filestr << std::endl;
		
		filestr << "CELL_TYPES " << nCells << std::endl;
		for( std::vector<ApproxFracture>::const_iterator it = M_approxFractureNetwork.begin();
			 it != M_approxFractureNetwork.end(); ++it )
		{
			filestr << CellType << " ";
		}
		filestr << std::endl;
		
		filestr.close();
		
		return 1;
	}
	
	bool Triangulation2D::exportApproxFractureNetworkVtk(const std::string & fileName) const
	{
		bool status=1;
		
		std::cout << std::endl <<
			" Exporting Triangulation2D::ApproxFaultNetwork in Vtk format... "
			<< std::endl;
		
		for( UInt i = 0; i < M_approxFractureNetwork.size(); ++i)
			status = status && this->exportApproxFractureVtk( fileName+"_ApproxFracture"+
							  boost::lexical_cast<std::string>(i)+
							  ".vtk", i );
		
		return status;
	}
	
	bool Triangulation2D::exportApproxFractureVtk( const std::string & fileName,
												   const UInt & f) const
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
		
		std::cout << std::endl <<
			" Exporting Triangulation2D::ApproxFracture[" << f << "] in Vtk format... "
			<< std::endl;
		
		UInt nPoints( M_approxFractureNetwork[f].size() );
		
		UInt nCells(nPoints);
		
		UInt CellType = 1; // for VTK_VERTEX
		
		// Header
		filestr << "# vtk DataFile Version 3.1" << std::endl;
		filestr << "this is a file created for Paraview" << std::endl;
		filestr << "ASCII" << std::endl;
		filestr << "DATASET UNSTRUCTURED_GRID" << std::endl;
		filestr << std::endl;	// The fifth line is empty.
		
		filestr << std::scientific << std::setprecision(10);
		
		// Pointdata
		filestr << "POINTS " << nPoints << " double" << std::endl;
		for( ApproxFracture::const_iterator it = M_approxFractureNetwork[f].begin();
		     it != M_approxFractureNetwork[f].end(); ++it )
		{
			filestr << (*it)->point().x() << " " << (*it)->point().y() << " 0" << std::endl;
			filestr << std::endl;
		}
		
		// Celldata
		UInt nVal(2*nCells);
		
		filestr << "CELLS " << nCells << " " << nVal << std::endl;
		
		UInt id(0);
		
		for( ApproxFracture::const_iterator it = M_approxFractureNetwork[f].begin();
		     it != M_approxFractureNetwork[f].end(); ++it )
		{
			filestr << "1 " << id << std::endl;
			id++;
		}
		
		filestr << std::endl;
		
		filestr << "CELL_TYPES " << nCells << std::endl;
		for( ApproxFracture::const_iterator it = M_approxFractureNetwork[f].begin();
		     it != M_approxFractureNetwork[f].end(); ++it )
		{
			filestr << CellType << " ";
		}
		filestr << std::endl;
		
		filestr.close();
		
		return 1;
	}
	
	Real Triangulation2D::approxFractureLength( const ApproxFracture & appF ) const
	{
		Real len(0);
		for( ApproxFracture::const_iterator it = appF.begin()+1; it != appF.end(); ++it )
		{
			len += std::sqrt( CGAL::squared_distance( (*(it-1))->point(),
													  (*it)->point() ) );
		}
		return len;
	}
	
	
} // namespace Geometry

