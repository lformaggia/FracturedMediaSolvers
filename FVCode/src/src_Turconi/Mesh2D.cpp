 /*!
 *	@file Mesh2D.cpp
 *	@brief Class for unstructured mesh in 2D space (definitions).
 *
 *	@author Luca Turconi <luca.turconi@moxoff.com>
 *
 */ 
 
#include "Mesh2D.hpp"
#include "GeomExportVtkTool.hpp"

#include<CGAL/squared_distance_2.h>

#include<boost/lexical_cast.hpp>

#include<utility>
#include<fstream>
#include<algorithm>
#include<cmath>

  
namespace Geometry{

	// --------------------   Class Mesh2D::Edge2D   --------------------
 
// ==================================================
// Constructors & Destructor
// ==================================================
 	Mesh2D::Edge2D::Edge2D() :
 		M_mesh(), M_idP1(), M_idP2(), M_separatedCells(), M_representedFractureIds() {}
 		
 	Mesh2D::Edge2D::Edge2D(Geometry::Mesh2D * const mesh) :
 		M_mesh(mesh), M_idP1(), M_idP2(), M_separatedCells(), M_representedFractureIds() {}
 	
  	Mesh2D::Edge2D::Edge2D(const Edge2D & edge) :
  		M_mesh(edge.getMesh()), M_idP1(edge.getP1()), M_idP2(edge.getP2()),
  		M_separatedCells(edge.getSeparatedCells()),
  		M_representedFractureIds(edge.getRepresentedFractureIds()) {}
 
	Mesh2D::Edge2D::Edge2D(Geometry::Mesh2D * const mesh, const UInt & p, const UInt & q) :
		M_mesh(mesh), M_separatedCells()
	{
		if(p<q)
		{ M_idP1 = p; M_idP2 = q; }

		if(p>q)
		{ M_idP1 = q; M_idP2 = p; }
	}
	

// ==================================================
// Methods
// ==================================================
 	const Geometry::Point2D & Mesh2D::Edge2D::getPointP1() const
		{ return M_mesh->getNodesVector()[M_idP1]; }

	const Geometry::Point2D & Mesh2D::Edge2D::getPointP2() const
		{ return M_mesh->getNodesVector()[M_idP2]; }
 
 	Geometry::Segment2D Mesh2D::Edge2D::segment() const
 	{
 		return Geometry::Segment2D( M_mesh->getNodesVector()[M_idP1],
 									M_mesh->getNodesVector()[M_idP2] );
 	}

	void Mesh2D::Edge2D::setP(const UInt & p, const UInt & q)
	{
		if(p<q)
		{ M_idP1 = p; M_idP2 = q; }

		if(p>q)
		{ M_idP1 = q; M_idP2 = p; }
	}

	Real Mesh2D::Edge2D::length() const
	{
		return CGAL::sqrt( CGAL::squared_distance(this->getPointP1(),this->getPointP2()) );
	}
	
	Real Mesh2D::Edge2D::equilateralLength() const
	{
		Real A1( M_mesh->getCellsMap().at( *M_separatedCells.begin() ).area() );
		Real equiLen( 2 * CGAL::sqrt( A1 / std::sqrt(3) ) );
		
		if(M_separatedCells.size() == 2)
		{
			Real A2( M_mesh->getCellsMap().at( *(++M_separatedCells.begin()) ).area() );
			Real l2( 2 * CGAL::sqrt( A2 / std::sqrt(3) ) );
			equiLen += l2;
			equiLen /= 2;
		}
		
		return equiLen;
	}


	void Mesh2D::Edge2D::showMe(std::ostream  & out) const
	{
		out << "Type = Edge2D :";
		out << "M_mesh = " << M_mesh << std::endl;
		out << "  -> M_idP1 = " << M_idP1;
		out << " , M_idP2 = " << M_idP2 << std::endl;
		out << "  -> M_separatedCells : size = " << M_separatedCells.size();
		out << "    [ ";
		if( !M_separatedCells.empty() )
			for( std::set<UInt>::const_iterator it = M_separatedCells.begin();
				 it != M_separatedCells.end(); ++it )
			out << *it << " ";
		out << " ] " << std::endl;
		out << "  -> M_representedFractureIds : [ ";
		if( !M_representedFractureIds.empty() )
			for( std::set<UInt>::const_iterator it = M_representedFractureIds.begin();
				 it != M_representedFractureIds.end(); ++it )
			out << *it << " ";
		out << " ] " << std::endl;
	}
 
	
	// --------------------   Class Mesh2D::Cell2D   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================	
	Mesh2D::Cell2D::Cell2D() :
		M_centroid( Geometry::Point2D(0,0) ),
		M_distanceCentroidFromFractureNetwork(0), M_refinementOrder(0) {}
	
	Mesh2D::Cell2D::Cell2D(Geometry::Mesh2D * const mesh) :
		M_mesh(mesh), M_centroid( Geometry::Point2D(0,0) ),
		M_distanceCentroidFromFractureNetwork(0), M_refinementOrder(0) {}
	
	Mesh2D::Cell2D::Cell2D(const Cell2D & cell) :
		M_mesh(cell.getMesh()), M_idVertexes(cell.getVertexesVector()),
		M_idNeighbors(cell.getNeighborsSet()), M_centroid(cell.getCentroid()),
		M_distanceCentroidFromFractureNetwork(cell.getDistanceFromFractureNetwork()),
		M_refinementOrder(cell.getRefinementOrder()) {}

	Mesh2D::Cell2D::Cell2D( Geometry::Mesh2D* const mesh,
			const std::vector<UInt> & idVertexes,
			const std::set<UInt> & idNeighbors ) :
		M_mesh(mesh), M_idVertexes(idVertexes), M_idNeighbors(idNeighbors),
		M_centroid( Geometry::Point2D(0,0) ),
		M_distanceCentroidFromFractureNetwork(0), M_refinementOrder(0) {}

// ==================================================
// Methods
// ==================================================
	Geometry::Vector2D Mesh2D::Cell2D::outerNormalToEdge( const UInt & idv1,
														  const UInt & idv2 ) const
	{
		bool foundv1(0), foundv2(0);
		Geometry::Vector2D VcellSide(0,0);
		
		// control that idv1 and idv2 are vertexes of this
		for( std::vector<UInt>::const_iterator it = M_idVertexes.begin();
			 it != M_idVertexes.end(); ++it )
		{
			if( *it == idv1 )
				foundv1 = 1;
			if( *it == idv2 )
				foundv2 = 1;
		}
	
		if( !foundv1 || !foundv2 )
		{
			std::cerr << " Error: the two points are not vertexes of this cell" << std::endl;
			return VcellSide;
		}
	
		// build the vector representing the edge (idv1,idv2)
		Geometry::Vector2D Vedge( M_mesh->getNodesVector()[idv1],
								  M_mesh->getNodesVector()[idv2] );
		
		// build a vector orthogonal to (idv1,idv2)
		Geometry::Vector2D VedgePerp( Vedge.perpendicular(CGAL::CLOCKWISE) );
		
		// build a vector (idv1,vertex) pointing into the cell
		bool found(0);
		for( std::vector<UInt>::const_iterator it = M_idVertexes.begin();
			 it != M_idVertexes.end() && !found; ++it )
		{
			if( *it != idv1 && *it != idv2 )
			{
				VcellSide = Geometry::Vector2D( M_mesh->getNodesVector()[idv1],
												M_mesh->getNodesVector()[*it] );
				found = 1;
			}
		}
		
		// if VedgePerp point into the cell -> normalize and change the sign before return
		if( VedgePerp*VcellSide > 0 )
			return - ( VedgePerp / CGAL::sqrt( VedgePerp.squared_length() ) );
		
		// if VedgePerp point outside the cell -> normalize and return
		return VedgePerp / CGAL::sqrt( VedgePerp.squared_length() );
	}

	bool Mesh2D::Cell2D::hasNeighborsThroughEdge( const UInt & vertex1, const UInt & vertex2,
							   			  UInt & idNeighbor) const
	{
		bool found(0);
		
		for( std::set<UInt>::const_iterator it = M_idNeighbors.begin();
			 it != M_idNeighbors.end(); ++it )
		{
			if( std::find( M_mesh->getCellsMap().at(*it).getVertexesVector().begin(),
					   M_mesh->getCellsMap().at(*it).getVertexesVector().end(),
					   vertex1 )
				!= M_mesh->getCellsMap().at(*it).getVertexesVector().end()
				&&
				std::find( M_mesh->getCellsMap().at(*it).getVertexesVector().begin(),
					   M_mesh->getCellsMap().at(*it).getVertexesVector().end(),
					   vertex2 )
				!= M_mesh->getCellsMap().at(*it).getVertexesVector().end() )
			{
				idNeighbor = *it;
				found = 1;
			}
		}
		
		return found;
	}

 	void Mesh2D::Cell2D::computeCentroid()
 	{
		Real x(0), y(0);
 		for(std::vector<UInt>::iterator it = M_idVertexes.begin();
 			it != M_idVertexes.end(); ++it)
 		{
 			x += M_mesh->getNodesVector()[*it].x();
 			y += M_mesh->getNodesVector()[*it].y();
 		}
 		
 		M_centroid = Geometry::Point2D( x/M_idVertexes.size(),
 						y/M_idVertexes.size() );
 	}
 
 	Real Mesh2D::Cell2D::maxEdge( std::pair<UInt,UInt> & maxEdge ) const
 	{
 		Real maxLength(0);
 		Real d;
 		for( std::vector<UInt>::const_iterator it = M_idVertexes.begin()+1;
 			 it != M_idVertexes.end(); ++it )
 		{
 			d = CGAL::squared_distance( M_mesh->getNodesVector()[*(it-1)],
										M_mesh->getNodesVector()[*it] );
 			if( d > maxLength )
 			{
 				maxLength = d;
 				maxEdge.first = (it-1) - M_idVertexes.begin();
 				maxEdge.second = it - M_idVertexes.begin();
 			}
 		}
 		
 		// last edge: (end()-1,begin())
 		d = CGAL::squared_distance( M_mesh->getNodesVector()[*(M_idVertexes.end()-1)],
									M_mesh->getNodesVector()[*M_idVertexes.begin()] );
 		if( d > maxLength )
 		{
 			maxLength = d;
 			maxEdge.first = M_idVertexes.size()-1;
 			maxEdge.second = 0;
 		}
 		
 		return CGAL::sqrt(maxLength);
 	}
 
	Real Mesh2D::Cell2D::signedArea() const
 	{
 		std::vector<Geometry::Point2D> borderPoints;
 		
 		for( std::vector<UInt>::const_iterator jt = M_idVertexes.begin();
				 jt != M_idVertexes.end(); ++jt )
		{
			borderPoints.push_back( M_mesh->getNodesVector()[*jt] );
		}
 	
 		Geometry::Polygon2D polygon(borderPoints.begin(), borderPoints.end());
 		
 		return polygon.area();
 	}
 
 	Real Mesh2D::Cell2D::area() const
 	{		
 		return CGAL::abs( this->signedArea() );
 	}
 	
 	Real Mesh2D::Cell2D::incircleDiameter() const
 	{
 		if( M_idVertexes.size() != 3 )
 		{
 			std::cerr << " Error: in Cell2D::incircleDiameter" << std::endl;
 			std::cerr << "        This cell is not a triangle" << std::endl;
 			return 0;
 		}
 		
 		Real a( CGAL::sqrt(
 					CGAL::squared_distance( M_mesh->getNodesVector()[M_idVertexes[0]],
 											M_mesh->getNodesVector()[M_idVertexes[1]] ) ) );
 		Real b( CGAL::sqrt(
 					CGAL::squared_distance( M_mesh->getNodesVector()[M_idVertexes[1]],
 											M_mesh->getNodesVector()[M_idVertexes[2]] ) ) );
 		Real c( CGAL::sqrt(
 					CGAL::squared_distance( M_mesh->getNodesVector()[M_idVertexes[2]],
 											M_mesh->getNodesVector()[M_idVertexes[0]] ) ) );
 		
 		return 4 * this->area() / (a+b+c);
 	}
 	
 	Real Mesh2D::Cell2D::circumcircleDiameter() const
	{
 		if( M_idVertexes.size() != 3 )
 		{
 			std::cerr << " Error: in Cell2D::circumcircleDiameter" << std::endl;
 			std::cerr << "        This cell is not a triangle" << std::endl;
 			return 0;
 		}
 		
 		Real a( CGAL::sqrt(
 					CGAL::squared_distance( M_mesh->getNodesVector()[M_idVertexes[0]],
 											M_mesh->getNodesVector()[M_idVertexes[1]] ) ) );
 		Real b( CGAL::sqrt(
 					CGAL::squared_distance( M_mesh->getNodesVector()[M_idVertexes[1]],
 											M_mesh->getNodesVector()[M_idVertexes[2]] ) ) );
 		Real c( CGAL::sqrt(
 					CGAL::squared_distance( M_mesh->getNodesVector()[M_idVertexes[2]],
 											M_mesh->getNodesVector()[M_idVertexes[0]] ) ) );
 		
 		return a*b*c / CGAL::sqrt( (a+b+c)*(b+c-a)*(a+c-b)*(a+b-c) );
 	}
 	
 	bool Mesh2D::Cell2D::has_small_angle( const Real & minAngle ) const
 	{
 		if( M_idVertexes.size() != 3 )
		{
			std::cerr << " Error: has_small_angle works only for triangular cells!" << std::endl;
 			return 0;
		}
		
		std::pair<UInt,UInt> maxEdgeVertexId;
		Real a = this->maxEdge( maxEdgeVertexId );
		
		UInt idOpposite(0);
		for( UInt i = 0; i < 3; ++i )
			if( i != maxEdgeVertexId.first && i != maxEdgeVertexId.second )
				idOpposite = i;
		
		// a1
		Geometry::Vector2D a1(
			M_mesh->getNodesVector()[ M_idVertexes[maxEdgeVertexId.first] ] -	
			M_mesh->getNodesVector()[ M_idVertexes[maxEdgeVertexId.second] ] );
		// b
		Geometry::Vector2D b(
			M_mesh->getNodesVector()[ M_idVertexes[idOpposite] ] -
			M_mesh->getNodesVector()[ M_idVertexes[maxEdgeVertexId.second] ] );
		// a2
		Geometry::Vector2D a2(
			M_mesh->getNodesVector()[ M_idVertexes[maxEdgeVertexId.second] ] -
			M_mesh->getNodesVector()[ M_idVertexes[maxEdgeVertexId.first] ] );
		// c
		Geometry::Vector2D c(
			M_mesh->getNodesVector()[ M_idVertexes[idOpposite] ] -
			M_mesh->getNodesVector()[ M_idVertexes[maxEdgeVertexId.first] ] );

		if( ( a1*b / ( a * CGAL::sqrt( b.squared_length() ) ) ) >
				std::cos( 4*std::atan(1) * CGAL::to_double(minAngle)/180 ) ||
			( a2*c / ( a * CGAL::sqrt( c.squared_length() ) ) ) >
				std::cos( 4*std::atan(1) * CGAL::to_double(minAngle)/180 ) )
		{
				return 1;
		}
		
		return 0;
 	}
 	
	bool Mesh2D::Cell2D::has_small_angle( const Real & minAngle,
						  				  Edge2D & edgeToBeFlipped ) const
	{
		if( M_idVertexes.size() != 3 )
		{
			std::cerr << " Error: has_small_angle works only for triangular cells!" << std::endl;
 			return 0;
		}
		
		std::pair<UInt,UInt> maxEdgeVertexId;
		Real a = this->maxEdge( maxEdgeVertexId );
		
		// check that maxEdgeVertexId is not a fracture edge
		Edge2D tmpEdge( M_mesh,
						M_idVertexes[maxEdgeVertexId.first],
						M_idVertexes[maxEdgeVertexId.second]);
		std::set<Edge2D>::iterator it = M_mesh->getEdgesSet().find( tmpEdge );
		
		if( !it->isFracture() // it is not a fracture edge
			&& it->getSeparatedCells().size() == 2 // it is not a border edge
			&& ( M_mesh->getCellsMap()[*it->getSeparatedCells().begin()].vertexesNumber() == 3
				 &&
				 M_mesh->getCellsMap()
				 	[*(++it->getSeparatedCells().begin())].vertexesNumber() == 3
			   ) // it separates two triangles
		  )
		{
			UInt idOpposite(0);
			for( UInt i = 0; i < 3; ++i )
				if( i != maxEdgeVertexId.first && i != maxEdgeVertexId.second )
					idOpposite = i;
		
			// a1
			Geometry::Vector2D a1(
				M_mesh->getNodesVector()[ M_idVertexes[maxEdgeVertexId.first] ] -	
				M_mesh->getNodesVector()[ M_idVertexes[maxEdgeVertexId.second] ] );
			// b
			Geometry::Vector2D b(
				M_mesh->getNodesVector()[ M_idVertexes[idOpposite] ] -
				M_mesh->getNodesVector()[ M_idVertexes[maxEdgeVertexId.second] ] );
			// a2
			Geometry::Vector2D a2(
				M_mesh->getNodesVector()[ M_idVertexes[maxEdgeVertexId.second] ] -
				M_mesh->getNodesVector()[ M_idVertexes[maxEdgeVertexId.first] ] );
			// b
			Geometry::Vector2D c(
				M_mesh->getNodesVector()[ M_idVertexes[idOpposite] ] -
				M_mesh->getNodesVector()[ M_idVertexes[maxEdgeVertexId.first] ] );

			if( ( a1*b / ( a * CGAL::sqrt( b.squared_length() ) ) ) >
					std::cos( 4*std::atan(1) * CGAL::to_double(minAngle)/180 ) ||
				( a2*c / ( a * CGAL::sqrt( c.squared_length() ) ) ) >
					std::cos( 4*std::atan(1) * CGAL::to_double(minAngle)/180 ) )
			{
				edgeToBeFlipped = *it;
				return 1;
			}
		}
		
		return 0;	
	}
 	
 	bool Mesh2D::Cell2D::is_degenerated_triangle() const
 	{
 		if( M_idVertexes.size() != 3 )
 		{
// 			std::cerr << " Error: is_degenerated works only for triangular cells!" << std::endl;
 			return 0;
 		}
 		
 		return CGAL::collinear( M_mesh->getNodesVector()[ M_idVertexes[0] ],
 								M_mesh->getNodesVector()[ M_idVertexes[1] ],
 								M_mesh->getNodesVector()[ M_idVertexes[2] ] );
 	}
 	
 	void Mesh2D::Cell2D::showMe(std::ostream  & out) const
	{
		out << "Type = Cell2D :";
		out << " M_mesh : " << M_mesh << std::endl;
		out << "  -> M_idVertexes : size = " << M_idVertexes.size();
		out << "    [ ";
		for( std::vector<UInt>::const_iterator it = M_idVertexes.begin();
			 it != M_idVertexes.end(); ++it )
			out << *it << " ";
		out << " ] " << std::endl;
		out << "  -> Edges :" << std::endl;
		Edge2D tmpEdge( M_mesh, M_idVertexes[0], M_idVertexes[1] );
//		tmpEdge.showMe();
		M_mesh->getEdgesSet().find(tmpEdge)->showMe();
		tmpEdge.setP( M_idVertexes[1], M_idVertexes[2] );
//		tmpEdge.showMe();
		M_mesh->getEdgesSet().find(tmpEdge)->showMe();
		tmpEdge.setP( M_idVertexes[2], M_idVertexes[0] );
//		tmpEdge.showMe();
		M_mesh->getEdgesSet().find(tmpEdge)->showMe();
		
		out << " M_idNeighbors : size = " << M_idNeighbors.size();
		out << "    [ ";
		for( std::set<UInt>::const_iterator it = M_idNeighbors.begin();
			 it != M_idNeighbors.end(); ++it )
			out << *it << " ";
		out << " ] " << std::endl;
		
		out << " Scalar Properties : " << std::endl;
		out << " -> M_centroid : " << M_centroid << std::endl;
		out << " -> M_distanceCentroidFromFractureNetwork : "
			<< M_distanceCentroidFromFractureNetwork << std::endl;
		out << " -> M_refinementOrder : " << M_refinementOrder << std::endl;
	}
 	
 	
	// --------------------   Class Mesh2D::NodesConnectivity   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================
 	Mesh2D::NodesConnectivity::NodesConnectivity() {}
 	
 	Mesh2D::NodesConnectivity::NodesConnectivity(const UInt & nGridNodes) :
 		M_vector(nGridNodes,0) {}
 		
 	Mesh2D::NodesConnectivity::~NodesConnectivity()
 	{
 		this->clear();
 	}
 
// ==================================================
// Methods
// ==================================================
	void Mesh2D::NodesConnectivity::clear()
	{
		for( std::vector<ConnectionsMap *>::iterator it = M_vector.begin();
			 it != M_vector.end(); ++it )
		{
			delete (*it);
		}
		M_vector.clear();
		M_nodePatchsMap.clear();
	}
 	
	void Mesh2D::NodesConnectivity::insert( const UInt & idPoint,
											const Connection & conn )
	{
		if( M_vector[idPoint] == 0 )
		{
			M_vector[idPoint] = new ConnectionsMap;
		}
			
		if(M_vector[idPoint]->insert(conn).second == 0)
		{
			std::cerr << " ERROR: element not inserted in Mesh2D::NodesConnectivity"
					  << std::endl;
		}
	}
 	
 	void Mesh2D::NodesConnectivity::updateNodePatchsMap()
 	{
 		M_nodePatchsMap.clear();
 		for( UInt p = 0; p < M_vector.size(); ++p )
 		{
 			for( ConnectionsMap::const_iterator it = M_vector[p]->begin();
 				 it != M_vector[p]->end(); ++it )
 			{
 				for( std::set<UInt>::const_iterator jt = it->second->getSeparatedCells().begin();
 					 jt != it->second->getSeparatedCells().end(); ++jt )
 				{
 						M_nodePatchsMap[p].insert(*jt);
 				}
 			}
 		}
 	}
 
	// --------------------   Class Mesh2D   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================
 	Mesh2D::Mesh2D() {}
 	
 	Mesh2D::Mesh2D(const Geometry::Mesh2D & mesh) :
 		M_fn(mesh.M_fn), M_nodes(mesh.getNodesVector()), M_edges(mesh.getEdgesSet()),
 		M_cells(mesh.getCellsMap()), M_approxFractureNetwork(mesh.M_approxFractureNetwork) {}
 		
 	Mesh2D::Mesh2D(const Real & Lx, const Real & Ly, const UInt & Nx, const UInt & Ny)
 	{
 		Real dx(Lx/Nx);
 		Real dy(Ly/Ny);
 		
 		// inserisco i nodi
 		M_nodes.resize((Nx+1)*(Ny+1));
 		for(UInt i=0; i<=Nx; ++i)
 		{
 			for(UInt j=0; j<=Ny; ++j)
 			{
 				M_nodes[i+j*(Nx+1)] = Geometry::Point2D(i*dx,j*dy);
 			}
 		}
 		
 		//inserisco le celle
 		std::vector<UInt> idVertexes;
 		idVertexes.resize(4);
 		std::set<UInt> idNeighbors;
 		
 		std::pair<std::set<Edge2D>::iterator,bool> insertion;
 		
 		for(UInt i=0; i<Nx; ++i)
 		{
 			for(UInt j=0; j<Ny; ++j)
 			{
 				// define the vertexes of the cell
 				idVertexes[0] = i+j*(Nx+1);
 				idVertexes[1] = (i+1)+j*(Nx+1);
 				idVertexes[2] = (i+1)+(j+1)*(Nx+1);
 				idVertexes[3] = i+(j+1)*(Nx+1);

 				// insert the edges
 					// edge 0 -> 1
 				insertion = M_edges.insert( 
 					Edge2D(this,idVertexes[0],idVertexes[1]) );
 				const_cast<std::set<UInt> & >
 					(insertion.first->getSeparatedCells()).insert(i+j*Nx);

 					// edge 1 -> 2
 				insertion = M_edges.insert( 
 					Edge2D(this,idVertexes[1],idVertexes[2]) );
 				const_cast<std::set<UInt> & >
 					(insertion.first->getSeparatedCells()).insert(i+j*Nx);

 					// edge 2 -> 3
 				insertion = M_edges.insert( 
 					Edge2D(this,idVertexes[2],idVertexes[3]) );
 				const_cast<std::set<UInt> & >
 					(insertion.first->getSeparatedCells()).insert(i+j*Nx);

 					// edge 3 -> 0
 				insertion = M_edges.insert( 
 					Edge2D(this,idVertexes[3],idVertexes[0]) );
 				const_cast<std::set<UInt> & >
 					(insertion.first->getSeparatedCells()).insert(i+j*Nx);

				// define the neighbors of the cell
 				if(j>0)
 					idNeighbors.insert(i+(j-1)*Nx);
 				if(i<Nx-1)
 					idNeighbors.insert((i+1)+j*Nx);
 				if(j<Ny-1)
 					idNeighbors.insert(i+(j+1)*Nx);
 				if(i>0)
 					idNeighbors.insert((i-1)+j*Nx);
 				
 				// insert the cell
 				M_cells[i+j*Nx] = Cell2D(this,idVertexes,idNeighbors);
 				
 				idNeighbors.clear();
 			}
 		}
 	}

	Mesh2D::Mesh2D(const Geometry::Triangulation2D & t)
	{
		// copy the fracture network
		M_fn = t.getFractureNetwork();
	
		// map Vertex_handle->NodeID used to search for already inserted
		std::map<Triangulation2D::Vertex_handle,UInt> mapHandleToId;
			
		// loop on all triangulation faces (elements)
		for( Triangulation2D::CDT::Finite_faces_iterator it =
				 t.getTriangulation().finite_faces_begin();
		     it != t.getTriangulation().finite_faces_end(); ++it )
		{
			// create the new cell
			Cell2D tmpCell(this);
			UInt idCell = M_cells.size();
			
			tmpCell.getVertexesVector().resize(3);
		
			// add cell vertexes to the vector of mesh nodes
			// and fill cell M_idVertexes vector
			for( UInt i = 0; i < 3; ++i )
			{
				std::pair<std::map<Triangulation2D::Vertex_handle,UInt>::iterator,bool>
					inserted =
						mapHandleToId.insert( std::make_pair( it->vertex(i), M_nodes.size() ) );
													  
				if(inserted.second)
				{
					tmpCell.getVertexesVector()[i] = M_nodes.size();
					M_nodes.push_back((it->vertex(i))->point());
				}else{
					tmpCell.getVertexesVector()[i] = (inserted.first)->second;
				}
			}
		
			
			// add cell edges to the set of the mesh edges
			// and fill the cell M_idNeighbors set
			Edge2D tmpEdge(this);
			std::pair<std::set<Edge2D>::iterator,bool> insertion;
		
			// manage the first edge (0->1)
			tmpEdge.setP( tmpCell.getVertexesVector()[0],
				      	  tmpCell.getVertexesVector()[1] );
			tmpEdge.getSeparatedCells().insert(idCell);
		
			insertion = M_edges.insert(tmpEdge);
		
			if(!insertion.second)
			{
				// fill the cell M_idNeighbors set
					// iterator to the first and unique elements in
					// tmpEdge.M_separatedCells set
				std::set<UInt>::const_iterator it =
					insertion.first->getSeparatedCells().begin();
					
				tmpCell.getNeighborsSet().insert(*it);
				M_cells.at(*it).getNeighborsSet().insert(idCell);
				
				const_cast<std::set<UInt> & >
					(insertion.first->getSeparatedCells()).insert(idCell);
			}

			// manage the first edge (1->2)
			tmpEdge.setP( tmpCell.getVertexesVector()[1],
				      	  tmpCell.getVertexesVector()[2] );
			tmpEdge.getSeparatedCells().insert(idCell);
		
			insertion = M_edges.insert(tmpEdge);
		
			if(!insertion.second)
			{
				// fill the cell M_idNeighbors set
					// iterator to the first and unique elements in
					// tmpEdge.M_separatedCells set
				std::set<UInt>::const_iterator it =
					insertion.first->getSeparatedCells().begin();
					
				tmpCell.getNeighborsSet().insert(*it);
				M_cells.at(*it).getNeighborsSet().insert(idCell);
				
				const_cast<std::set<UInt> & >
					(insertion.first->getSeparatedCells()).insert(idCell);
			}

			// manage the final edge (2->0)
			tmpEdge.setP( tmpCell.getVertexesVector()[2],
				      	  tmpCell.getVertexesVector()[0] );
			tmpEdge.getSeparatedCells().insert(idCell);
		
			insertion = M_edges.insert(tmpEdge);
		
			if(!insertion.second)
			{
				// fill the cell M_idNeighbors set
					// iterator to the first and unique elements in
					// tmpEdge.M_separatedCells set
				std::set<UInt>::const_iterator it =
					insertion.first->getSeparatedCells().begin();
					
				tmpCell.getNeighborsSet().insert(*it);
				M_cells.at(*it).getNeighborsSet().insert(idCell);
				
				const_cast<std::set<UInt> & >
					(insertion.first->getSeparatedCells()).insert(idCell);
			}
			
			// add the cell to the mesh (at the end of the map)
			if( M_cells.size() == 0 )
				M_cells[0] = tmpCell;
			else
				M_cells[ (--M_cells.end())->first + 1 ] = tmpCell;

		}

//		// DEBUG
//		for( std::map<UInt,Cell2D>::const_iterator it = M_cells.begin();
//			 it != M_cells.end(); ++it )
//		{
//			it->second.showMe();
//		}
//		// ----------------------------------------------

// 		// DEBUG: print mapHandleToId
//		std::cout << " mapHandleToId:" << std::endl;
//		for( std::map<Triangulation2D::Vertex_handle,UInt>::const_iterator
//				 it = mapHandleToId.begin();
//			 it != mapHandleToId.end(); ++it )
//		{
//			std::cout << " pair: ( " << &(*it->first) << " , " << it->second << " )" << std::endl;
//		}

		M_approxFractureNetwork.clear();

		// mark each mesh edges with the fracture it represents
		// and fill M_approxFractureNetwork with nodes ids
		for( UInt i = 0; i < t.getApproxFractureNetwork().size(); ++i )
		{
//	std::cout << "ApproxFracture : " << i << std::endl;
			ApproxFracture appFrac;
			appFrac.push_back(
				mapHandleToId.at( *t.getApproxFractureNetwork()[i].begin() ) );

			for( Triangulation2D::ApproxFracture::const_iterator jt =
					 t.getApproxFractureNetwork()[i].begin()+1;
				 jt != t.getApproxFractureNetwork()[i].end(); ++jt )
			{
				appFrac.push_back( mapHandleToId.at(*jt) );
							
				Edge2D tmpEdge(this);
//	// DEBUG
//	std::cout << " *(jt-1) = " << &(*(*(jt-1))) << " *jt = " << &(*(*(jt)))
//			  << std::endl;
//	std::cout << " mapHandleToId.at(*(jt-1)) = " << mapHandleToId.at(*(jt-1))
//			  << " mapHandleToId.at(*(jt)) = " << mapHandleToId.at(*(jt))
//			  << std::endl;
// ----------------------------------------------
						  
				tmpEdge.setP( mapHandleToId.at(*(jt-1)), mapHandleToId.at(*jt) );

//	// DEBUG
//	std::cout << " tmpEdge.showMe()" << std::endl;
//	tmpEdge.showMe();
//	std::cout << " --------------------------" << std::endl;
// ----------------------------------------------
				
				std::set<Edge2D>::iterator found = M_edges.find(tmpEdge);
				
				if( found != M_edges.end() )
					const_cast<std::set<UInt> & >
						(found->getRepresentedFractureIds()).insert(i);
				else{
					std::cerr << " WARNING: fracture edge not found!!!" << std::endl;
//					std::vector<Geometry::Point2D> edgeVertex;
//					edgeVertex.push_back( tmpEdge.getPointP1() );
//					edgeVertex.push_back( tmpEdge.getPointP2() );
//					Geometry::exportSetOfPointsVtk( edgeVertex,
//													"./vtk/simpleTest1/non_inserted_edge.vtk" );
				}
			}
			
			M_approxFractureNetwork.push_back(appFrac);
			
		}
	}
	
	Mesh2D::~Mesh2D() {}

// ==================================================
// Methods
// ==================================================
	void Mesh2D::computeDistanceFieldFromFractureNetwork(const Geometry::FractureNetwork2D & fn)
	{
		for( std::map<UInt,Cell2D>::iterator it = M_cells.begin();
		     it != M_cells.end(); ++it )
		{
			it->second.computeCentroid();
//TODO tolto da me			it->second.computeCentroidDistanceFromFractureNetwork(fn);
		}
	}

	void Mesh2D::addDisjointCell(const std::vector<Geometry::Point2D> & cellVertexes)
	{
		// check for negative area cells
		Polygon2D cellShape(cellVertexes.begin(),cellVertexes.end());
		if( cellShape.area() < 0 )
		{
			std::cout << " *** WARNING: Mesh2D::addDisjointCell() will insert "
					  << "a negative area cell! ***"
					  << std::endl;
		}
		
		// create the new cell
		Cell2D tmpCell(this);
		UInt idCell = M_cells.size();
			
		tmpCell.getVertexesVector().resize(cellVertexes.size());
		
		// add cell vertexes to the vector of mesh nodes (without any control)
		// and fill cell M_idVertexes vector
		for( UInt i = 0; i < cellVertexes.size(); ++i )
		{
			tmpCell.getVertexesVector()[i] = M_nodes.size();
			M_nodes.push_back(cellVertexes[i]);
		}
		
		// add cell edges to the set of the mesh edges
		// and fill the cell M_idNeighbors set
		Edge2D tmpEdge(this);
		std::pair<std::set<Edge2D>::iterator,bool> insertion;
		
		for( UInt i = 0; i < cellVertexes.size()-1; ++i )
		{
			tmpEdge.setP( tmpCell.getVertexesVector()[i],
				      tmpCell.getVertexesVector()[i+1] );
			tmpEdge.getSeparatedCells().insert(idCell);
			
			insertion = M_edges.insert(tmpEdge);
			
			if(!insertion.second)
			{
				// fill the cell M_idNeighbors set
					// iterator to the first and unique elements in
					// tmpEdge.M_separatedCells set
				std::set<UInt>::const_iterator it =
					insertion.first->getSeparatedCells().begin();
				
				tmpCell.getNeighborsSet().insert(*it);
			
				const_cast<std::set<UInt> & >
					(insertion.first->getSeparatedCells()).insert(idCell);
			}
			
		}

			// manage the final edge (from last vertex to the first)
		tmpEdge.setP( tmpCell.getVertexesVector()[cellVertexes.size()-1],
			      tmpCell.getVertexesVector()[0] );
		tmpEdge.getSeparatedCells().insert(idCell);
		
		insertion = M_edges.insert(tmpEdge);
		
		if(!insertion.second)
		{
			// fill the cell M_idNeighbors set
				// iterator to the first and unique elements in
				// tmpEdge.M_separatedCells set
			std::set<UInt>::const_iterator it =
				insertion.first->getSeparatedCells().begin();
				
			tmpCell.getNeighborsSet().insert(*it);
			
			const_cast<std::set<UInt> & >
				(insertion.first->getSeparatedCells()).insert(idCell);
		}
		
		// add the cell to the mesh (at the end of the map)
		//M_cells.push_back(tmpCell);
		if( M_cells.size() == 0 )
			M_cells[0] = tmpCell;
		else
			M_cells[ (--M_cells.end())->first + 1 ] = tmpCell;

		
	}
	
	void Mesh2D::addDisjointCell(const Cell2D & cell)
	{
		std::vector<Geometry::Point2D> cellVertexes;
		for( std::vector<UInt>::const_iterator it = cell.getVertexesVector().begin();
		     it != cell.getVertexesVector().end(); ++it )
		{
			cellVertexes.push_back( cell.getMesh()->getNodesVector()[*it] );
		}
	
		this->addDisjointCell(cellVertexes);
	}
	
	void Mesh2D::addCell(const std::vector<Geometry::Point2D> & cellVertexes)
	{
		// check for negative area cells
		Polygon2D cellShape(cellVertexes.begin(),cellVertexes.end());
		if( cellShape.area() < 0 )
		{
			std::cout << " *** WARNING: Mesh2D::addDisjointCell() will insert "
					  << "a negative area cell! ***"
					  << std::endl;
		}
	
		// create the new cell
		Cell2D tmpCell(this);
		UInt idCell = M_cells.size();
			
		tmpCell.getVertexesVector().resize(cellVertexes.size());
		
		// add cell vertexes to the vector of mesh nodes (without any control)
		// and fill cell M_idVertexes vector
		for( UInt i = 0; i < cellVertexes.size(); ++i )
		{
			std::vector<Geometry::Point2D>::iterator found =
				std::find(M_nodes.begin(),M_nodes.end(),cellVertexes[i]);
		
			if(found == M_nodes.end())
			{
				tmpCell.getVertexesVector()[i] = M_nodes.size();
				M_nodes.push_back(cellVertexes[i]);
			}else{
				tmpCell.getVertexesVector()[i] = found - M_nodes.begin();
			}
		}
		
		// add cell edges to the set of the mesh edges
		// and fill the cell M_idNeighbors set
		Edge2D tmpEdge(this);
		std::pair<std::set<Edge2D>::iterator,bool> insertion;
		
		for( UInt i = 0; i < cellVertexes.size()-1; ++i )
		{
			tmpEdge.setP( tmpCell.getVertexesVector()[i],
				      	  tmpCell.getVertexesVector()[i+1] );
			tmpEdge.getSeparatedCells().insert(idCell);
			
			insertion = M_edges.insert(tmpEdge);
			
			if(!insertion.second)
			{
				// fill the cell M_idNeighbors set
					// iterator to the first and unique elements in
					// tmpEdge.M_separatedCells set
				std::set<UInt>::const_iterator it =
					insertion.first->getSeparatedCells().begin();
				
				tmpCell.getNeighborsSet().insert(*it);
			
				const_cast<std::set<UInt> & >
					(insertion.first->getSeparatedCells()).insert(idCell);
			}
			
		}

			// manage the final edge (from last vertex to the first)
		tmpEdge.setP( tmpCell.getVertexesVector()[cellVertexes.size()-1],
			      	  tmpCell.getVertexesVector()[0] );
		tmpEdge.getSeparatedCells().insert(idCell);
		
		insertion = M_edges.insert(tmpEdge);
		
		if(!insertion.second)
		{
			// fill the cell M_idNeighbors set
				// iterator to the first and unique elements in
				// tmpEdge.M_separatedCells set
			std::set<UInt>::const_iterator it =
				insertion.first->getSeparatedCells().begin();
				
			tmpCell.getNeighborsSet().insert(*it);
			
			const_cast<std::set<UInt> & >
				(insertion.first->getSeparatedCells()).insert(idCell);
		}
		
		// add the cell to the mesh (at the end of the map)
		//M_cells.push_back(tmpCell);
		if( M_cells.size() == 0 )
			M_cells[0] = tmpCell;
		else
			M_cells[ (--M_cells.end())->first + 1 ] = tmpCell;


	}
	
	void Mesh2D::addCell(const Cell2D & cell)
	{
		std::vector<Geometry::Point2D> cellVertexes;
		for( std::vector<UInt>::const_iterator it = cell.getVertexesVector().begin();
		     it != cell.getVertexesVector().end(); ++it )
		{
			cellVertexes.push_back( cell.getMesh()->getNodesVector()[*it] );
		}
	
		this->addCell(cellVertexes);
	}

	void Mesh2D::addTriangulation(const Geometry::Triangulation2D & t)
	{
		// copy the fracture network
		M_fn = t.getFractureNetwork();
	
		std::map<Triangulation2D::Vertex_handle,UInt> mapHandleToId;
		
		// loop on all triangulation faces (elements)
		for( Triangulation2D::CDT::Finite_faces_iterator it =
				 t.getTriangulation().finite_faces_begin();
		     it != t.getTriangulation().finite_faces_end(); ++it )
		{
			// create the new cell
			Cell2D tmpCell(this);
			UInt idCell = M_cells.size();
			
			tmpCell.getVertexesVector().resize(3);
		
			// add cell vertexes to the vector of mesh nodes
			// and fill cell M_idVertexes vector
			for( UInt i = 0; i < 3; ++i )
			{
				std::vector<Geometry::Point2D>::iterator found =
					std::find(M_nodes.begin(),M_nodes.end(),(it->vertex(i))->point());
		
				if(found == M_nodes.end())
				{
					tmpCell.getVertexesVector()[i] = M_nodes.size();
					M_nodes.push_back((it->vertex(i))->point());
				}else{
					tmpCell.getVertexesVector()[i] = found - M_nodes.begin();
				}
				
				mapHandleToId.insert( std::make_pair( it->vertex(i),
													  tmpCell.getVertexesVector()[i] ) );
			}
		
			
			// add cell edges to the set of the mesh edges
			// and fill the cell M_idNeighbors set
			Edge2D tmpEdge(this);
			std::pair<std::set<Edge2D>::iterator,bool> insertion;
		
			// manage the first edge (0->1)
			tmpEdge.setP( tmpCell.getVertexesVector()[0],
				      	  tmpCell.getVertexesVector()[1] );
			tmpEdge.getSeparatedCells().insert(idCell);
		
			insertion = M_edges.insert(tmpEdge);
		
			if(!insertion.second)
			{
				// fill the cell M_idNeighbors set
					// iterator to the first and unique elements in
					// tmpEdge.M_separatedCells set
				std::set<UInt>::const_iterator it =
					insertion.first->getSeparatedCells().begin();
					
				tmpCell.getNeighborsSet().insert(*it);
				M_cells.at(*it).getNeighborsSet().insert(idCell);
				
				const_cast<std::set<UInt> & >
					(insertion.first->getSeparatedCells()).insert(idCell);
			}

			// manage the first edge (1->2)
			tmpEdge.setP( tmpCell.getVertexesVector()[1],
				      	  tmpCell.getVertexesVector()[2] );
			tmpEdge.getSeparatedCells().insert(idCell);
		
			insertion = M_edges.insert(tmpEdge);
		
			if(!insertion.second)
			{
				// fill the cell M_idNeighbors set
					// iterator to the first and unique elements in
					// tmpEdge.M_separatedCells set
				std::set<UInt>::const_iterator it =
					insertion.first->getSeparatedCells().begin();
					
				tmpCell.getNeighborsSet().insert(*it);
				M_cells.at(*it).getNeighborsSet().insert(idCell);
				
				const_cast<std::set<UInt> & >
					(insertion.first->getSeparatedCells()).insert(idCell);
			}

			// manage the final edge (2->0)
			tmpEdge.setP( tmpCell.getVertexesVector()[2],
				      	  tmpCell.getVertexesVector()[0] );
			tmpEdge.getSeparatedCells().insert(idCell);
		
			insertion = M_edges.insert(tmpEdge);
		
			if(!insertion.second)
			{
				// fill the cell M_idNeighbors set
					// iterator to the first and unique elements in
					// tmpEdge.M_separatedCells set
				std::set<UInt>::const_iterator it =
					insertion.first->getSeparatedCells().begin();
					
				tmpCell.getNeighborsSet().insert(*it);
				M_cells.at(*it).getNeighborsSet().insert(idCell);
				
				const_cast<std::set<UInt> & >
					(insertion.first->getSeparatedCells()).insert(idCell);
			}
			
			// add the cell to the mesh (at the end of the map)
			if( M_cells.size() == 0 )
				M_cells[0] = tmpCell;
			else
				M_cells[ (--M_cells.end())->first + 1 ] = tmpCell;

		}

//		// DEBUG
//		for( std::map<UInt,Cell2D>::const_iterator it = M_cells.begin();
//			 it != M_cells.end(); ++it )
//		{
//			it->second.showMe();
//		}
//		// ----------------------------------------------

// 		// DEBUG: print mapHandleToId
//		std::cout << " mapHandleToId:" << std::endl;
//		for( std::map<Triangulation2D::Vertex_handle,UInt>::const_iterator
//				 it = mapHandleToId.begin();
//			 it != mapHandleToId.end(); ++it )
//		{
//			std::cout << " pair: ( " << &(*it->first) << " , " << it->second << " )" << std::endl;
//		}

		M_approxFractureNetwork.clear();

		// mark each mesh edges with the fracture it represents
		// and fill M_approxFractureNetwork with nodes ids
		for( UInt i = 0; i < t.getApproxFractureNetwork().size(); ++i )
		{
//	std::cout << "ApproxFracture : " << i << std::endl;
			ApproxFracture appFrac;
			appFrac.push_back(
				mapHandleToId.at( *t.getApproxFractureNetwork()[i].begin() ) );

			for( Triangulation2D::ApproxFracture::const_iterator jt =
					 t.getApproxFractureNetwork()[i].begin()+1;
				 jt != t.getApproxFractureNetwork()[i].end(); ++jt )
			{
				appFrac.push_back( mapHandleToId.at(*jt) );
							
				Edge2D tmpEdge(this);
//	// DEBUG
//	std::cout << " *(jt-1) = " << &(*(*(jt-1))) << " *jt = " << &(*(*(jt)))
//			  << std::endl;
//	std::cout << " mapHandleToId.at(*(jt-1)) = " << mapHandleToId.at(*(jt-1))
//			  << " mapHandleToId.at(*(jt)) = " << mapHandleToId.at(*(jt))
//			  << std::endl;
// ----------------------------------------------
						  
				tmpEdge.setP( mapHandleToId.at(*(jt-1)), mapHandleToId.at(*jt) );

//	// DEBUG
//	std::cout << " tmpEdge.showMe()" << std::endl;
//	tmpEdge.showMe();
//	std::cout << " --------------------------" << std::endl;
// ----------------------------------------------
				
				std::set<Edge2D>::iterator found = M_edges.find(tmpEdge);
				
				if( found != M_edges.end() )
					const_cast<std::set<UInt> & >
						(found->getRepresentedFractureIds()).insert(i);
				else{
					std::cerr << " WARNING: fracture edge not found!!!" << std::endl;
//					std::vector<Geometry::Point2D> edgeVertex;
//					edgeVertex.push_back( tmpEdge.getPointP1() );
//					edgeVertex.push_back( tmpEdge.getPointP2() );
//					Geometry::exportSetOfPointsVtk( edgeVertex,
//													"./vtk/simpleTest1/non_inserted_edge.vtk" );
				}
			}
			
			M_approxFractureNetwork.push_back(appFrac);
			
		}
	}

	void Mesh2D::splitMesh( Geometry::Mesh2D & inner, Geometry::Mesh2D & outer,
							const Real & cutoff) const
	{
		for( std::map<UInt,Cell2D>::const_iterator it = M_cells.begin();
		     it != M_cells.end(); ++it )
		{
			if(it->second.getDistanceFromFractureNetwork() <= cutoff)
			{
				inner.addCell(it->second);
			}else{
				outer.addCell(it->second);
			}
		}
	}

	bool Mesh2D::findCellContainingPoint( const Geometry::Point2D & point, UInt & idCell ) const
	{
		std::vector<Geometry::Point2D> borderPoints;
		
		for( std::map<UInt,Cell2D>::const_iterator it = M_cells.begin();
			 it != M_cells.end(); ++it )
		{
			borderPoints.clear();
			borderPoints.reserve( it->second.getVertexesVector().size() );
			// build a vector<Geometry::Point2D> containing the polygon vertexes
			// (counterclockwise convention)
			for( std::vector<UInt>::const_iterator jt = it->second.getVertexesVector().begin();
				 jt != it->second.getVertexesVector().end(); ++jt )
			{
				borderPoints.push_back( M_nodes[*jt] );
			}
			
			// create a Polygon_2 object
			Geometry::Polygon2D polygon(borderPoints.begin(), borderPoints.end());
			
			// control if polygon.is_simple() (precondition for bounded_side() method)
			if( !polygon.is_simple() )
			{
				std::cerr << " Error: cell " << it->first
						  << " is a non simple polygon!!! " << std::endl;
				return 0;
			}
			
			// test for position of points
			if( ! (polygon.bounded_side( point ) == CGAL::ON_UNBOUNDED_SIDE) )
			{
				idCell = it->first;
				return 1;
			}
			
		}
		
		return 0;
	}

	void Mesh2D::genericUnifyTriangle()
	{
		std::map<UInt,Cell2D>::iterator it = M_cells.begin();
		while( it != M_cells.end() )
		{
			// check for triangular cells
			if( it->second.getVertexesVector().size() == 3 )
			{
				// find longest edge in the cell
				std::pair<UInt,UInt> edgeVertex;
				Real d1( it->second.maxEdge( edgeVertex ) );
				
				// find the edge in M_edges
				Edge2D tmpEdge( this,
								it->second.getVertexesVector()[ edgeVertex.first ],
								it->second.getVertexesVector()[ edgeVertex.second ] );
				
				std::set<Edge2D>::iterator et = M_edges.find( tmpEdge );
				
//				if( !et->isFracture() && et->getSeparatedCells().size() == 2 )
				if( et->getSeparatedCells().size() == 2 )
				{
					// find id of cell2
					UInt idCell2;
					if( *et->getSeparatedCells().begin() != it->first )
					{
						idCell2 = *et->getSeparatedCells().begin();
					}else{
						idCell2 = *(++et->getSeparatedCells().begin());
					}
					
					if( M_cells.at(idCell2).getVertexesVector().size() == 3 )
					{
						// find the opposite vertex in cell1
						UInt oppositeC1;
						for( UInt i = 0; i<3; ++i )
							if( i != edgeVertex.first && i != edgeVertex.second )
								oppositeC1 = i;
					
						// find the edge vertexes and the opposite vertex in cell2
						UInt oppositeC2, firstInC2, secondInC2;
						for( UInt i = 0; i<3; ++i )
						{
							if( M_cells.at(idCell2).getVertexesVector()[i] == 
								it->second.getVertexesVector()[ edgeVertex.first ] )
							{
								firstInC2 = i;
							}else{
								if( M_cells.at(idCell2).getVertexesVector()[i] == 
									it->second.getVertexesVector()[ edgeVertex.second ] )
								{
									secondInC2 = i;
								}else{
									oppositeC2 = i;
								}
							}
						}
						
						// compute d2: dist(oppositeC1,oppositeC2)
						Real d2( CGAL::sqrt(
								   CGAL::squared_distance(
									 M_nodes[ it->second.getVertexesVector()[oppositeC1] ],
									 M_nodes[ M_cells.at(idCell2).getVertexesVector()
									 										[oppositeC2] ]
								   ) ) );
						
						// check the ratio of the two diagonals
//						if( d1/d2 > 0.8 && d1/d2 < 1.25 )

						UInt nearestOpposite;
						UInt nearestOppositeCell;
						std::set<Edge2D>::iterator ne1t, ne2t;
						
						if( et->isFracture() )
						{
							// find the nearest opposite
							if( CGAL::squared_distance(							
									M_nodes[ M_cells.at(it->first).getVertexesVector()
																				[oppositeC1] ],
									et->segment() ) <
								CGAL::squared_distance(							
									M_nodes[ M_cells.at(idCell2).getVertexesVector()
																				[oppositeC2] ],
									et->segment() )
							  )
							{
								nearestOpposite =
									M_cells.at(it->first).getVertexesVector()[oppositeC1];
								nearestOppositeCell = it->first;
							}else{
								nearestOpposite =
									M_cells.at(idCell2).getVertexesVector()[oppositeC2];
								nearestOppositeCell = idCell2;
							}
								// find the two nearest edges
							tmpEdge.setP( nearestOpposite,
									  M_cells.at(idCell2).getVertexesVector()[firstInC2] );
							ne1t = M_edges.find( tmpEdge );
													
							tmpEdge.setP( nearestOpposite,
									  M_cells.at(idCell2).getVertexesVector()[secondInC2] );
							ne2t = M_edges.find( tmpEdge );
						}	
						
						if( ( !et->isFracture() && d1/d2 > 0.8 && d1/d2 < 1.25 )
							||
							( et->isFracture() && d1/d2 > 0.8 && d1/d2 < 1.25 &&
							  M_cells.at(nearestOppositeCell).has_small_angle(11) ) )
							// considero triangoli BAD quelli con angolo minimo inferiore a:
							// 		acos(1/(1+0.02)) *180 / pi = 11 gradi
							// questa formula assicura che l'errore commesso nell'approx
							// del "lato frattura" eliminato sia inferiore a 0.05 (5%)
						{
							// unify the two cell in cell1 (it)
						
							if( et->isFracture() )
							{
								// transfer the "fracture" properties of the erased edge
								// to the two nearest edges:			
								// change the M_representedFractureIds of ne1t and ne2t
								const_cast<std::set<UInt> & >
									(ne1t->getRepresentedFractureIds()).insert(
										et->getRepresentedFractureIds().begin(),
										et->getRepresentedFractureIds().end() );
								
								const_cast<std::set<UInt> & >
									(ne2t->getRepresentedFractureIds()).insert(
										et->getRepresentedFractureIds().begin(),
										et->getRepresentedFractureIds().end() );

								// change the ApproxFracture vectors of fractures in
								// M_representedFractureIds by inserting the nearestOpposite id
								// between the two ids of the extreme points of the erased edge
								for( std::set<UInt>::const_iterator nt = 
										 et->getRepresentedFractureIds().begin();
									 nt != et->getRepresentedFractureIds().end(); ++nt )
								{
									std::vector<UInt>::iterator appFt =
										find( M_approxFractureNetwork[*nt].begin(),
											  M_approxFractureNetwork[*nt].end(),
											  et->getP1() );
									
									if( *(appFt+1) == et->getP2() )
									{
										// if P2 follows P1:
										// insert nearestOpposite before P2
										M_approxFractureNetwork[*nt]
											.insert(appFt+1,nearestOpposite);
									}else{
										// else = P2 precedes P1
										// insert nearestOpposite before P1
										M_approxFractureNetwork[*nt]
											.insert(appFt,nearestOpposite);
									}
								}
							}
							
								// erase the edge from M_edges set
							M_edges.erase( et );
							
								// add oppositeC2 in cell1 between the two vertexes of edge
							std::vector<UInt>::iterator pos =
								it->second.getVertexesVector().begin() + edgeVertex.second;
		
							it->second.getVertexesVector().insert( pos,
								M_cells.at(idCell2).getVertexesVector()[oppositeC2] );
							
								//remove cell2 from the neighbor of cell1
							it->second.getNeighborsSet().erase(idCell2);
								
								// add to c1 the neighbor of cell2 (!=cell1)
							for( std::set<UInt>::iterator kt =
									 M_cells.at(idCell2).getNeighborsSet().begin();
								 kt != M_cells.at(idCell2).getNeighborsSet().end(); ++kt )
							{	
								if( *kt != it->first )
									it->second.getNeighborsSet().insert(*kt);
							}
								
								// update the edges between cell2 and its neighbor (!=cell1)
								// and the associated neighbor cell
									// oppositeC2 -> firstInC2
										// erase cell2 from the edge M_separatedCells
							tmpEdge.setP( M_cells.at(idCell2).getVertexesVector()[oppositeC2],
										  M_cells.at(idCell2).getVertexesVector()[firstInC2] );
							et = M_edges.find( tmpEdge );
							const_cast<std::set<UInt> & >
								(et->getSeparatedCells()).erase( idCell2 );
							
										// update the neighbor cell (if exist)
							if( !et->getSeparatedCells().empty() )
							{
								UInt neighbor = *et->getSeparatedCells().begin();
								M_cells.at(neighbor).getNeighborsSet().erase(idCell2);
								M_cells.at(neighbor).getNeighborsSet().insert(it->first);
							}
										// insert cell1 in the edge M_separatedCells
							const_cast<std::set<UInt> & >
								(et->getSeparatedCells()).insert( it->first );
							
									// oppositeC2 -> secondInC2
										// erase cell2 from the edge M_separatedCells
							tmpEdge.setP( M_cells.at(idCell2).getVertexesVector()[oppositeC2],
										  M_cells.at(idCell2).getVertexesVector()[secondInC2] );
							et = M_edges.find( tmpEdge );
							const_cast<std::set<UInt> & >
								(et->getSeparatedCells()).erase( idCell2 );
								
										// update the neighbor cell (if exist)
							if( !et->getSeparatedCells().empty() )
							{
								UInt neighbor = *et->getSeparatedCells().begin();
								M_cells.at(neighbor).getNeighborsSet().erase(idCell2);
								M_cells.at(neighbor).getNeighborsSet().insert(it->first);
							}
										// insert cell1 in the edge M_separatedCells
							const_cast<std::set<UInt> & >
								(et->getSeparatedCells()).insert( it->first );
							
								// erase c2
							M_cells.erase( idCell2 );
							
							// reset the iterator on M_cells (changed after the erasement)
							it = M_cells.begin();
							
						}else{
							++it;
						}
					}else{
						++it;
					}
				}else{
					++it;
				}
			}else{	
				++it;
			}
		}
	}

	void Mesh2D::buildNodesConnectivity()
	{
		M_nodesConnectivity.clear();
		M_nodesConnectivity.setDimension( M_nodes.size() );
		
		for( std::set<Edge2D>::iterator it = M_edges.begin();
			 it != M_edges.end(); ++it )
		{
			M_nodesConnectivity.insert( it->getP1(),
				std::make_pair(it->getP2(),const_cast<Edge2D * >(&(*it))) );
			M_nodesConnectivity.insert( it->getP2(),
				std::make_pair(it->getP1(),const_cast<Edge2D * >(&(*it))) );
		}
		
		M_nodesConnectivity.updateNodePatchsMap();
	}

	Real Mesh2D::averageCellsQualityIndex() const
	{
		UInt count(0);
		Real sum(0);
		for( std::map<UInt,Cell2D>::const_iterator it = M_cells.begin();
			 it != M_cells.end(); ++it )
		{
			if( it->second.getVertexesVector().size() == 3 )
			{
				sum += it->second.qualityIndex();
				++count;
			}
		}
		
		return sum / count;
	}
	
	std::pair<UInt,Real> Mesh2D::maximumCellsQualityIndex() const
	{
		std::pair<UInt,Real> result;
		
		for( std::map<UInt,Cell2D>::const_iterator it = M_cells.begin();
			 it != M_cells.end(); ++it )
		{
			if( it->second.getVertexesVector().size() == 3 )
			{
				if( it->second.qualityIndex() > result.second )
				{
				 	result.second = it->second.qualityIndex();
					result.first = it->first;
				}
			}
		}
		
		return result;
	}

	bool Mesh2D::checkCellsConsistency( std::vector<UInt> & degeneratedCellsIds ) const
	{
		bool foundDegeneratedCells(0);
		for( std::map<UInt,Cell2D>::const_iterator it = M_cells.begin();
			 it != M_cells.end(); ++it )
		{
			if( it->second.is_degenerated_triangle() || it->second.signedArea() <= 0  )
			{
				degeneratedCellsIds.push_back( it->first );
				foundDegeneratedCells = 1;
			}
		}
	
		return foundDegeneratedCells;
	}

	bool Mesh2D::exportVtk(const std::string & fileName) const
	{
		std::fstream filestr;
		
		filestr.open (fileName.c_str(), std::ios_base::out);
	
		if (filestr.is_open())
		{
			std::cout << std::endl << " File: " << fileName << ", successfully opened";
		}
		else
		{
			std::cerr << std::endl
					  << " *** Error: file not opened *** " << std::endl << std::endl;
			return  0;
		}
		
		std::cout << std::endl << " Exporting Mesh2D in Vtk format... " << std::endl;
		
		UInt nCells = M_cells.size();
		UInt nPoints = M_nodes.size();
//		UInt CellType = 3; // for VTK_LINE
		
		// Header
		filestr << "# vtk DataFile Version 3.1" << std::endl;
		filestr << "this is a file created for Paraview" << std::endl;
		filestr << "ASCII" << std::endl;
		filestr << "DATASET UNSTRUCTURED_GRID" << std::endl;
		filestr << std::endl;	// The fifth line is empty.
		
		filestr << std::scientific << std::setprecision(10);
		
		// Pointdata
		filestr << "POINTS " << nPoints << " double" << std::endl;
		for( std::vector<Geometry::Point2D>::const_iterator it = M_nodes.begin();
		     it != M_nodes.end(); ++it )
		{
			filestr << it->x() << " " << it->y()
				<< " 0" << std::endl;
			filestr << std::endl;
		}
		
		// Celldata
		UInt nVal = 0;
		for( std::map<UInt,Cell2D>::const_iterator it = M_cells.begin();
		     it != M_cells.end(); ++it )
			nVal += it->second.vertexesNumber() +1;
		
		filestr << "CELLS " << nCells << " " << nVal << std::endl;
		
		for( std::map<UInt,Cell2D>::const_iterator it = M_cells.begin();
		     it != M_cells.end(); ++it )
		{
			filestr << it->second.vertexesNumber();
			for( std::vector<UInt>::const_iterator jt = it->second.getVertexesVector().begin();
			     jt != it->second.getVertexesVector().end(); ++jt )
			{
				filestr << " " << *jt;
			}
			filestr << std::endl;
		}
		filestr << std::endl;
		
		filestr << "CELL_TYPES " << nCells << std::endl;
		for( std::map<UInt,Cell2D>::const_iterator it = M_cells.begin();
		     it != M_cells.end(); ++it )
		{
			if(it->second.vertexesNumber() == 3)
				filestr << "5" << std::endl;
			if(it->second.vertexesNumber() == 4)
				filestr << "9" << std::endl;
			if(it->second.vertexesNumber() != 3 && it->second.vertexesNumber() != 4)
			{
				std::cerr << "   *** Error *** : Mesh2D::exportVtk() failed, ";
				std::cerr << "unexpected cell geometry";
				return 0;
			}
			filestr << std::endl;
		}
		
		filestr << "CELL_DATA " << nCells << std::endl;
		filestr << std::endl;
		filestr << "SCALARS cellID int 1" << std::endl;
		filestr << "LOOKUP_TABLE default" << std::endl;
		
		filestr << std::scientific << std::setprecision(0);
		for( std::map<UInt,Cell2D>::const_iterator it = M_cells.begin();
		     it != M_cells.end(); ++it )
		{
			filestr << it->first << " ";
		}
		filestr << std::endl;
		
		filestr.close();
		
		return 1;
	}

	bool Mesh2D::exportCentroidsDistanceFromFractureNetworkVtk(const std::string & fileName) const
	{
		std::fstream filestr;
		
		filestr.open (fileName.c_str(), std::ios_base::out);
	
		if (filestr.is_open())
		{
			std::cout << std::endl << " File: " << fileName << ", successfully opened";
		}
		else
		{
			std::cerr << std::endl
					  << " *** Error: file not opened *** " << std::endl << std::endl;
			return  0;
		}
		
		std::cout << std::endl
				  << " Exporting Mesh2D::CentroidsDistanceFromFaultNetwork in Vtk format... "
				  << std::endl;
		
		UInt nCells = M_cells.size();
		UInt nPoints = M_nodes.size();
//		UInt CellType = 3; // for VTK_LINE
		
		// Header
		filestr << "# vtk DataFile Version 3.1" << std::endl;
		filestr << "this is a file created for Paraview" << std::endl;
		filestr << "ASCII" << std::endl;
		filestr << "DATASET UNSTRUCTURED_GRID" << std::endl;
		filestr << std::endl;	// The fifth line is empty.
		
		filestr << std::scientific << std::setprecision(10);
		
		// Pointdata
		filestr << "POINTS " << nPoints << " double" << std::endl;
		for( std::vector<Geometry::Point2D>::const_iterator it = M_nodes.begin();
		     it != M_nodes.end(); ++it )
		{
			filestr << it->x() << " " << it->y() << " 0" << std::endl;
			filestr << std::endl;
		}
		
		// Celldata
		UInt nVal = 0;
		for( std::map<UInt,Cell2D>::const_iterator it = M_cells.begin();
		     it != M_cells.end(); ++it )
			nVal += it->second.vertexesNumber() +1;
		
		filestr << "CELLS " << nCells << " " << nVal << std::endl;
		
		for( std::map<UInt,Cell2D>::const_iterator it = M_cells.begin();
		     it != M_cells.end(); ++it )
		{
			filestr << it->second.vertexesNumber();
			for( std::vector<UInt>::const_iterator jt = it->second.getVertexesVector().begin();
			     jt != it->second.getVertexesVector().end(); ++jt )
			{
				filestr << " " << *jt;
			}
			filestr << std::endl;
		}
		filestr << std::endl;
		
		// Cell types
		filestr << "CELL_TYPES " << nCells << std::endl;
		for( std::map<UInt,Cell2D>::const_iterator it = M_cells.begin();
		     it != M_cells.end(); ++it )
		{
			if(it->second.vertexesNumber() == 3)
				filestr << "5" << std::endl;
			if(it->second.vertexesNumber() == 4)
				filestr << "9" << std::endl;
			if(it->second.vertexesNumber() != 3 && it->second.vertexesNumber() != 4)
			{
				std::cerr << "   *** Error *** : Mesh2D::exportVtk() failed, ";
				std::cerr << "unexpected cell geometry" << std::endl;
				return 0;
			}
			filestr << std::endl;
		}
		
		// Cell scalar field
		filestr << "CELL_DATA " << nCells << std::endl;
		filestr << std::endl;
		filestr << "SCALARS fn_distance double 1" << std::endl;
		filestr << "LOOKUP_TABLE default" << std::endl;
		for( std::map<UInt,Cell2D>::const_iterator it = M_cells.begin();
		     it != M_cells.end(); ++it )
		{
			filestr << it->second.getDistanceFromFractureNetwork() << " ";
		}
		filestr << std::endl;		
		
		filestr.close();
		
		return 1;
	} 
 
	bool Mesh2D::exportFractureEdgesVtk(const std::string & fileName) const
	{
		std::fstream filestr;
		
		filestr.open (fileName.c_str(), std::ios_base::out);
	
		if (filestr.is_open())
		{
			std::cout << std::endl << " File: " << fileName << ", successfully opened";
		}
		else
		{
			std::cerr << std::endl << " *** Error: file not opened *** "
					  << std::endl << std::endl;
			return  0;
		}
		
		std::cout << std::endl <<
			" Exporting Fracture Edges in Vtk format... "
			<< std::endl;
		
		// count the Fracture edges
		UInt nCells(0);
		for( std::set<Edge2D>::const_iterator it = M_edges.begin();
			 it != M_edges.end(); ++it )
		{
			if(it->isFracture())
				++nCells;
		}
		
		UInt nPoints(2*nCells);
		UInt CellType = 3; // for VTK_LINE
		
		// Header
		filestr << "# vtk DataFile Version 3.1" << std::endl;
		filestr << "this is a file created for Paraview" << std::endl;
		filestr << "ASCII" << std::endl;
		filestr << "DATASET UNSTRUCTURED_GRID" << std::endl;
		filestr << std::endl;	// The fifth line is empty.
		
		filestr << std::scientific << std::setprecision(10);
		
		// Pointdata
		filestr << "POINTS " << nPoints << " double" << std::endl;
		for( std::set<Edge2D>::const_iterator it = M_edges.begin();
			 it != M_edges.end(); ++it )
		{
			if(it->isFracture())
			{
				filestr << it->getPointP1().x() << " " << it->getPointP1().y()
						<< " 0" << std::endl;
				filestr << it->getPointP2().x() << " " << it->getPointP2().y()
						<< " 0" << std::endl;
			}
		}
		filestr << std::endl;
		
		// Celldata
		UInt id(0);
		filestr << "CELLS " << nCells << " " << 3*nCells << std::endl;
		for( std::set<Edge2D>::const_iterator it = M_edges.begin();
			 it != M_edges.end(); ++it )
		{
			if(it->isFracture())
			{
				filestr << "2 " << id << " ";
				++id;
				filestr << id << std::endl;
				++id;
			}
		}
		filestr << std::endl;
		
		filestr << "CELL_TYPES " << nCells << std::endl;
		for( std::set<Edge2D>::const_iterator it = M_edges.begin();
			 it != M_edges.end(); ++it )
		{
			if(it->isFracture())
			{
				filestr << CellType << std::endl;
			}
		}
		filestr << std::endl;
		
		filestr.close();
		return 1;
	}
 
	bool Mesh2D::exportCellsVtk(const std::string & fileName, const std::vector<UInt> & idCells) const
	{
		std::fstream filestr;
		
		filestr.open (fileName.c_str(), std::ios_base::out);
	
		if (filestr.is_open())
		{
			std::cout << std::endl << " File: " << fileName << ", successfully opened";
		}
		else
		{
			std::cerr << std::endl
					  << " *** Error: file not opened *** " << std::endl << std::endl;
			return  0;
		}
		
		std::cout << std::endl << " Exporting Mesh2D in Vtk format... " << std::endl;
		
		UInt nCells = idCells.size();
		
		UInt nPoints(0);
		for( std::vector<UInt>::const_iterator it = idCells.begin();
			 it != idCells.end(); ++it )
		{
			nPoints += M_cells.at(*it).vertexesNumber();
		}
		
		UInt CellType = 7; // for VTK_POLYGON
		
		// Header
		filestr << "# vtk DataFile Version 3.1" << std::endl;
		filestr << "this is a file created for Paraview" << std::endl;
		filestr << "ASCII" << std::endl;
		filestr << "DATASET UNSTRUCTURED_GRID" << std::endl;
		filestr << std::endl;	// The fifth line is empty.
		
		filestr << std::scientific << std::setprecision(10);
		
		// Pointdata
		filestr << "POINTS " << nPoints << " double" << std::endl;
		for( std::vector<UInt>::const_iterator it = idCells.begin();
		     it != idCells.end(); ++it )
		{
			for( std::vector<UInt>::const_iterator jt =
					 M_cells.at(*it).getVertexesVector().begin();
				 jt != M_cells.at(*it).getVertexesVector().end(); ++jt )
			{
				filestr << M_nodes[*jt].x() << " " << M_nodes[*jt].y()
						<< " 0" << std::endl;
				filestr << std::endl;
			
			}
		}
		
		// Celldata
		UInt nVal = nPoints + nCells;
		
		filestr << "CELLS " << nCells << " " << nVal << std::endl;
		
		UInt idPoint(0);
		for( std::vector<UInt>::const_iterator it = idCells.begin();
		     it != idCells.end(); ++it )
		{
			filestr << M_cells.at(*it).vertexesNumber();
			for( std::vector<UInt>::const_iterator jt =
					 M_cells.at(*it).getVertexesVector().begin();
			     jt != M_cells.at(*it).getVertexesVector().end(); ++jt )
			{
				filestr << " " << idPoint++;
			}
			filestr << std::endl;
		}
		filestr << std::endl;
		
		filestr << "CELL_TYPES " << nCells << std::endl;
		for( std::vector<UInt>::const_iterator it = idCells.begin();
		     it != idCells.end(); ++it )
		{
			filestr << CellType << std::endl;
		}
		
		filestr << "CELL_DATA " << nCells << std::endl;
		filestr << std::endl;
		filestr << "SCALARS cellID int 1" << std::endl;
		filestr << "LOOKUP_TABLE default" << std::endl;
		
		filestr << std::scientific << std::setprecision(0);
		for( std::vector<UInt>::const_iterator it = idCells.begin();
		     it != idCells.end(); ++it )
		{
			filestr << *it << " ";
		}
		filestr << std::endl;
		
		filestr.close();
		
		return 1;
	}

 	bool Mesh2D::exportApproxFractureNetworkUniqueVtk(const std::string & fileName) const
	{
		std::fstream filestr;
		
		filestr.open (fileName.c_str(), std::ios_base::out);
	
		if (filestr.is_open())
		{
			std::cout << std::endl << " File: " << fileName << ", successfully opened";
		}
		else
		{
			std::cerr << std::endl << " *** Error: file not opened *** " 
					  << std::endl << std::endl;
			return  0;
		}
		
		std::cout << std::endl <<
			" Exporting Mesh2D::ApproxFractureNetwork in Vtk format... "
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
				filestr << M_nodes[*jt].x() << " " << M_nodes[*jt].y() << " 0" << std::endl;
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
	
	bool Mesh2D::exportApproxFractureNetworkVtk(const std::string & fileName) const
	{
		bool status=1;
		
		std::cout << std::endl <<
			" Exporting CMesh2D::ApproxFaultNetwork in Vtk format... "
			<< std::endl;
		
		for( UInt i = 0; i < M_approxFractureNetwork.size(); ++i)
			status = status && this->exportApproxFractureVtk( fileName+"_ApproxFracture"+
							  boost::lexical_cast<std::string>(i)+
							  ".vtk", i );
		
		return status;
	}
	
	bool Mesh2D::exportApproxFractureVtk( const std::string & fileName,
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
			std::cerr << std::endl << " *** Error: file not opened *** "
					  << std::endl << std::endl;
			return  0;
		}
		
		std::cout << std::endl <<
			" Exporting CMesh2D::ApproxFracture[" << f << "] in Vtk format... "
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
			filestr << M_nodes[*it].x() << " " << M_nodes[*it].y() << " 0" << std::endl;
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
 
 
	// --------------------   Overloading operator< (for Edge2D)  --------------------
 
	bool operator<(const Geometry::Mesh2D::Edge2D & e1, const Geometry::Mesh2D::Edge2D & e2)
	{
		if(e1.getP1()<e2.getP1())
			return 1;
		if(e1.getP1()>e2.getP1())
			return 0;
		if(e1.getP2()<e2.getP2())
			return 1;
		return 0;
	}
 
} // namespace Geometry
