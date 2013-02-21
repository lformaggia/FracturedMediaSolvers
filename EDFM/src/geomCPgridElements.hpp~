 /*!
 *	@file geomCPgridElements.hpp
 *	@brief Classes for Corner Point Grid elements: pillar and cell.
 *
 *	@author Luca Turconi <lturconi@gmail.com>
 *  @date 22-09-2012
 *
 */ 

#ifndef GEOMCPGRIDELEMENTS_HPP_
#define GEOMCPGRIDELEMENTS_HPP_

#include<vector>
#include<iostream>
#include<string>
#include<fstream>

#include "TypeDefinition.hpp"
#include "geomPoint3D.hpp"
#include "geomLine.hpp"
#include "geomSegment.hpp"
#include "geomBilinearSurface.hpp"


namespace Geometry
{	
class CPgrid;	// Forward declaration

/*!
		@class CPpillar
		
		@author Luca Turconi <lturconi@gmail.com>
		
    	This class implements the concept of pillar grid.
    	
    	A pillar grid is a straight line associated with a Corner Point grid
    	and identified by two indexes i and j. According with this point of view
    	a pillar is build inheriting the basic funtionality available for a simple line.
    	
    	This class allows to save a pillar independently of a Corner Point grid.
    	It is also possible to get the segment of the pillar directly associated
    	with the grid domain.
    	
    	The method showMe allows to print the class informations.
    	
    */
class CPpillar: public Line{
public:
	//! @name Constructor & Destructor
	//@{
		
	//! Empty constructor
	CPpillar();
	
	//! Constructor, getting coordinates the generating points and grid indexes
	/*!
	 * @param xA The first point x-coordinate
	 * @param yA The first point y-coordinate
	 * @param zA The first point z-coordinate
	 * @param xB The second point x-coordinate
	 * @param yB The second point y-coordinate
	 * @param zB The second point z-coordinate
	 * @param i Index i in the grid numeration
	 * @param j Index j in the grid numeration
	 */
	CPpillar(const Real & xA, const Real & yA, const Real & zA,
			 const Real & xB, const Real & yB, const Real & zB,
			 const UInt & i=0, const UInt & j=0);
	
	//! Constructor, getting the generating points and grid indexes
	/*!
	 * @param a The first point
	 * @param b The second point
	 * @param i Index i in the grid numeration
	 * @param j Index j in the grid numeration
	 */
	CPpillar(const Point3D & a, const Point3D & b, const UInt & i=0, const UInt & j=0);
	
	//! Destructor
	virtual ~CPpillar();
	
	//@}
	
	//! @name Get Methods
	//@{
		
	//! Get index i
	/*!
	 * @return The grid index i
	 */
	inline UInt getI() const
	{ return M_i; }
	
	//! Get index j
	/*!
	 * @return The grid index j
	 */
	inline UInt getJ() const
	{ return M_j; }
	
	//@}
	
	//! @name Methods
	//@{
		
	//! Get segment associated with grid domain
	/*!
	 * It computes the segment of the pillar associated with the given grid.
	 * @param grid The grid
	 * @return The segment of the pillar associated with the grid domain
	 */	
	Segment getLimitedPillar(const CPgrid & grid) const;
	
	//! Display general information about the content of the class
	/*!
	 * List of things displayed in the class
	 * @param out Specify the output format (std::cout by default)
	 */	
	void showMe(std::ostream  & out=std::cout) const;
	
	//@}
	
private:
	UInt M_i, M_j;
};


/*!
		@class CPcell
		
		@author Luca Turconi <lturconi@gmail.com>
		
    	This class implements the concept of Corner Point cell.
    	
    	A CPcell is an exahedral cell defined by its 8 vertices.
    	
    	This class allows to save a cell independently of a Corner Point grid.
    	It is also possible to get the vertices and the edges of the cell.
    	
    	The method showMe allows to print the class informations.
    	
    	For 3D visualization use the exportVtk method.
    	
    */
class CPcell{
public:
	//! @name Constructor & Destructor
	//@{
		
	//! Empty constructor
	CPcell();
	
	//! Constructor, getting vertices and indexes
	/*!
	 * @param v A vector containing the 8 vertices
	 * @param i The index i
	 * @param j The index j
	 * @param k The index k
	 * @param actnum The ACTNUM feature
	 */
	CPcell(const std::vector<Point3D> & v,
		   const UInt & i=0, const UInt & j=0, const UInt & k=0, const UInt & actnum=1);
	
	//! Copy constructor
	/*!
	 * @param c The Corner Point cell copied in the new object
	 */
	CPcell(const CPcell & c);
	
	//!Destructor
	~CPcell();
	
	//@}
	
	//! @name Get Methods
	//@{
		
	//! Get index i
	/*!
	 * @return The index i
	 */
	inline UInt i() const { return M_i; }
	
	//! Get index j
	/*!
	 * @return The index j
	 */
	inline UInt j() const { return M_j; }
	
	//! Get index k
	/*!
	 * @return The index k
	 */
	inline UInt k() const { return M_k; }
	
	//! Get ACTNUM feature
	/*!
	 * @return The ACTNUM feature
	 */
	inline UInt getActnum() const { return M_actnum; }
	
	//! Get vector of vertices
	/*!
	 * @return A references to the vector storing the vertices
	 */
	inline std::vector<Point3D> & getVerticesVector()
		{ return M_vertex; }
	
	//! Get vector of vertices (const version)
	/*!
	 * @return A constant references to the vector storing the vertices
	 */
	inline const std::vector<Point3D> & getVerticesVector() const
		{ return M_vertex; }
	
		
	//! Get single vertex
	/*!
	 * The vertices of the cell are numbered from 1 to 8.\n
	 * The convention adopted is the following:\n
	 * \n
	 * .......3----4......................... \n
	 * ....../|.../|............z............ \n
	 * .....1----2.|............|...y........ \n
	 * .....| 7--|-8............|../......... \n
	 * .....|/...|/.............|./.......... \n
	 * .....5----6..............o------x..... \n
	 *
	 * @return The selected vertex
	 */
    Point3D getVertex(const UInt & i) const;
		
	//! Get edge
	/*!
	 * The edges of the cell are numbered from 1 to 12. \n
	 * The convention adopted is the following: \n
	 * \n
	 * .........................................edge 1 = vertices 1->2 \n
	 * .........................o--3--o.........edge 2 = vertices 2->4 \n
	 * ......................../|..../|.........edge 3 = vertices 4->3 \n
	 * .......................4.8...2.7.........edge 4 = vertices 3->1 \n
	 * ....................../..|../..|.........edge 5 = vertices 1->5 \n
	 * .....z...............o--1--o...|.........edge 6 = vertices 2->6 \n
	 * .....|...y...........|...o-|11-o.........edge 7 = vertices 4->8 \n
	 * .....|../............5../..6../..........edge 8 = vertices 3->7 \n
	 * .....|./.............|.12..|.10..........edge 9 = vertices 5->6 \n
	 * .....o------x........|/....|/............edge 10 = vertices 6->8 \n
	 * .....................o--9--o.............edge 11 = vertices 8->7 \n
	 * .........................................edge 12 = vertices 7->5 \n
	 *
	 * @return The selected edge
	 */
	Segment getEdge(const UInt & i) const;
	
	//@}
	
	//! @name Methods
	//@{
		
	//! Export cell in vtk format
	/*!
	 * It generates the vtk file, for the 3D visualization of the Corner Point Cell.
	 * (Use paraview to open the vtk file)
	 * @param filename The name of the vtk file created by this method
	 * @return TRUE -> operation ended correctly
				FALSE -> an error occurred
	 */	
	bool exportVtk(const std::string & filename) const;
	
	
	//! Display general information about the content of the class
	/*!
	 * List of things displayed in the class
	 * @param out Specify the output format (std::cout by default)
	 */
	void showMe(std::ostream  & out=std::cout) const ;
	
	bool isIn(Point3D & ) const;

	bool intersectTheFace(Segment ,  int ,   Point3D & ) const ;

	std::vector<Point3D> segmentIntersectCell(Segment) const;

	double set_volume() const;

	//@}
	
private:
    std::vector<Point3D> M_vertex;
	UInt M_actnum;
	UInt M_i, M_j, M_k;
};

} // namespace Geometry

#endif /* GEOMCPGRIDELEMENTS_HPP_ */
