 /*!
 *	@file IntersectionSolver.hpp
 *	@brief Wrapper class for solve intersection problems.
 *
 *	@author Luca Turconi <lturconi@gmail.com>
 *  @date 22-09-2012
 *
 */


#ifndef INTERSECTIONSOLVER_HPP_
#define INTERSECTIONSOLVER_HPP_

#include "TypeDefinition.hpp"
#include "geomPoint3D.hpp"
#include "geomSegment.hpp"
#include "geomCPgridElements.hpp"
#include "geomCPgrid.hpp"
#include "geomFault.hpp"
#include "interCellIntersections.hpp"
#include "interGridIntersections.hpp"

using namespace Geometry;

	/*!
		@enum IntersectionAlgorithm
    	Enum type for algorithm selection:
    	 - NEWTON : exact intersection
    	 - APPROX : 2 triangles surface approximation
    	 - APPROXNEWTON : 2 triangles approximated test + newton solver
    	 - APPROXREFINED : intersection with surface triangulation
    */
enum IntersectionAlgorithm {NEWTON, APPROX, APPROXNEWTON, APPROXREFINED};

	/*!
		@enum GridStrategy
    	Enum type to chose the grid scanning strategy:
    	 - FOR3 : 3 for loops, cell centered (each intersection is computed several times)
    	 - EDGEMAP : uses auxiliary containers to avoid repetitions
    */
enum GridStrategy {FOR3, FOR3OPT, EDGEMAP};

	/*!
		@class IntersectionSolver
		
		@author Luca Turconi <lturconi@gmail.com>
		
		This class is a wrapper for different algorithm.
		It allows to call fault methods simply and intuitively.
    	
    	The method showMe prints the class informations.
    	
    */
class IntersectionSolver {
public:
	//! @name Constructors & Destructor
	//@{
		
	//! Constructor, getting algorithm type
	IntersectionSolver(IntersectionAlgorithm alg, GridStrategy strategy);
	
	//! Destructor
	~IntersectionSolver();
	
	//@}
	
	//! @name Get Methods
	//@{
	
	//! Get tolerance
	/*!
	 * @return The tolerance
	 */
	inline Real getTolerance() const
		{ return M_toll; }
	
	//! Get max iterations number
	/*!
	 * @return The max iterations number
	 */
	inline UInt getMaxIteration() const
		{ return M_maxIter; }
	
	//! Get stdDivision
	/*!
	 * @return The stdDivision
	 */
	inline bool getStdDivision() const
		{ return M_stdDivision; }
	
	//@}
	
	//! @name Set Methods
	//@{
	
	//! Set tolerance
	/*!
	 * @param toll The new tolerance value
	 */
	inline void setTolerance(const Real & toll)
		{ M_toll = toll; }
	
	//! Set max iterations number
	/*!
	 * @param maxIter The new M_maxIter value
	 */
	inline void setMaxIteration(const UInt & maxIter)
		{ M_maxIter = maxIter; }
	
	//! Set stdDivision
	/*!
	 * @param stdDivision The new stdDivision value
	 */
	inline void setStdDivision(const bool & stdDivision)
		{ M_stdDivision = stdDivision; }

	//@}
	
	//! @name Methods
	//@{
		
	//! Display general information about the content of the class
	/*!
	 * List of things displayed in the class
	 * @param out Specify the output format (std::cout by default)
	 */
	void showMe(std::ostream  & out=std::cout);
	
	//@}
	
	//! @name Operators
	//@{
	
	//! operator() for solve Fault-Grid intersections
	/*!
	 * Compute the Fault-Grid intersections using the selected algorithm
	 * @param f The fault
	 * @param g The grid
	 * @param gridInter Where the intersections are saved
	 */
	void operator()(Fault & f, const CPgrid & g, Intersect::GridIntersections & gridInter) const;
	
	//! operator() for solve Fault-Cell intersections
	/*!
	 * Compute the Fault-Cell intersections using the selected algorithm
	 * @param f The fault
	 * @param c The cell
	 * @param cellInter Where the intersections are saved
	 */
	void operator()(Fault & f, const CPcell & c, Intersect::CellIntersections & cellInter) const;
	
	//! operator() for solve Fault-Segment intersection
	/*!
	 * Compute the Fault-segment intersection using the selected algorithm
	 * @param f The fault
	 * @param s The segment
	 * @return The intersection point
	 */
	Point3D operator()(Fault & f, const Segment & s) const;
	
	//@}
	
private:
	IntersectionAlgorithm M_alg;
	GridStrategy M_strategy;
	Real M_toll;
	UInt M_maxIter;
	bool M_stdDivision;
};

 
#endif /* INTERSECTIONSOLVER_HPP_ */
