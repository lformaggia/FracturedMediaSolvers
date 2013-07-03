/*!
 *	@file IntersectionSolver.cpp
 *	@brief Wrapper class for solve intersection problems.
 *
 *	@author Luca Turconi <lturconi@gmail.com>
 *  @date 22-09-2012
 *
 */


#include "IntersectionSolver.hpp"

// --------------------   Class IntersectionSolver   --------------------

// ==================================================
// Constructors & Destructor
// ==================================================
IntersectionSolver::IntersectionSolver(IntersectionAlgorithm alg, GridStrategy strategy) :
  M_alg(alg), M_strategy(strategy), M_toll(1e-8), M_maxIter(60), M_stdDivision(1){}
		
IntersectionSolver::~IntersectionSolver() {}

// ==================================================
// Methods
// ==================================================
void IntersectionSolver::showMe(std::ostream  & out)
{
	out << " INTERSECTION SOLVER PARAMETER: "<< std::endl;
	
	switch(M_alg)
	{
		case NEWTON :
			out << "  Algorithm : NEWTON" << std::endl;
			break;
		case APPROX :
			out << "  Algorithm : APPROX" << std::endl;
			break;
		case APPROXNEWTON :
			out << "  Algorithm : APPROXNEWTON" << std::endl;
			break;
		case APPROXREFINED :
			out << "  Algorithm : APPROXREFINED" << std::endl;
			break;
	}
	
	
switch(M_strategy)
	{
		case FOR3 :
			out << "  GridStrategy : FOR3" << std::endl;
			break;
		case EDGEMAP :
			out << "  GridStrategy : EDGEMAP" << std::endl;
			break;
		case FOR3OPT :
			out << "  GridStrategy : FOR3OPT" << std::endl;
			break;
	}
	
	out << "  Tolerance = " << M_toll << std::endl;
	out << "  Max Iteration = " << M_maxIter << std::endl;
	out << "  stdDivision = " << M_stdDivision << std::endl;
	
}

// ==================================================
// Operators
// ==================================================
void IntersectionSolver::operator()
	(Fault & f, const CPgrid & g, Intersect::GridIntersections & gridInter) const
{
  // For the optimised option we need the tree search structure
	switch(M_alg)
	{
		case NEWTON :
			switch(M_strategy)
			{
				case FOR3 :
					f.newtonIntersectionWithGrid_FOR3(g,gridInter,M_toll,M_maxIter);
					break;
			        case FOR3OPT :
				  /* The first time the method is called we need to set up
				     the range search tree.  
				     @warning If we operate on a new grid we must create a new IntersectionSolver object!
				  */
				    f.newtonIntersectionWithGrid_FOR3OPT(g,gridInter,M_toll,M_maxIter);
					break;
				case EDGEMAP :
					f.newtonIntersectionWithGrid_EDGEMAP(g,gridInter,M_toll,M_maxIter);
					break;
			}
			break;
		case APPROX :
			switch(M_strategy)
			{
				case FOR3 :
					f.approxIntersectionWithGrid_FOR3(g,gridInter,M_stdDivision);
					break;
				case EDGEMAP :
					f.approxIntersectionWithGrid_EDGEMAP(g,gridInter,M_stdDivision);
					break;
			}
			break;
		case APPROXNEWTON :
			switch(M_strategy)
			{
				case FOR3 :
					f.approxNewtonIntersectionWithGrid_FOR3(g,gridInter,M_toll,M_maxIter);
					break;
				case EDGEMAP :
					f.approxNewtonIntersectionWithGrid_EDGEMAP(g,gridInter,M_toll,M_maxIter);
					break;
			}
			break;
		case APPROXREFINED :
			switch(M_strategy)
			{
				case FOR3 :
					f.approxRefinedIntersectionWithGrid_FOR3(g,gridInter,M_toll);
					break;
				case EDGEMAP :
					f.approxRefinedIntersectionWithGrid_EDGEMAP(g,gridInter,M_toll);
					break;
			}
			break;
	}
}
		
void IntersectionSolver::operator()
	(Fault & f, const CPcell & c, Intersect::CellIntersections & cellInter) const
{
	switch(M_alg)
	{
		case NEWTON :
			f.newtonIntersectionWithCell(c,cellInter,M_toll,M_maxIter);
			break;
		case APPROX :
			f.approxIntersectionWithCell(c,cellInter,M_stdDivision);
			break;
		case APPROXNEWTON :
			f.approxNewtonIntersectionWithCell(c,cellInter,M_toll,M_maxIter);
			break;
		case APPROXREFINED :
			f.approxRefinedIntersectionWithCell(c,cellInter,M_toll);
			break;
	}
}

Point3D IntersectionSolver::operator()
	(Fault & f, const Segment & s) const
{
	Point3D p;
	switch(M_alg)
	{
		case NEWTON :
			p = f.newtonIntersectionWith(s,M_toll,M_maxIter);
			break;
		case APPROX :
			p = f.approxIntersectionWith(s,M_stdDivision);
			break;
		case APPROXNEWTON :
			p = f.approxNewtonIntersectionWith(s,M_toll,M_maxIter);
			break;
		case APPROXREFINED :
			p = f.approxRefinedIntersectionWith(s,M_toll);
			break;
	}
	return p;
}
