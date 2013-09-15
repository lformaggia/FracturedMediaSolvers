/*!
 *	@file DarcyQuadratureRules.hpp
 *	@brief Classes for quadrature rules.
 *
 *	@author Francesco Della Porta 
 *
 */ 
#ifndef __DARCYQUADRATURERULE_HPP__
#define __DARCYQUADRATURERULE_HPP__

#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>
#include "DarcyTypeDefinitions.hpp"
	

namespace Darcy{

/*!
	@class QuadratureRule
	This class is a base class. The aim of this class is to be a base from which derive classes that are quadrature rules
*/
template <class T> 
class QuadratureRule
{
public:
typedef typename std::unique_ptr<QuadratureRule<T> > QuadratureRuleHandler;
//! Method clone for the class
/*!
	It is possible to pass the QuadratureRule to the class Quadrature without having problems
*/
virtual std::unique_ptr<QuadratureRule<T>> clone() const =0;
//! Method apply
/*!
	@param pointVector A vector with the vertexes of a cell
	@param _volume The volume of the cell
	@param Integrand The function we want to integrate
	@return the approximation of the integral of Integrand on the passed cell
*/
virtual D_Real apply(const std::vector<Geometry::Point2D>& pointVector, const D_Real _volume, const std::function<D_Real(Geometry::Point2D)>& Integrand) const = 0;
};


/*!
	@class Triangle2D
	This class is a quadrature rule for triangular meshes in 2D. It has order 2
*/
class Triangle2D  : public QuadratureRule<Geometry::Dimension<2> >
{
public:
	//! constructor
	Triangle2D() = default;
	//! clone method
    QuadratureRuleHandler clone() const;
	//! destructor
    ~Triangle2D(){};
	//! apply method
	D_Real apply (const std::vector<Geometry::Point2D>& pointVector, const D_Real _volume, const std::function<D_Real(Geometry::Point2D)>& Integrand) const;
};


/*!
	@class CentroidQuadrature
	This class is a quadrature rule for arbitrary polygons in 2D or 3D
*/
template<class T>
class CentroidQuadrature : public QuadratureRule<T>
{
public:
	typedef typename T::Generic_Point Generic_Point;
	//! constructor
	CentroidQuadrature() = default;
	//! clone method
    std::unique_ptr<QuadratureRule<T> > clone() const;
	//! destructor
    ~CentroidQuadrature(){};
	//! apply method
	D_Real apply (const std::vector<Generic_Point>& pointVector, const D_Real _volume, const std::function<D_Real(Generic_Point)>& Integrand) const;
};



// --------------------   Class Triangle2D   --------------------
 
// ==================================================
// Methods
// ==================================================

QuadratureRule<Geometry::Dimension<2> >::QuadratureRuleHandler Triangle2D::clone()const { return 
   QuadratureRuleHandler(new Triangle2D(*this));}


D_Real Triangle2D::apply (const std::vector<Geometry::Point2D>& pointVector, const D_Real _volume, const std::function<D_Real(Geometry::Point2D)>& Integrand) const
{
	D_Real q_integral = 0;
	Geometry::Point2D baricenter = CGAL::ORIGIN;

	UInt bit = 2;
	for (UInt it = 0; it < pointVector.size(); ++it)
	{	
		baricenter = baricenter + (pointVector[bit]-CGAL::ORIGIN)/2;
		baricenter = baricenter + (pointVector[it]-CGAL::ORIGIN)/2.;
		q_integral += Integrand(baricenter)/3.;
		baricenter = CGAL::ORIGIN;
		bit = it;
	}
	return q_integral*_volume;
}


// --------------------   Class CentroidQuadrature   --------------------
 
// ==================================================
// Methods
// ==================================================

template<class T>
std::unique_ptr<QuadratureRule<T> > CentroidQuadrature<T>::clone()const { return 
   std::unique_ptr<QuadratureRule<T> >(new CentroidQuadrature(*this));}


template<class T>
D_Real CentroidQuadrature<T>::apply(const std::vector<Generic_Point>& pointVector, const D_Real _volume, const std::function<D_Real(Generic_Point)>& Integrand) const
{
	D_Real q_integral = 0;
	Generic_Point sum = CGAL::ORIGIN;
	UInt N = pointVector.size();

	for ( UInt j = 0; j < N; ++j)
		sum = sum + (pointVector[j]-CGAL::ORIGIN)/N;

	q_integral += Integrand(sum);

	return q_integral*_volume;
}

}
#endif
