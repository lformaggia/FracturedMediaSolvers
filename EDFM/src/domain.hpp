/**
 *  \file domain.hpp
 *  A modified version of the code produced by Pigoli Davide and Prada Daniele
 *  \author Luca Formaggia 2013
 */

#ifndef DOMAIN_HPP_
#define DOMAIN_HPP_

#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include "geomCPgrid.hpp"

namespace ADT
{
/** \class Domain
 *  \brief It defines geometric limits of a set of points.
 */
template<class Shape>
class Domain
{
protected:
    /// Origin of the object's bounding box = min(coord(*,1:number of points)).
    std::vector<double> origin;
    /// Scaling factors = 1./(max(coord(*,1:number of points)) - min(coord(*,1:number of points))).
    std::vector<double> scalingfactors;
    /// Tolerance being applied to the object's bounding box.
    static double tolerance;
    /// Minimum difference between coordinates allowed.
    static double mindiff;
public:
    /** Default constructor.
     *
     *  It's fundamental in creating an ADTree object from a MeshFile::ff2dmesh or a MeshFile::ff3dmesh object.
     */
    Domain() {};
    /** Another constructor.
     *
     *  \param[in] coord Point coordinates organized in a matrix (coord[i, j] is the i-th coordinate of the j-th point). \n
     *  This matrix is created through a vector of vectors in order to use standard algorithms.
     *
     *  It finds geometric limits of a set of points first. Then adds a tolerance.
     *  Repeats the limits if the tree dimension is 2 * physical space dimension.
     *  This is an useful trick. For example, when you have to scale dimensions.
     */
    //    Domain (std::vector<std::vector<double> > const& coord);
    /*! Constructor that takes a corner point grid.
      Generic template version. Is activated if Shape=CPgrid
     */
    Domain (Shape const& cpgrid);
    /// Sets the tolerance being applied to the object's bounding box.
    inline static void settolerance (double const& tol)
    {
        tolerance = tol;
    }
    /// Gets the tolerance being applied to the object's bounding box.
    inline static double gettolerance()
    {
        return tolerance;
    }
    /// Sets the minimum difference between coordinates allowed.
    inline static void setmindiff (double const& md)
    {
        mindiff = md;
    }
    /// Gets the minimum difference between coordinates allowed.
    inline static double getmindiff()
    {
        return mindiff;
    }
    /// Gets the i-th coordinate of the origin of the domain's bounding box.
    std::vector<double> const& orig () const
    {
        return origin;
    }
    /// Gets the i-th scaling factor of the domain's bounding box.
    std::vector<double> const& scal () const
    {
        return scalingfactors;
    }
    /** Output operator.
     *
     *  It outputs the bounding box.
     */
    template<class S>
    friend std::ostream& operator<< (std::ostream&, Domain<S> const&);
};

template<class T>
double Domain<T>::tolerance = 1.e-3;

template<class T>
double Domain<T>::mindiff = std::numeric_limits<double>::min();

// template<class T>
// Domain<T>::Domain (std::vector<std::vector<double> > const& coord)
// {
//     int ndimp = T::dp();
//     origin.resize (T::dt() );
//     scalingfactors.resize (T::dt() );

//     /* Find geometric limits.
//      *
//      * If loops are put outside for loops in order to improve performance.
//      */
//     if (ndimp == int (coord.size() ) )
//     {
//         // T is equal to Point<NDIMP>, Triangle<NDIMP> or Tetrahedron.
//         if (T::dp() == T::dt() )
//         {
//             // T is equal to Point<NDIMP>
//             for (int i = 0; i < ndimp; ++i)
//             {
//                 origin[i] = * (std::min_element (coord[i].begin(), coord[i].end() ) );
//                 scalingfactors[i] = * (std::max_element (coord[i].begin(), coord[i].end() ) );

//                 // Add the tolerance.
//                 double delta = scalingfactors[i] - origin[i];
//                 origin[i] -= delta * gettolerance();
//                 scalingfactors[i] += delta * gettolerance();

//                 delta = scalingfactors[i] - origin[i];
//                 scalingfactors[i] = 1. / std::max (delta, getmindiff() );
//             }
//         }
//         else
//         {
//             // T is equal to Triangle<NDIMP> or Tetrahedron.
//             for (int i = 0; i < ndimp; ++i)
//             {
//                 origin[i] = * (std::min_element (coord[i].begin(), coord[i].end() ) );
//                 scalingfactors[i] = * (std::max_element (coord[i].begin(), coord[i].end() ) );

//                 // Add the tolerance.
//                 double delta = scalingfactors[i] - origin[i];
//                 origin[i] -= delta * gettolerance();
//                 scalingfactors[i] += delta * gettolerance();

//                 delta = scalingfactors[i] - origin[i];
//                 scalingfactors[i] = 1. / std::max (delta, getmindiff() );

//                 /* Repeat the limits because tree dimension is in fact 2 * physical space dimension
//                  * because the tree contains triangle or tetrahedron bounding boxes.
//                  */
//                 origin[i + ndimp] = origin[i];
//                 scalingfactors[i + ndimp] = scalingfactors[i];
//             }
//         }
//     }
//     else
//     {
//         // T is equal to Box<NDIMP>.
//         for (int i = 0; i < ndimp; ++i)
//         {
//             origin[i] = * (std::min_element (coord[i].begin(), coord[i].end() ) );

//             /* This statement is necessary when representing a rectangle with the corner having
//              * minimum coordinates and the opposite one.
//              */
//             scalingfactors[i] = * (std::max_element (coord[i + ndimp].begin(), coord[i + ndimp].end() ) );

//             // Add the tolerance.
//             double delta = scalingfactors[i] - origin[i];
//             origin[i] -= delta * gettolerance();
//             scalingfactors[i] += delta * gettolerance();

//             delta = scalingfactors[i] - origin[i];
//             scalingfactors[i] = 1. / std::max (delta, getmindiff() );

//             /* Repeat the limits because tree dimension is in fact 2 * physical space dimension
//              * because because a box is defined by two points.
//              */
//             origin[i + ndimp] = origin[i];
//             scalingfactors[i + ndimp] = scalingfactors[i];
//         }
//     }
// }




template<class T>
std::ostream& operator<< (std::ostream& ostr, Domain<T> const& d)
{
    ostr << std::endl << std::endl;
    ostr << "Domain" << std::endl;
    ostr << "------" << std::endl << std::endl;

    int dimp (3);
    for (int i = 0; i < dimp; ++i)
    {
        ostr << "x" << i + 1 << "_min = " << d.origin[i] << std::endl;
    }
    ostr << std::endl;

    for (int i = 0; i < dimp; ++i)
    {
        ostr << "x" << i + 1 << "_max = " << d.origin[i] + 1. / d.scalingfactors[i] << std::endl;
    }
    ostr << std::endl;

    return ostr;
}
}// end namespace ADT


#endif /* DOMAIN_HPP_ */
