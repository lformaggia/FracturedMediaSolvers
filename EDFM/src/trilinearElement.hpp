#ifndef HH_TRILINEARELEMENT_HH
#define HH_TRILINEARELEMENT_HH
#include "trilinearShape.hpp"
#include "geomCPgridElements.hpp"
#include "gmm/gmm.h"
#include <iosfwd>
namespace Geometry{
  //! Options for the inverse map problem.
  struct InvMapOption
  {
    //! Tolerance on the iterates.
    /*!
      The iteration stops when \f$||x^{k+1}-x^{k}||<tol\f$.
     */
    double iteratesTolerance;
    /*! Tolerance on the residual.
      The iteration stops when \f$||F(x^{k})||<tol\f$.
     */
    double residualTolerance;
    //! Max number of iterations.
    unsigned int maxIter;
    
    /*! Bounding box tolerance.
      If the seeked point is outside a bounding box the newton method
      is simply stopped
     */
    double bboxTolerance;
  };
    

  //! The result of the algorithm for the inverse mapping.
  struct InvMapResult{
    //! The point in the reference space.
    Point3D finalPoint;
    //! True if doint is inside the reference hexa
    bool inside;
    //! indicates the direction of the target point w.r.t. element.
    /*!
      direction[i]>0 Point is on the right of element along direction i
      direction[i]=0 Point is inside the element along direction i
      direction[i]<0 Point is on the left of the element along direction i
    */
    int direction[3];
    bool converged;
    bool numIter;
  };
  //! To help debugging
  std::ostream & operator <<(std::ostream & out, InvMapResult const & res);

  class InverseMapping;

    /*! A wrapper class for a CPCell.
      It is required since the numbering of the nodes in a CPcell
      is different from that assumed in a TrilinearElement (the latter was required
      to simplify the use of expression templates.

    */
    class TrilinearElement
    {
    public:
      //! To do the inverse mapping I rely on a proxy.
      friend class InverseMapping;
      TrilinearElement(const CPcell & cell);
      //! Points given with the convention of a trilinear element (not CPcell!).
      /*!
	See the documentation of trilinearShape.
       */
      template<class LinearContainer>
      TrilinearElement(const LinearContainer points);
      //! The map from reference to physical space.
      /*! 
       * We use xx,yy and zz to indicate coordinates in the reference space.
       */
      Point3D map(double const & xx, double const & yy, double const & zz)const;
      //! The Jacobian in a point in the reference space.
      /*!
	\f$ J_{ij}=\frac{\partial x_i}{\partial xx_j} \f$
       */
      gmm::dense_matrix<double> const & J(double const & xx, double const & yy, double const & zz) const;
    private:
      Point3D M_points[8];
      //! I store the Jacobian to save memory copy
      mutable gmm::dense_matrix<double> my_J;
    };
  
    template<class LinearContainer>
    TrilinearElement::TrilinearElement(const LinearContainer points):my_J(3,3)
    {
      for (unsigned int i=0;i<8;++i)M_points[i]=points[i];
    }
  
  //! Class for inverse mapping.
  /*!  The methods for inverse mapping have not been stored in
    TrilinearElement in order to avoid loading the former class of too
    meny methods.  In fact inverse mapping requires the solution of a
    non-linear system, which has been enucleated in this class.  
    It is
    however a wrapper arount a TrilinearElement, it means that it
    requires a TrilinearElement object to exist!
    
    How to use it:
    \code
    TrilinearElement foo;
    ...
    InverseMapping proxy(foo);
    InvMapResult result=proxy(0.5,6.7,8.9);
    \endcode
    
  */
  class InverseMapping{
  public:
    InverseMapping(TrilinearElement const & element);
    static InvMapOption getInvMapOption()
    {
      return M_options;
    }
    static void  setInvMapOption(InvMapOption const & option)
    {
      M_options=option;
      }
    InvMapResult operator () (double const & x, double const & y, double const & z)const;
  private:
    TrilinearElement const & M_element;
    double M_BBmin[3];
    double M_BBmax[3];
    Point3D M_target;
    void F(const Point3D & x, const Point3D & target, Point3D & Fx)const;
  public:
    static InvMapOption M_options;
  };
  
  }//end namespace Geometry
#endif
