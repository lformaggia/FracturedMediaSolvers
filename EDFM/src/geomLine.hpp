/*!
 *  @file geomLine.hpp
 *  @brief Class for Line in 3D space.
 *
 *  @author Luca Turconi <lturconi@gmail.com>
 *  @date 22-09-2012
 *
 */

#ifndef GEOMLINE_HPP_
#define GEOMLINE_HPP_

#include<iostream>

#include "TypeDefinition.hpp"
#include "geomPoint3D.hpp"
#include "geomSegment.hpp"


namespace Geometry
{

  /*!
    @class Line

    @author Luca Turconi <lturconi@gmail.com>

      This class implement the concept of Straight Line.
      It stores the two points generating the line.

      The class provides methods to get the point of the line
      with a fixed coordinate.

      The method showMe allows to print the class informations.
      It is also possible to print the attributes of the line by using << operator.

    */
  class Line
  {
  public:
    //! @name Constructor & Destructor
    //@{

    //! Empty constructor
    Line();

    //! Constructor, getting coordinates the generating points
    /*!
     * @param xA The first point x-coordinate
     * @param yA The first point y-coordinate
     * @param zA The first point z-coordinate
     * @param xB The second point x-coordinate
     * @param yB The second point y-coordinate
     * @param zB The second point z-coordinate
     */
    Line (const Real& xA, const Real& yA, const Real& zA,
          const Real& xB, const Real& yB, const Real& zB);

    //! Constructor, getting the generating points
    /*!
     * @param a The first point
     * @param b The second point
     */
    Line (const Point3D& a, const Point3D& b);

    //! Constructor, getting the generating segment
    /*!
     * @param s The segment copied in the new object
     */
    Line (const Segment& s);

    //! Copy constructor
    /*!
     * @param l The line copied in the new object
     */
    Line (const Line& l);

    //! Destructor
    virtual ~Line();

    //@}

    //! @name Get Methods
    //@{

    //! Get point A
    /*!
     * @return The generating point A
     */
    inline Point3D A() const
    {
      return M_pA;
    }

    //! Get point B
    /*!
     * @return The generating point B
     */
    inline Point3D B() const
    {
      return M_pB;
    }

    //@}

    //! @name Set Methods
    //@{

    //! Set point A
    /*!
     * @param p The new value for the generating point A
     */
    inline void setA (const Point3D& p)
    {
      M_pA = p;
    }

    //! Set point B
    /*!
     * @param p The new value for the generating point B
     */
    inline void setB (const Point3D& p)
    {
      M_pB = p;
    }

    //@}

    //! @name Methods
    //@{

    //! Get line point at x
    /*!
     * Given a x-coordinate value, this method return the corresponding point on the line.
     * @param x The x-coordinate value
     * @return The point of the line corresponding to the x-coordinate value
     */
    Point3D getPointAtX (const Real& x) const;

    //! Get line point at y
    /*!
     * Given a y-coordinate value, this method return the corresponding point on the line.
     * @param y The y-coordinate value
     * @return The point of the line corresponding to the y-coordinate value
     */
    Point3D getPointAtY (const Real& y) const;

    //! Get line point at z
    /*!
     * Given a z-coordinate value, this method return the corresponding point on the line.
     * @param z The z-coordinate value
     * @return The point of the line corresponding to the z-coordinate value
     */
    Point3D getPointAtZ (const Real& z) const;

    //! Display general information about the content of the class
    /*!
     * List of things displayed in the class
     * @param out Specify the output format (std::cout by default)
     */
    virtual void showMe (std::ostream& out = std::cout) const;

    //@}

  protected:
    Point3D M_pA;
    Point3D M_pB;
  };

  //! @name External Operators
  //@{

  //! The insertion operator
  /*!
   * @param ostr The stream object on which the action is performed.
   *      This is the first parameter of the global functions,
   *      and represents the object to the left of the operator,
   *      i.e. the object on which the extraction operation is performed.
   * @param l The segment inserted on the stream.
   * @return The stream object on which the action is performed (ostr)
   */
  std::ostream& operator<< (std::ostream& ostr, const Line& l);

  //@}


} // namespace Geometry

#endif /* GEOMLINE_HPP_ */