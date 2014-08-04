#ifndef FVCODE3D_BOUNDINGBOX_HPP_
#define FVCODE3D_BOUNDINGBOX_HPP_

#include <FVCode3D/geometry/Point3D.hpp>

namespace FVCode3D
{

//! Bounding Box that contains the data
/*!
 * @Class BoundingBox
 * (It is also known as BB)
 * This class contains the parameters that define the bounding box.
 * All the data is contained in the extremes of the bounding box.
 */
class BoundingBox
{
public:

    //! Empty Constructor
    BoundingBox();

    //! Default copy-constructor
    BoundingBox(const BoundingBox &) = default;

    //! Constructor
    /*!
     * Set the extremes of the bounding box by 6 values
     */
    BoundingBox(const Real xMin, const Real xMax, const Real yMin, const Real yMax, const Real zMin, const Real zMax);

    //! Constructor
    /*!
     * Set the extremes of the bounding box by two points
     */
    BoundingBox(const Point3D & pMin, const Point3D & pMax);

    //! Constructor
    /*!
     * Set the extremes of the bounding box by a vector of nodes
     */
    BoundingBox(const std::vector<Point3D> & nodes);

    //! @name Get extremes of BB
    //@{
    Real xMin() const { return M_xMin; }
    Real xMax() const { return M_xMax; }
    Real yMin() const { return M_yMin; }
    Real yMax() const { return M_yMax; }
    Real zMin() const { return M_zMin; }
    Real zMax() const { return M_zMax; }
    //@}

    //! Get diagonal of the BB
    /*!
     * @return Diagonal of the BB
     */
    Real diagonal() const { return M_L; }

    //! @name Get scaling parameters
    //@{
    Real mx() const { return M_mx; }
    Real my() const { return M_my; }
    Real mz() const { return M_mz; }
    Real qx() const { return M_qx; }
    Real qy() const { return M_qy; }
    Real qz() const { return M_qz; }
    //@}

    //! @name Set extremes of BB
    //@{
    void setXMin(const Real x) { M_xMin = x; }
    void setXMax(const Real x) { M_xMax = x; }
    void setYMin(const Real y) { M_yMin = y; }
    void setYMax(const Real y) { M_yMax = y; }
    void setZMin(const Real z) { M_zMin = z; }
    void setZMax(const Real z) { M_zMax = z; }
    void setExtremes( const Real xmin, const Real xmax,
                      const Real ymin, const Real ymax,
                      const Real zmin, const Real zmax);
    //@}

    //! Set the diagonal of the BB
    /*!
     * Compute the diagonal of the BB
     */
    void computeDiagonal();

    //! Set the scaling parameters
    void computeScalingParameters();

    //! Scale node such that it is in the unit BB
    /*!
     * @param node point to be scaled
     * @pre set the BB
     */
    void scaleNodesToUnit(Point3D & node) const;

    //! Scale nodes such that all of them are in the unit BB
    /*!
     * @param nodes vector of points to be scaled
     * @pre set the BB
     */
    void scaleNodesToUnit(std::vector<Point3D> & nodes) const;

    //! Scale nodes such that it is in the original BB
    /*!
     * @param node point to be scaled
     * @pre set the BB
     */
    void scaleNodesToPhysical(Point3D & node) const;

    //! Scale nodes such that all of them are in the original BB
    /*!
     * @param nodes vector of points to be scaled
     * @pre set the BB
     */
    void scaleNodesToPhysical(std::vector<Point3D> & nodes) const;

    //! Destructor
    ~BoundingBox() = default;

private:

    //! @name Bounding box dimensions
    //@{
    Real M_xMin;
    Real M_xMax;
    Real M_yMin;
    Real M_yMax;
    Real M_zMin;
    Real M_zMax;
    //@}

    //! Diagonal of the bounding box
    Real M_L;

    //! @name Scaling parameters
    //@{
    Real M_mx;
    Real M_my;
    Real M_mz;
    Real M_qx;
    Real M_qy;
    Real M_qz;
    //@}
};

} // namespace FVCode3D

#endif // FVCODE3D_BOUNDINGBOX_HPP_
