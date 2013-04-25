#ifndef HH_BBOX_HPP_HH
#define HH_BBOX_HPP_HH
#include <vector>
#include <iostream>
namespace ADT
{
/*! A simple structure that defines a bounding box.
 *  \author Luca Formaggia 2013
 */
template <unsigned int PhysicalDimensions>
class BBox
{
public:
    //!@{ Constructors
    BBox (double const* xmin, double const* xmax)
    {
        for (unsigned int i = 0; i < PhysicalDimensions; ++i)
        {
            M_coordinates[i] = xmin[i];
            M_coordinates[PhysicalDimensions + i] = xmax[i];
        }
    }
    BBox (const std::vector<double> xmin, const std::vector<double> xmax) :
        BBox (& (xmin[0]), & (xmax[0]) ) {}

    BBox()
    {
        for (unsigned int i = 0; i < 2 * PhysicalDimensions; ++i)
        {
            M_coordinates[i] = 0.0;
        }
    }

    BBox (BBox const& bb)
    {
        this->set (bb.M_coordinates, bb.M_coordinates + PhysicalDimensions);
    }
    //!@}
    //! EXtract lowest left corner.
    std::vector<double> min() const;
    //! Extract upper right corner
    std::vector<double> max() const;
    //! Set values;
    inline void set (double const* xmin, double const* xmax);
    //! Extract value circulating with level.
    double operator [] (unsigned int level) const
    {
        return
            M_coordinates[level % (2 * PhysicalDimensions)];
    }
    //! Extract value circulating with level (non const version).
    double & operator [] (unsigned int level)
    {
        return
            M_coordinates[level % (2 * PhysicalDimensions)];
    }
    /*! Transform the bounding box according to a given transformation
      \f$ S\mathbf{x}\leftarrow S(\mathbf{x}-\mathbf{o})
    */
    void transform (double const* origin, double const* scalingFactor);
    BBox& operator = (BBox const& rhs)
    {
        this->set (rhs.M_coordinates, rhs.M_coordinates + PhysicalDimensions);
        return *this;
    }
    //! Checks if two bboxes intersect each other
    bool intersect (BBox const& rhs) const;

private:
    //! xmin,ymix,zmin,xmax...
    double M_coordinates[2 * PhysicalDimensions];
};

template <unsigned int PhysicalDimensions>
std::vector<double> BBox<PhysicalDimensions>::min() const
{
    return
        std::vector<double> (M_coordinates,
                             M_coordinates + PhysicalDimensions);
}

template <unsigned int PhysicalDimensions>
std::vector<double> BBox<PhysicalDimensions>::max() const
{
    return
        std::vector<double> (M_coordinates + PhysicalDimensions,
                             M_coordinates + 2 * PhysicalDimensions);
}

template <unsigned int PhysicalDimensions>
void BBox<PhysicalDimensions>::set (double const* xmin, double const* xmax)
{
    for (unsigned int i = 0; i < PhysicalDimensions; ++i)
    {
        M_coordinates[i] = xmin[i];
        M_coordinates[PhysicalDimensions + i] = xmax[i];
    }
}

template <unsigned int PhysicalDimensions>
void BBox<PhysicalDimensions>::transform (double const* origin, double const* scalingFactor)
{
    for (unsigned int i = 0; i < PhysicalDimensions; ++i)
    {
        M_coordinates[i] = (M_coordinates[i]
                            - origin[i]) * scalingFactor[i];
        M_coordinates[PhysicalDimensions + i] = (M_coordinates[PhysicalDimensions + i]
                                                 - origin[i]) * scalingFactor[i];
    }
}

template<unsigned int PhysicalDimensions>
bool BBox<PhysicalDimensions>::intersect (BBox const& rhs) const
{
    for (unsigned int i = 0; i < PhysicalDimensions; ++i)
    {
        if (M_coordinates[i] > rhs.M_coordinates[i + PhysicalDimensions])
        {
            return false;
        }
        if (rhs.M_coordinates[i] > M_coordinates[i + PhysicalDimensions])
        {
            return false;
        }
    }
    return true;
}

template <unsigned int PhysicalDimensions>
std::ostream& operator << (std::ostream& out, BBox<PhysicalDimensions>
                           const& box)
{
    out << " ********** BOUNDING  BOX *************" << std::endl;
    out << " Lower left: ";
    for (unsigned int i = 0; i < PhysicalDimensions; ++i)
    {
        out << box[i] << " ";
    }
    out << std::endl;
    out << " Upper right: ";
    for (unsigned int i = 0; i < PhysicalDimensions; ++i)
    {
        out << box[i + PhysicalDimensions] << " ";
    }
    out << std::endl;
    return out;
}

} // end namespace
#endif
