#include "domain.hpp"
#include "geomCPgrid.hpp"
namespace ADT
{
template<>
Domain<Geometry::CPgrid>::Domain (Geometry::CPgrid const& cPGrid)
{
    int ndimp = 3;
    origin.resize (6);
    scalingfactors.resize (6);
    // T is equal to Box<NDIMP>.
    const std::vector<Real>& cPCoord = cPGrid.coord();
    unsigned int numPoints = cPCoord.size() / 3;
    std::vector<Real> coord[3];
    for (unsigned int i = 0; i < 3; ++i)
    {
        coord[i].reserve (numPoints);
    }
    for (unsigned int i = 0; i < numPoints; ++i)
    {
        unsigned int pos = 3 * i;
        coord[0].push_back (cPCoord[pos + 0]);
        coord[1].push_back (cPCoord[pos + 1]);
        coord[2].push_back (cPCoord[pos + 2]);
    }
    for (int i = 0; i < ndimp; ++i)
    {
        origin[i] = * (std::min_element (coord[i].begin(), coord[i].end() ) );

        /* This statement is necessary when representing a
         * rectangle with the corner having minimum coordinates
         * and the opposite one.
         */
        scalingfactors[i] = * (std::max_element (coord[i].begin(), coord[i].end() ) );

        // Add the tolerance.
        double delta = scalingfactors[i] - origin[i];
        origin[i] -= delta * this->gettolerance();
        scalingfactors[i] += delta * this->gettolerance();

        delta = scalingfactors[i] - origin[i];
        scalingfactors[i] = 1. / std::max (delta, this->getmindiff() );

        /* Repeat the limits because tree dimension is in fact 2 *
         * physical space dimension because because a box is
         * defined by two points. Only for convenience
         */
        origin[i + ndimp] = origin[i];
        scalingfactors[i + ndimp] = scalingfactors[i];
    }
}

}// end namespace ADT
