/*!
 * @file SolverHandler.cpp
 * @brief Factory that handles the creation of the Solver classes (definitions).
 */
#include <FVCode3D/solver/SolverHandler.hpp>

namespace FVCode3D
{

//! Initializes the static pointer to singleton
std::unique_ptr<SolverHandler> SolverHandler::theOnlySolverHandler = 0;

SolverHandler & SolverHandler::Instance()
{
    if (theOnlySolverHandler == 0)
    {
        theOnlySolverHandler.reset(new SolverHandler());
        theOnlySolverHandler->registration();
    }

    return *theOnlySolverHandler;
} // SolverHandler::Instance

void SolverHandler::registration()
{
    SolverProxy<EigenCholesky>  SolverChol("EigenCholesky");
    SolverProxy<EigenLU>        SolverLU("EigenLU");
#ifdef FVCODE3D_HAS_UMFPACK
    SolverProxy<EigenUmfPack>   SolverUmfPack("EigenUmfPack");
#endif // FVCODE3D_HAS_UMFPACK
    SolverProxy<EigenCG>        SolverCG("EigenCG");
    SolverProxy<EigenBiCGSTAB>  SolverBiCGSTAB("EigenBiCGSTAB");
#ifdef FVCODE3D_HAS_SAMG
    SolverProxy<Samg>           SolverSamg("Samg");
#endif // FVCODE3D_HAS_SAMG
} // SolverHandler::registration

void SolverHandler::addProduct(const std::string productName, const SolverBuilder & builder)
{
    M_productList.insert(std::make_pair(productName, builder));
} // SolverHandler::addProduct

SolverPtr_Type SolverHandler::getProduct(const std::string productName) const
{
    std::map<std::string,SolverBuilder>::const_iterator it = M_productList.find(productName);

    return (it == M_productList.end()) ? SolverPtr_Type() : it->second();
} // SolverHandler::getProduct

} // namespace FVCode3D
