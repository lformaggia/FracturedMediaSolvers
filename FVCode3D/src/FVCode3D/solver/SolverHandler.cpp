/*!
 * @file SolverHandler.cpp
 * @brief Factory that handles the creation of the Solver classes (definitions).
 */
#include <FVCode3D/solver/SolverHandler.hpp>

namespace FVCode3D
{

SolverHandler & SolverHandler::Instance()
{
    static SolverHandler theOnlySolverHandler;
    theOnlySolverHandler.registration();
    return theOnlySolverHandler;
} // SolverHandler::Instance

void SolverHandler::registration()
{
    if(M_productList.size() > 0)
        return;
    SolverProxy<EigenCholesky>  SolverChol("EigenCholesky");
    SolverProxy<EigenLU>        SolverLU("EigenLU");
    SolverProxy<EigenUmfPack>   SolverUmfPack("EigenUmfPack");
    SolverProxy<EigenCG>        SolverCG("EigenCG");
    SolverProxy<EigenBiCGSTAB>  SolverBiCGSTAB("EigenBiCGSTAB");
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
