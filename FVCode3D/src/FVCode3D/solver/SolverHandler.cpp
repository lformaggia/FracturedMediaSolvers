#include <FVCode3D/solver/SolverHandler.hpp>

namespace FVCode3D
{

//! Initializes the static pointer to singleton
SolverHandler * SolverHandler::theOnlySolverHandler = 0;

SolverHandler & SolverHandler::Instance()
{
    if (theOnlySolverHandler == 0)
    {
        theOnlySolverHandler = new SolverHandler();
        theOnlySolverHandler->registration();
    }

    return *theOnlySolverHandler;
}

void SolverHandler::registration()
{
    SolverProxy<EigenCholesky>  SolverChol("EigenCholesky");
    SolverProxy<EigenLU>        SolverLU("EigenLU");
    SolverProxy<EigenUmfPack>   SolverUmfPack("EigenUmfPack");
    SolverProxy<EigenCG>        SolverCG("EigenCG");
    SolverProxy<EigenBiCGSTAB>  SolverBiCGSTAB("EigenBiCGSTAB");
}

void SolverHandler::addProduct(const std::string productName, const SolverBuilder & builder)
{
    M_productList.insert(std::make_pair(productName, builder));
}

SolverPtr_Type SolverHandler::getProduct(const std::string productName) const
{
    std::map<std::string,SolverBuilder>::const_iterator it = M_productList.find(productName);

    return (it == M_productList.end()) ? SolverPtr_Type() : it->second();
}

} // namespace FVCode3D
