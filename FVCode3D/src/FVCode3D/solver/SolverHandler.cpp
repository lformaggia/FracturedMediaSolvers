/*!
 * @file SolverHandler.cpp
 * @brief Factory that handles the creation of the Solver classes (definitions).
 */
#include <FVCode3D/solver/SolverHandler.hpp>
#include <FVCode3D/utility/StringManipolator.hpp>

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
    SolverProxy<EigenCholesky>  SolverChol( toUpper( "EigenCholesky" ) );
    SolverProxy<EigenLU>        SolverLU( toUpper( "EigenLU" ) );
#ifdef FVCODE3D_HAS_UMFPACK
    SolverProxy<EigenUmfPack>   SolverUmfPack( toUpper( "EigenUmfPack" ) );
#endif // FVCODE3D_HAS_UMFPACK
    SolverProxy<imlCG>        SolverCG( toUpper( "imlCG" ) );
    SolverProxy<imlBiCGSTAB>    SolverBiCGSTAB( toUpper( "imlBiCGSTAB" ) );
    SolverProxy<imlGMRES>       SolverGMRES( toUpper( "imlGMRES" ) );
#ifdef FVCODE3D_HAS_SAMG
    SolverProxy<SamgSym>        SolverSamgSym( toUpper( "SamgSym" ) );
    SolverProxy<SamgNotSym>     SolverSamgNotSym( toUpper( "SamgNotSym" ) );
#endif // FVCODE3D_HAS_SAMG
} // SolverHandler::registration

void SolverHandler::addProduct(const std::string productName, const SolverBuilder & builder)
{
    M_productList.insert(std::make_pair(productName, builder));
} // SolverHandler::addProduct

SolverPtr_Type SolverHandler::getProduct(const std::string productName) const
{
    std::map<std::string,SolverBuilder>::const_iterator it = M_productList.find( toUpper( productName ) );

    return (it == M_productList.end()) ? SolverPtr_Type() : it->second();
} // SolverHandler::getProduct

} // namespace FVCode3D
