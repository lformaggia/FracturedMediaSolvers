/*!
 * @file preconHandler.cpp
 * @brief Factory that handles the creation of the preconditioner classes (definitions).
 */
 
#include <FVCode3D/preconditioner/preconditioner.hpp>
#include <FVCode3D/preconditioner/preconHandler.hpp>
#include <FVCode3D/utility/StringManipolator.hpp>

namespace FVCode3D
{

//! Initializes the static pointer to singleton
std::unique_ptr<preconHandler> preconHandler::theOnlypreconHandler = 0;

preconHandler & preconHandler::Instance()
{
    if (theOnlypreconHandler == 0)
    {
        theOnlypreconHandler.reset(new preconHandler());
        theOnlypreconHandler->registration();
    }

    return *theOnlypreconHandler;
} // preconHandler::Instance

void preconHandler::registration()
{
    preconProxy<identity_preconditioner>          preconId( toUpper( "Identity" ) );
    preconProxy<diagonal_preconditioner>          preconDiag( toUpper( "Diagonal" ) );
    preconProxy<BlockTriangular_preconditioner>   preconBT( toUpper( "BlockTriangular" ) );
    preconProxy<ILU_preconditioner>               preconILU( toUpper( "ILU" ) );
    preconProxy<HSS_preconditioner>               preconHSS( toUpper( "HSS" ) );
    preconProxy<BlockDiagonal_preconditioner>     preconBDiag(toUpper( "BlockDiagonal" ) );
} // preconHandler::registration

void preconHandler::addProduct(const std::string productName, const preconBuilder & builder)
{
    M_productList.insert(std::make_pair(productName, builder));
} // preconHandler::addProduct

preconPtr_Type preconHandler::getProduct(const std::string productName) const
{
    std::map<std::string,preconBuilder>::const_iterator it = M_productList.find( toUpper( productName ) );

    return (it == M_productList.end()) ? preconPtr_Type() : it->second();
} // preconHandler::getProduct

} // namespace FVCode3D
