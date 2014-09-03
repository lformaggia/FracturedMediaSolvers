#include <FVCode3D/property/PermeabilityFactory.hpp>

namespace FVCode3D
{

//! Initializes the static pointer to singleton
PermeabilityFactory * PermeabilityFactory::theOnlyPermeabilityFactory = 0;

PermeabilityFactory & PermeabilityFactory::Instance()
{
    if (theOnlyPermeabilityFactory == 0)
    {
        theOnlyPermeabilityFactory = new PermeabilityFactory();
        theOnlyPermeabilityFactory->registration();
    }

    return *theOnlyPermeabilityFactory;
}

void PermeabilityFactory::registration()
{
    PermeabilityProxy<PermeabilityScalar> PermScal("ScalarPermeability");
    PermeabilityProxy<PermeabilityDiagonal> PermDiag("DiagonalPermeability");
    PermeabilityProxy<PermeabilitySymTensor> PermSymT("SymTensorPermeability");
    PermeabilityProxy<PermeabilityFullTensor> PermFullT("FullTensorPermeability");
}

void PermeabilityFactory::addProduct(const std::string productName, const PermeabilityBuilder & builder)
{
    M_productList.insert(std::make_pair(productName, builder));
}

PermPtr_Type PermeabilityFactory::getProduct(const std::string productName) const
{
    std::map<std::string,PermeabilityBuilder>::const_iterator it = M_productList.find(productName);

    return (it == M_productList.end()) ? PermPtr_Type() : it->second();
}

} // namespace FVCode3D
