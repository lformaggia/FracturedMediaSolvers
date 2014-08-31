/*
 * @file PermeabilityFactory.hpp
 */
#ifndef PERMEABILITYFACTORY_HPP_
#define PERMEABILITYFACTORY_HPP_

#include <FVCode3D/core/TypeDefinition.hpp>
#include <FVCode3D/property/Permeability.hpp>

namespace FVCode3D
{

/*!
 * Builder
 */
typedef PermPtr_Type (* PermeabilityBuilder)();


/*!
 * Factory
 */
class PermeabilityFactory
{
public:

    static PermeabilityFactory & Instance();

    void addProduct(const std::string productName, const PermeabilityBuilder & builder);

    PermPtr_Type getProduct(const std::string productName) const;

    ~PermeabilityFactory() = default;

private:

    //! No default constrcutor
    PermeabilityFactory() = default;

    //! No copy constructor
    PermeabilityFactory(const PermeabilityFactory &) = delete;

    //! No assignment operator
    PermeabilityFactory & operator=(const PermeabilityFactory &) = delete;

    //! Register the products
    /*!
     * Fill the map with the products
     */
    void registration();

    //! Map that contains the products
    std::map<std::string,PermeabilityBuilder> M_productList;

    //! Pointer to the singleton instance
    static PermeabilityFactory * theOnlyPermeabilityFactory;
};


/*!
 * Proxy
 */
template<typename T>
class PermeabilityProxy
{
public:

    PermeabilityProxy(const std::string & proxyName);

    ~PermeabilityProxy() = default;

    static PermPtr_Type Build();

private:

    //! No copy constructor
    PermeabilityProxy(const PermeabilityProxy &) = delete;

    //! No assignment operator
    PermeabilityProxy & operator=(const PermeabilityProxy &) = delete;
};


/*!
 * IMPLEMENTATION
 */

template<typename T>
PermeabilityProxy<T>::PermeabilityProxy(const std::string & proxyName)
{
    PermeabilityFactory & theFactory(PermeabilityFactory::Instance());
    theFactory.addProduct(std::string(proxyName),PermeabilityProxy<T>::Build);
}

template<typename T>
PermPtr_Type PermeabilityProxy<T>::Build()
{
    return PermPtr_Type(new T);
}

} // namespace FVCode3D

#endif /* PERMEABILITYFACTORY_HPP_ */
