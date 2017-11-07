/*!
 * @file PermeabilityFactory.hpp
 * @brief Factory that handles the creation of the Permeability classes.
 */

#ifndef PERMEABILITYFACTORY_HPP_
#define PERMEABILITYFACTORY_HPP_

#include <FVCode3D/core/TypeDefinition.hpp>
#include <FVCode3D/property/Permeability.hpp>

namespace FVCode3D
{

//! Typedef for PermeabilityBuilder function
/*!
 * @typedef PermeabilityBuilder
 * It is the Builder: takes nothing, returns a PermPtr_Type
 */
typedef PermPtr_Type (* PermeabilityBuilder)();

//! Factory that produces Permeability classes
/*!
 * @class PermeabilityFactory
 * This class implements a Factory that produces Permeability classes.
 */
class PermeabilityFactory
{
public:

    //! Static method that return an instance of this class.
    /*!
     * @return PermeabilityFactory instance
     */
    static PermeabilityFactory & Instance();

    //! Add a product to the factory
    /*!
     * @param productName name of the product
     * @param builder builder function
     */
    void addProduct(const std::string productName, const PermeabilityBuilder & builder);

    //! Get a product
    /*!
     * @param productName product name to get
     * @return a pointer to the PemeabilityBase class of the product
     */
    PermPtr_Type getProduct(const std::string productName) const;

    //! Destructor
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


//! Proxy for PermeabilityFactory
/*!
 * @class PermeabilityProxy
 * This class implements the Proxy that actually create a new product
 */
template<typename T>
class PermeabilityProxy
{
public:

    //! Constructor
    /*!
     * @param proxyName name of the product
     */
    PermeabilityProxy(const std::string & proxyName);

    //! Destructor
    ~PermeabilityProxy() = default;

    //! Static function that creates a new product
    /*!
     * @return a pointer to a PermeabilityBase class
     */
    static PermPtr_Type Build();

private:

    //! No default constructor
    PermeabilityProxy() = delete;

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
