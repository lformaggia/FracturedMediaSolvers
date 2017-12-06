/*!
 * @file preconHandler.hpp
 * @brief Factory that handles the creation of the preconditioner classes.
 */

#ifndef PRECONDITIONERHANDLER_HPP_
#define PRECONDITIONERHANDLER_HPP_

#include <FVCode3D/core/TypeDefinition.hpp>
#include <FVCode3D/preconditioner/preconditioner.hpp>

namespace FVCode3D
{

//! Typedef for preconBuilder function
/*!
 * @typedef preconBuilder
 * It is the Builder: takes nothing, returns a preconPtr_Type
 */
typedef preconPtr_Type (* preconBuilder)();


//! Factory that produces preconditioner classes
/*!
 * @class preconHandler
 * This class implements a Factory that produces preconditioner classes.
 */
class preconHandler
{
public:

    //! Static method that return an instance of this class.
    /*!
     * @return preconHandler instance
     */
    static preconHandler & Instance();

    //! Add a product to the factory
    /*!
     * @param productName name of the product
     * @param builder builder function
     */
    void addProduct(const std::string productName, const preconBuilder & builder);

    //! Get a product
    /*!
     * @param productName product name to get
     * @return a pointer to the preconditioner class of the product
     */
    preconPtr_Type getProduct(const std::string productName) const;

    //! Destructor
    ~preconHandler() = default;

private:

    //! private default constructor
    preconHandler() = default;

    //! No copy constructor
    preconHandler(const preconHandler &) = delete;

    //! No assignment operator
    preconHandler & operator=(const preconHandler &) = delete;

    //! Register the products
    /*!
     * Fill the map with the products
     */
    void registration();

    //! Map that contains the products
    std::map<std::string,preconBuilder> M_productList;

    //! Pointer to the singleton instance
    static std::unique_ptr<preconHandler> theOnlypreconHandler;
};


//! Proxy for preconHandler
/*!
 * @class preconProxy
 * This class implements the Proxy that actually create a new product
 */
template<typename T>
class preconProxy
{
public:

    //! Constructor
    /*!
     * @param proxyName name of the product
     */
    preconProxy(const std::string & proxyName);

    //! Destructor
    ~preconProxy() = default;

    //! Static function that creates a new product
    /*!
     * @return a pointer to a preconditioner class
     */
    static preconPtr_Type Build();

private:

    //! No default constructor
    preconProxy() = delete;

    //! No copy constructor
    preconProxy(const preconProxy &) = delete;

    //! No assignment operator
    preconProxy & operator=(const preconProxy &) = delete;
};


/*!
 * IMPLEMENTATION
 */

template<typename T>
preconProxy<T>::preconProxy(const std::string & proxyName)
{
    preconHandler & theHandler(preconHandler::Instance());
    theHandler.addProduct(std::string(proxyName),preconProxy<T>::Build);
}

template<typename T>
preconPtr_Type preconProxy<T>::Build()
{
    return preconPtr_Type(new T);
}

} // namespace FVCode3D

#endif // PRECONDITIONERHANDLER_HPP_
