/*!
 * @file SolverHandler.hpp
 * @brief Factory that handles the creation of the Solver classes.
 */

#ifndef SOLVERHANDLER_HPP_
#define SOLVERHANDLER_HPP_

#include <FVCode3D/core/TypeDefinition.hpp>
#include <FVCode3D/solver/Solver.hpp>

namespace FVCode3D
{

//! Typedef for SolverBuilder function
/*!
 * @typedef SolverBuilder
 * It is the Builder: takes nothing, returns a SolverPtr_Type
 */
typedef SolverPtr_Type (* SolverBuilder)();


//! Factory that produces Solver classes
/*!
 * @class SolverHandler
 * This class implements a Factory that produces Solver classes.
 */
class SolverHandler
{
public:

    //! Static method that return an instance of this class.
    /*!
     * @return SolverHandler instance
     */
    static SolverHandler & Instance();

    //! Add a product to the factory
    /*!
     * @param productName name of the product
     * @param builder builder function
     */
    void addProduct(const std::string productName, const SolverBuilder & builder);

    //! Get a product
    /*!
     * @param productName product name to get
     * @return a pointer to the Solver class of the product
     */
    SolverPtr_Type getProduct(const std::string productName) const;

    //! Destructor
    ~SolverHandler() = default;

private:

    //! No default constructor
    SolverHandler() = default;

    //! No copy constructor
    SolverHandler(const SolverHandler &) = delete;

    //! No assignment operator
    SolverHandler & operator=(const SolverHandler &) = delete;

    //! Register the products
    /*!
     * Fill the map with the products
     */
    void registration();

    //! Map that contains the products
    std::map<std::string,SolverBuilder> M_productList;
};


//! Proxy for SolverHandler
/*!
 * @class SolverProxy
 * This class implements the Proxy that actually create a new product
 */
template<typename T>
class SolverProxy
{
public:

    //! Constructor
    /*!
     * @param proxyName name of the product
     */
    SolverProxy(const std::string & proxyName);

    //! Destructor
    ~SolverProxy() = default;

    //! Static function that creates a new product
    /*!
     * @return a pointer to a Solver class
     */
    static SolverPtr_Type Build();

private:

    //! No default constructor
    SolverProxy() = delete;

    //! No copy constructor
    SolverProxy(const SolverProxy &) = delete;

    //! No assignment operator
    SolverProxy & operator=(const SolverProxy &) = delete;
};


/*!
 * IMPLEMENTATION
 */

template<typename T>
SolverProxy<T>::SolverProxy(const std::string & proxyName)
{
    SolverHandler & theHandler(SolverHandler::Instance());
    theHandler.addProduct(std::string(proxyName),SolverProxy<T>::Build);
}

template<typename T>
SolverPtr_Type SolverProxy<T>::Build()
{
    return SolverPtr_Type(new T);
}

} // namespace FVCode3D

#endif // SOLVERHANDLER_HPP_
