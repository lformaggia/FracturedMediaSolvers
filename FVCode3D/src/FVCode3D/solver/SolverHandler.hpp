#ifndef SOLVERHANDLER_HPP_
#define SOLVERHANDLER_HPP_

#include <FVCode3D/core/TypeDefinition.hpp>
#include <FVCode3D/solver/Solver.hpp>

namespace FVCode3D
{

/*!
 * Builder
 */
typedef SolverPtr_Type (* SolverBuilder)();


/*!
 * Handler
 */
class SolverHandler
{
public:

    static SolverHandler & Instance();

    void addProduct(const std::string productName, const SolverBuilder & builder);

    SolverPtr_Type getProduct(const std::string productName) const;

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

    //! Pointer to the singleton instance
    static SolverHandler * theOnlySolverHandler;
};


/*!
 * Proxy
 */
template<typename T>
class SolverProxy
{
public:

    SolverProxy(const std::string & proxyName);

    ~SolverProxy() = default;

    static SolverPtr_Type Build();

private:

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
