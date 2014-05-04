/*
 * PermeabilityFactory.hpp
 *
 *  Created on: 04/mag/2014
 *      Author: viskius
 */

#ifndef PERMEABILITYFACTORY_HPP_
#define PERMEABILITYFACTORY_HPP_

#include "core/TypeDefinition.hpp"

class PermeabilityBase;

/*!
 * Builder
 */
typedef std::shared_ptr<PermeabilityBase> (* PermeabilityBuilder)();


/*!
 * Factory
 */
class PermeabilityFactory
{
public:

	inline static PermeabilityFactory & Instance();

	inline void addProduct(const std::string productName, const PermeabilityBuilder & builder);

	inline std::shared_ptr<PermeabilityBase> getProduct(const std::string productName) const;

	inline ~PermeabilityFactory();

private:

	//! No default constrcutor
	inline PermeabilityFactory();

	//! No copy constructor
	PermeabilityFactory(const PermeabilityFactory &);

	//! No assignment operator
	PermeabilityFactory & operator=(const PermeabilityFactory &);

	//! Map that contains the products
	std::map<std::string,PermeabilityBuilder> M_productList;
};


/*!
 * Proxy
 */
template<typename T>
class PermeabilityProxy
{
public:

	PermeabilityProxy(char const * const & proxyName);

	~PermeabilityProxy();

	static std::shared_ptr<PermeabilityBase> Build();

private:

	//! No copy constructor
	PermeabilityProxy(const PermeabilityProxy &);

	//! No assignment operator
	PermeabilityProxy & operator=(const PermeabilityProxy &);
};


/*!
 * IMPLEMENTATION
 */

PermeabilityFactory::PermeabilityFactory()
{}

PermeabilityFactory & PermeabilityFactory::Instance()
{
	static PermeabilityFactory theFactory;
	return theFactory;
}

void PermeabilityFactory::addProduct(const std::string productName, const PermeabilityBuilder & builder)
{
	M_productList.insert(std::make_pair(productName, builder));
}

std::shared_ptr<PermeabilityBase> PermeabilityFactory::getProduct(const std::string productName) const
{
	std::map<std::string,PermeabilityBuilder>::const_iterator it = M_productList.find(productName);

	return (it == M_productList.end()) ? std::shared_ptr<PermeabilityBase>() : it->second();
}

PermeabilityFactory::~PermeabilityFactory()
{}


template<typename T>
PermeabilityProxy<T>::PermeabilityProxy(char const * const & proxyName)
{
	PermeabilityFactory & theFactory(PermeabilityFactory::Instance());
	theFactory.addProduct(std::string(proxyName),PermeabilityProxy<T>::Build);
}

template<typename T>
std::shared_ptr<PermeabilityBase> PermeabilityProxy<T>::Build()
{
	return std::shared_ptr<PermeabilityBase>(new T);
}

template<typename T>
PermeabilityProxy<T>::~PermeabilityProxy()
{}

#endif /* PERMEABILITYFACTORY_HPP_ */
