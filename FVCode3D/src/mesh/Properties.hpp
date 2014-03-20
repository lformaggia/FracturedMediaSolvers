/*!
 *	@file Properties.hpp
 *	@brief Classes that handle the properties.
 */

#ifndef PROPERTIES_HPP_
#define PROPERTIES_HPP_

#include "core/TypeDefinition.hpp"

namespace Geometry{

//! Permeability
/*!
 * @class Permeability
 * This is an abstract class that allows to handle the permeability property.
 * Its derived classes define the permeability as a scalar, diagonal tensor, symmetric tensor and full tensor.
 */
class Permeability
{
public:

	//! Constructor
	/*!
	 * Create the permeability as a scalar(=1), diagonal tensor(3), symmetric tensor(6) or full tensor(9)
	 * @param size this parameter allows to select the type of the permeability. Default 1.
	 */
	Permeability(UInt size = 1):M_permeability(size){};

	//! Get the permeability at position (i,j) of the 3x3 tensor
	/*!
	 * @param i row
	 * @param j column
	 * @return the permeability at (i,j)
	 */
	virtual Real operator()(UInt i=0, UInt j=0) const = 0;

	//! Get the permeability at position (i) of the vector associated with the tensor, reading by row
	/*!
	 * @param i position of the permeability into the vector
	 * @return the permeability at i
	 */
	virtual Real operator()(UInt i=0) const = 0;

	//! Set the permeability at position (i,j) of the 3x3 tensor
	/*!
	 * @param perm permeability at (i,j)
	 * @param i row
	 * @param j column
	 */
	virtual void operator()(Real perm, UInt i=0, UInt j=0) = 0;

	//! Get the permeability at position (i) of the vector associated with the tensor, reading by row
	/*!
	 * @param perm permeability at i
	 * @param i position of the permeability into the vector
	 */
	virtual void operator()(Real perm, UInt i=0) = 0;

	//! Size of the permeability
	/*!
	 * @return the size of the permeability: 1,3,6,9
	 */
	virtual UInt size()
		{return M_size;};

	//! Destructor
	virtual ~Permeability(){};

protected:

	//! Vector that represents the permeability tensor
	std::vector<Real> M_permeability;

private:
	//! Size of the permeability
	static const UInt M_size = 0;

	//! No default constructor
	Permeability();

	//! No copy-constructor
	Permeability(const Permeability &);

	//! No assignment operator
	Permeability & operator=(const Permeability &);

};

//! Scalar Permeability
/*!
 * @class PermeabilityScalar
 * This class defines the permeability property as a scalar.
 * It's a derived class of Permeability.
 */
class PermeabilityScalar : public Permeability
{
public:

	//! Default constructor
	PermeabilityScalar():Permeability(M_size){};

	//! Get the permeability at position (i,j) of the 3x3 tensor
	/*!
	 * @param i row
	 * @param j column
	 * @return the permeability at (i,j)
	 */
	virtual Real operator()(UInt i=0, UInt j=0) const
		{ return i==j ? M_permeability[0] : 0; };

	//! Get the permeability at position (i) of the vector associated with the tensor, reading by row
	/*!
	 * @param i position of the permeability into the vector
	 * @return the permeability at i
	 */
	virtual Real operator()(UInt i=0) const
		{ return i%4==0 ? M_permeability[0] : 0; };

	//! Set the permeability at position (i,j) of the 3x3 tensor
	/*!
	 * @param perm permeability at (i,j)
	 * @param i row
	 * @param j column
	 */
	virtual void operator()(Real perm, UInt i=0, UInt j=0)
		{ M_permeability[0] = i==j ? perm : M_permeability[0]; };

	//! Get the permeability at position (i) of the vector associated with the tensor, reading by row
	/*!
	 * @param perm permeability at i
	 * @param i position of the permeability into the vector
	 */
	virtual void operator()(Real perm, UInt i=0)
		{ M_permeability[0] = i%4==0 ? perm : M_permeability[0]; };

	//! Destructor
	virtual ~PermeabilityScalar(){};

private:
	//! Size of the permeability
	static const UInt M_size = 1;

	//! No copy-constructor
	PermeabilityScalar(const PermeabilityScalar &);

	//! No assignment operator
	PermeabilityScalar & operator=(const PermeabilityScalar &);
};

//! Diagonal tensor Permeability
/*!
 * @class PermeabilityDiagonal
 * This class defines the permeability property as a diagonal tensor.
 * It's a derived class of Permeability.
 */
class PermeabilityDiagonal : public Permeability
{
public:

	//! Default constructor
	PermeabilityDiagonal():Permeability(M_size){};

	//! Get the permeability at position (i,j) of the 3x3 tensor
	/*!
	 * @param i row
	 * @param j column
	 * @return the permeability at (i,j)
	 */
	virtual Real operator()(UInt i=0, UInt j=0) const
		{ return i==j ? M_permeability[i] : 0; };

	//! Get the permeability at position (i) of the vector associated with the tensor, reading by row
	/*!
	 * @param i position of the permeability into the vector
	 * @return the permeability at i
	 */
	virtual Real operator()(UInt i=0) const
		{ return i%4==0 ? M_permeability[i/4] : 0; };

	//! Set the permeability at position (i,j) of the 3x3 tensor
	/*!
	 * @param perm permeability at (i,j)
	 * @param i row
	 * @param j column
	 */
	virtual void operator()(Real perm, UInt i=0, UInt j=0)
		{ M_permeability[i] = i==j ? perm : M_permeability[i];};

	//! Get the permeability at position (i) of the vector associated with the tensor, reading by row
	/*!
	 * @param perm permeability at i
	 * @param i position of the permeability into the vector
	 */
	virtual void operator()(Real perm, UInt i=0)
		{ M_permeability[i/4] = i%4==0 ? perm : M_permeability[i/4];};

	//! Destructor
	virtual ~PermeabilityDiagonal(){};

private:
	//! Size of the permeability
	static const UInt M_size = 3;

	//! No copy-constructor
	PermeabilityDiagonal(const PermeabilityDiagonal &);

	//! No assignment operator
	PermeabilityDiagonal & operator=(const PermeabilityDiagonal &);
};

//! Symmetric tensor Permeability
/*!
 * @class PermeabilitySymTensor
 * This class defines the permeability property as a symmetric tensor.
 * It's a derived class of Permeability.
 */
class PermeabilitySymTensor : public Permeability
{
public:

	//! Default constructor
	PermeabilitySymTensor():Permeability(M_size){};

	//! Get the permeability at position (i,j) of the 3x3 tensor
	/*!
	 * @param i row
	 * @param j column
	 * @return the permeability at (i,j)
	 */
	virtual Real operator()(UInt i=0, UInt j=0) const
		{ return j>=i ? M_permeability[3*i+j-i-(i>1)] : operator()(j,i); };

	//! Get the permeability at position (i) of the vector associated with the tensor, reading by row
	/*!
	 * @param i position of the permeability into the vector
	 * @return the permeability at i
	 */
	virtual Real operator()(UInt i=0) const
		{ return operator()(i/3, i%3); };

	//! Set the permeability at position (i,j) of the 3x3 tensor
	/*!
	 * @param perm permeability at (i,j)
	 * @param i row
	 * @param j column
	 */
	virtual void operator()(Real perm, UInt i=0, UInt j=0)
		{ M_permeability[3*i+j-i-(i>1)] = j>=i? perm : M_permeability[3*i+j-i-(i>1)]; };

	//! Get the permeability at position (i) of the vector associated with the tensor, reading by row
	/*!
	 * @param perm permeability at i
	 * @param i position of the permeability into the vector
	 */
	virtual void operator()(Real perm, UInt i=0)
		{ operator()(perm, i/3, i%3); };

	//! Destructor
	virtual ~PermeabilitySymTensor(){};

private:
	//! Size of the permeability
	static const UInt M_size = 6;

	//! No copy-constructor
	PermeabilitySymTensor(const PermeabilitySymTensor &);

	//! No assignment operator
	PermeabilitySymTensor & operator=(const PermeabilitySymTensor &);
};

//! Full tensor Permeability
/*!
 * @class PermeabilityFullTensor
 * This class defines the permeability property as a full tensor.
 * It's a derived class of Permeability.
 */
class PermeabilityFullTensor : public Permeability
{
public:

	//! Default constructor
	PermeabilityFullTensor():Permeability(M_size){};

	//! Get the permeability at position (i,j) of the 3x3 tensor
	/*!
	 * @param i row
	 * @param j column
	 * @return the permeability at (i,j)
	 */
	virtual Real operator()(UInt i=0, UInt j=0) const
		{ return M_permeability[3*i+j];};

	//! Get the permeability at position (i) of the vector associated with the tensor, reading by row
	/*!
	 * @param i position of the permeability into the vector
	 * @return the permeability at i
	 */
	virtual Real operator()(UInt i=0) const
		{ return M_permeability[i]; };

	//! Set the permeability at position (i,j) of the 3x3 tensor
	/*!
	 * @param perm permeability at (i,j)
	 * @param i row
	 * @param j column
	 */
	virtual void operator()(Real perm, UInt i=0, UInt j=0)
		{ M_permeability[3*i+j] = perm; };

	//! Get the permeability at position (i) of the vector associated with the tensor, reading by row
	/*!
	 * @param perm permeability at i
	 * @param i position of the permeability into the vector
	 */
	virtual void operator()(Real perm, UInt i=0)
		{ M_permeability[i] = perm; };

	//! Destructor
	virtual ~PermeabilityFullTensor(){};

private:
	//! Size of the permeability
	static const UInt M_size = 9;

	//! No copy-constructor
	PermeabilityFullTensor(const PermeabilityFullTensor &);

	//! No assignment operator
	PermeabilityFullTensor & operator=(const PermeabilityFullTensor &);
};

//! Properties (aperture, porosity, permeability) of a single fracture/matrix.
/*!
 * @class Properties
 * This class contains the properties for a specific zone.
 */
class Properties
{
public:

	//! Empty constructor
	Properties():
		M_aperture(0.), M_porosity(0.), M_permeability(0.) {}

	//! Constructor
	/*!
	 * @param aperture Aperture of the fracture. If cell, should be = 1
	 * @param porosity Porosity of the matrix. If fracture, should be = 1
	 * @param permeability Permeability of the fracture/matrix
	 */
	Properties(const Real aperture, const Real porosity, const Real permeability):
		M_aperture(aperture), M_porosity(porosity), M_permeability(permeability) {}

	//! Set properties
	/*!
	 * Set properties of a single fracture/matrix.
	 * @param aperture Aperture of the fracture. If cell, should be = 1
	 * @param porosity Porosity of the matrix. If fracture, should be = 1
	 * @param permeability Permeability of the fracture/matrix
	 */
	void setProperties(const Real aperture, const Real porosity, const Real permeability)
		{ M_aperture = aperture; M_porosity = porosity;
		  M_permeability = permeability; };

	//! Aperture of the fracture. If cell, should be = 1
	Real M_aperture;
	//! Porosity of the matrix. If fracture, should be = 1
	Real M_porosity;
	//! Permeability of the fracture/matrix
	Real M_permeability;
};

//! Class that handles the properties for all the reservoir
/*!
 * @class PropertiesMap
 * This class implements a map that contains the properties of all zones
 */
class PropertiesMap
{
public:

	//! Empty constructor
	PropertiesMap():M_mobility(0.), M_compressibility(0.) {};

	//! Constructor
	/*!
	 * @param mobility mobility of the fluid
	 */
	PropertiesMap(const Real mobility, const Real compressibility = 0.):
		M_mobility(mobility), M_compressibility(compressibility) {};

	//! Get the number of zones
	/*!
	 * @return the number of zones
	 */
	UInt getNumberOfZone() const
		{ return M_properties.size(); };

	//! Get the mobility of the fluid
	/*!
	 * @return the mobility of the fluid
	 */
	Real getMobility() const
		{ return M_mobility; };

	//! Get the compressbility
	/*!
	 * @return the compressibility
	 */
	Real getCompressibility() const
		{ return M_compressibility; };

	//! Get the map of properties (const)
	/*!
	 * @return the map of the properties. First -> zone code; Second -> Geometry::Properties
	 */
	const std::map<UInt, Properties> & getProperties() const
		{ return M_properties; };

	//! Get the map of properties
	/*!
	 * @return the map of the properties. First -> zone code; Second -> Geometry::Properties
	 */
	std::map<UInt, Properties> & getProperties()
		{ return M_properties; };

	//! Get the properties of a selected zone (const)
	/*!
	 * @param zone the id of the zone code
	 * @return a reference of the related Geometry::Properties
	 */
	const Properties & getProperties(const UInt zone) const
		{ return M_properties.at(zone); };

	//! Get the properties of a selected zone
	/*!
	 * @param zone the id of the zone code
	 * @return a reference of the related Geometry::Properties
	 */
	Properties & getProperties(const UInt zone)
		{ return M_properties[zone]; };

	//! Set the properties of a selected zone
	/*!
	 * @param zone the id of the zone code
	 * @param prop a reference to a Geometry::Properties
	 */
	void setZone(const UInt zone, const Properties & prop)
		{ M_properties[zone] = prop; };

	//! Set the mobility of the fluid
	/*!
	 * @param the mobility of the fluid
	 */
	void setMobility(const Real mobility)
		{ M_mobility = mobility; };

	//! Set the compressibility
	/*!
	 * @param the compressibility
	 */
	void setCompressibility(const Real compressibility)
		{ M_compressibility = compressibility; };

	//! Destructor
	~PropertiesMap() {};

private:

	//! Map that link each zone code with the property
	std::map<UInt, Properties> M_properties;

	//! Mobility of the fluid.
	Real M_mobility;

	//! Compressibility.
	Real M_compressibility;
};

} //namespace Geometry

#endif /* PROPERTIES_HPP_ */
