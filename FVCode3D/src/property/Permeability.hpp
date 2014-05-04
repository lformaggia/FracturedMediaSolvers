/*!
 *	@file Permeability.hpp
 *	@brief Classes that handle the permeability.
 */

#ifndef PERMEABILITY_HPP_
#define PERMEABILITY_HPP_

#include "core/TypeDefinition.hpp"

namespace Geometry{

//! PermeabilityBase
/*!
 * @class PermeabilityBase
 * This is an abstract class that allows to handle the permeability property.
 * Its derived classes define the permeability as a scalar, diagonal tensor, symmetric tensor and full tensor.
 */
class PermeabilityBase
{
public:

	//! Constructor
	/*!
	 * Create the permeability as a scalar(=1), diagonal tensor(3), symmetric tensor(6) or full tensor(9)
	 * @param size this parameter allows to select the type of the permeability. Default 1.
	 */
	PermeabilityBase(UInt size = 1):M_permeability(size){};

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

	//! Vector-tensor product
	/*!
	 * @param vector vector
	 * @return the product between vector^T and the permeability tensor
	 */
	virtual Point3D vectorTensorProduct(const Point3D & vector) const = 0;

	//! Tensor-vector product
	/*!
	 * @param vector vector
	 * @return the product between the permeability tensor and vector
	 */
	virtual Point3D tensorVectorProduct(const Point3D & vector) const = 0;

	//! Vector-tensor product
	/*!
	 * @param vector vector
	 * @param tensor permeability tensor
	 * @return the product between vector^T and tensor
	 */
	friend Point3D operator*(const Point3D & vector, const PermeabilityBase & tensor);

	//! Tensor-vector product
	/*!
	 * @param tensor permeability tensor
	 * @param vector vector
	 * @return the product between tensor and vector
	 */
	friend Point3D operator*(const PermeabilityBase & tensor, const Point3D & vector);

	//! Destructor
	virtual ~PermeabilityBase(){};

protected:

	//! Vector that represents the permeability tensor
	std::vector<Real> M_permeability;

private:
	//! Size of the permeability
	static const UInt M_size = 0;

	//! No default constructor
	PermeabilityBase();

	//! No copy-constructor
	PermeabilityBase(const PermeabilityBase &);

	//! No assignment operator
	PermeabilityBase & operator=(const PermeabilityBase &);

};

//! Scalar Permeability
/*!
 * @class PermeabilityScalar
 * This class defines the permeability property as a scalar.
 * It's a derived class of PermeabilityBase.
 */
class PermeabilityScalar : public PermeabilityBase
{
public:

	//! Default constructor
	PermeabilityScalar():PermeabilityBase(M_size){};

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

	//! Vector-tensor product
	/*!
	 * @param vector vector
	 * @return the product between vector^T and the permeability tensor
	 */
	virtual Point3D vectorTensorProduct(const Point3D & vector) const;

	//! Tensor-vector product
	/*!
	 * @param vector vector
	 * @return the product between the permeability tensor and vector
	 */
	virtual Point3D tensorVectorProduct(const Point3D & vector) const;

	//! Vector-tensor product
	/*!
	 * @param vector vector
	 * @param tensor permeability tensor
	 * @return the product between vector^T and tensor
	 */
	friend Point3D operator*(const Point3D & vector, const PermeabilityScalar & tensor);

	//! Tensor-vector product
	/*!
	 * @param tensor permeability tensor
	 * @param vector vector
	 * @return the product between tensor and vector
	 */
	friend Point3D operator*(const PermeabilityScalar & tensor, const Point3D & vector);

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
 * It's a derived class of PermeabilityBase.
 */
class PermeabilityDiagonal : public PermeabilityBase
{
public:

	//! Default constructor
	PermeabilityDiagonal():PermeabilityBase(M_size){};

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

	//! Vector-tensor product
	/*!
	 * @param vector vector
	 * @return the product between vector^T and the permeability tensor
	 */
	virtual Point3D vectorTensorProduct(const Point3D & vector) const;

	//! Tensor-vector product
	/*!
	 * @param vector vector
	 * @return the product between the permeability tensor and vector
	 */
	virtual Point3D tensorVectorProduct(const Point3D & vector) const;

	//! Vector-tensor product
	/*!
	 * @param vector vector
	 * @param tensor permeability tensor
	 * @return the product between vector^T and tensor
	 */
	friend Point3D operator*(const Point3D & vector, const PermeabilityDiagonal & tensor);

	//! Tensor-vector product
	/*!
	 * @param tensor permeability tensor
	 * @param vector vector
	 * @return the product between tensor and vector
	 */
	friend Point3D operator*(const PermeabilityDiagonal & tensor, const Point3D & vector);

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
 * It's a derived class of PermeabilityBase.
 */
class PermeabilitySymTensor : public PermeabilityBase
{
public:

	//! Default constructor
	PermeabilitySymTensor():PermeabilityBase(M_size){};

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

	//! Vector-tensor product
	/*!
	 * @param vector vector
	 * @return the product between vector^T and the permeability tensor
	 */
	virtual Point3D vectorTensorProduct(const Point3D & vector) const;

	//! Tensor-vector product
	/*!
	 * @param vector vector
	 * @return the product between the permeability tensor and vector
	 */
	virtual Point3D tensorVectorProduct(const Point3D & vector) const;

	//! Vector-tensor product
	/*!
	 * @param vector vector
	 * @param tensor permeability tensor
	 * @return the product between vector^T and tensor
	 */
	friend Point3D operator*(const Point3D & vector, const PermeabilitySymTensor & tensor);

	//! Tensor-vector product
	/*!
	 * @param tensor permeability tensor
	 * @param vector vector
	 * @return the product between tensor and vector
	 */
	friend Point3D operator*(const PermeabilitySymTensor & tensor, const Point3D & vector);

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
 * It's a derived class of PermeabilityBase.
 */
class PermeabilityFullTensor : public PermeabilityBase
{
public:

	//! Default constructor
	PermeabilityFullTensor():PermeabilityBase(M_size){};

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

	//! Vector-tensor product
	/*!
	 * @param vector vector
	 * @return the product between vector^T and the permeability tensor
	 */
	virtual Point3D vectorTensorProduct(const Point3D & vector) const;

	//! Tensor-vector product
	/*!
	 * @param vector vector
	 * @return the product between the permeability tensor and vector
	 */
	virtual Point3D tensorVectorProduct(const Point3D & vector) const;

	//! Vector-tensor product
	/*!
	 * @param vector vector
	 * @param tensor permeability tensor
	 * @return the product between vector^T and tensor
	 */
	friend Point3D operator*(const Point3D & vector, const PermeabilityFullTensor & tensor);

	//! Tensor-vector product
	/*!
	 * @param tensor permeability tensor
	 * @param vector vector
	 * @return the product between tensor and vector
	 */
	friend Point3D operator*(const PermeabilityFullTensor & tensor, const Point3D & vector);

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

} // namespace Geometry

#endif /* PERMEABILITY_HPP_ */
