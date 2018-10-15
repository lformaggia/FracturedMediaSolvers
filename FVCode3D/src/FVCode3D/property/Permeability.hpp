/*!
 *  @file Permeability.hpp
 *  @brief Classes that handle the permeability.
 */

#ifndef PERMEABILITY_HPP_
#define PERMEABILITY_HPP_

#include <FVCode3D/core/TypeDefinition.hpp>

namespace FVCode3D
{

class PermeabilityBase;

//! Typedef for PermPtr_Type
/*!
 * @typedef PermPtr_Type
 * This type definition permits to handle a std::shared_ptr<PermeabilityBase> as a PermPtr_Type.
 */
typedef std::shared_ptr<PermeabilityBase> PermPtr_Type;

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
    PermeabilityBase(const UInt size = 1u): M_permeability(size, 0.) {};

    //! Get the permeability at position (i,j) of the 3x3 tensor
    /*!
     * @param i row
     * @param j column
     * @return the permeability at (i,j)
     */
    virtual Real operator()(const UInt i, const UInt j) const = 0;

    //! Get the permeability at position (i) of the vector associated with the tensor, reading by row
    /*!
     * @param i position of the permeability into the vector
     * @return the permeability at i
     */
    virtual Real operator()(const UInt i) const = 0;

    //! Get the permeability at position [i] of the vector
    /*!
     * @param i position of the permeability into the vector
     * @return the permeability at i
     */
    virtual Real operator[](const UInt i) const
        { return M_permeability[i]; }

    //! Set the permeability at position (i,j) of the 3x3 tensor
    /*!
     * @param perm permeability at (i,j)
     * @param i row
     * @param j column
     */
    virtual void setPermeability(const Real perm, const UInt i, const UInt j) = 0;

    //! Get the permeability at position (i) of the vector associated with the tensor, reading by row
    /*!
     * @param perm permeability at i
     * @param i position of the permeability into the vector
     */
    virtual void setPermeability(const Real perm, const UInt i) = 0;

    //! Get the permeability at position [i] of the vector
    /*!
     * @param perm permeability at i
     * @param i position of the permeability into the vector
     */
    virtual Real & operator[](const UInt i)
        { return M_permeability[i]; }

    //! Compute the norm of the permeability tensor
    /*!
     * @return the norm of the tensor
     */
    virtual Real norm() const;

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

    //! Vector-tensor product
    /*!
     * @param vector vector
     * @param tensor permeability tensor
     * @return the product between vector^T and tensor
     */
    friend Point3D operator*(const Point3D & vector, const PermPtr_Type & tensor);

    //! Tensor-vector product
    /*!
     * @param tensor permeability tensor
     * @param vector vector
     * @return the product between tensor and vector
     */
    friend Point3D operator*(const PermPtr_Type & tensor, const Point3D & vector);

    //! Show the permeability tensor
    /*!
     * @param os output stream
     * @return output stream
     */
    virtual std::ostream & showMe(std::ostream & os = std::cout) const = 0;

    //! Show the permeability tensor
    /*!
     * @param os output stream
     * @param perm shared pointer to PermeabilityBase
     * @return output stream
     */
    friend std::ostream & operator<<(std::ostream & os, const PermPtr_Type & tensor);

    //! Destructor
    virtual ~PermeabilityBase() = default;

protected:

    //! Vector that represents the permeability tensor
    std::vector<Real> M_permeability;

private:
    //! Size of the permeability
    static const UInt M_size = 0;

    //! No default constructor
    PermeabilityBase() = delete;

    //! No copy-constructor
    PermeabilityBase(const PermeabilityBase &) = delete;

    //! No assignment operator
    PermeabilityBase & operator=(const PermeabilityBase &) = delete;

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
    PermeabilityScalar(): PermeabilityBase(PermeabilityScalar::M_size) {};

    //! Get the permeability at position (i,j) of the 3x3 tensor
    /*!
     * @param i row
     * @param j column
     * @return the permeability at (i,j)
     */
    Real operator()(const UInt i, const UInt j) const override
        { return i==j ? M_permeability[0] : 0; };

    //! Get the permeability at position (i) of the vector associated with the tensor, reading by row
    /*!
     * @param i position of the permeability into the vector
     * @return the permeability at i
     */
    Real operator()(const UInt i) const override
        { return i%4==0 ? M_permeability[0] : 0; };

    //! Set the permeability at position (i,j) of the 3x3 tensor
    //! Since it is a scalar I admit only i== j. If not we do nothing
    /*!
     * @param perm permeability at (i,j)
     * @param i row
     * @param j column
     */
    void setPermeability(const Real perm, const UInt i, const UInt j) override
        { M_permeability[0] = i==j ? perm : M_permeability[0]; };

    //! Set the permeability at position (i) of the vector associated with the tensor, reading by row
    /*!
     * @param perm permeability at i
     * @param i position of the permeability into the vector
     */
    void setPermeability(const Real perm, const UInt i) override
        { M_permeability[0] = i%4==0 ? perm : M_permeability[0]; };

    //! Vector-tensor product
    /*!
     * @param vector vector
     * @return the product between vector^T and the permeability tensor
     */
    Point3D vectorTensorProduct(const Point3D & vector) const override;

    //! Tensor-vector product
    /*!
     * @param vector vector
     * @return the product between the permeability tensor and vector
     */
    Point3D tensorVectorProduct(const Point3D & vector) const override;

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

    //! Vector-tensor product
    /*!
     * @param vector vector
     * @param tensor permeability tensor
     * @return the product between vector^T and tensor
     */
    friend Point3D operator*(const Point3D & vector, const std::shared_ptr<PermeabilityScalar> & tensor);

    //! Vector-tensor product
    /*!
     * @param vector vector
     * @param tensor permeability tensor
     * @return the product between vector^T and tensor
     */
    friend Point3D operator*(const std::shared_ptr<PermeabilityScalar> & tensor, const Point3D & vector);

    //! Show the permeability tensor
    /*!
     * @param os output stream
     * @return output stream
     */
    std::ostream & showMe(std::ostream & os = std::cout) const override;

    //! Show the permeability tensor
    /*!
     * @param os output stream
     * @param perm shared pointer to PermeabilityBase
     * @return output stream
     */
    friend std::ostream & operator<<(std::ostream & os, const std::shared_ptr<PermeabilityScalar> & perm);

    //! Destructor
    virtual ~PermeabilityScalar() = default;

private:
    //! Size of the permeability
    static const UInt M_size = 1;

public:
    //! Copy-constructor
    PermeabilityScalar(const PermeabilityScalar &) = default;

    //! Assignment operator
    PermeabilityScalar & operator=(const PermeabilityScalar &) = default;
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
    PermeabilityDiagonal():PermeabilityBase(PermeabilityDiagonal::M_size){};

    //! Get the permeability at position (i,j) of the 3x3 tensor
    /*!
     * @param i row
     * @param j column
     * @return the permeability at (i,j)
     */
    Real operator()(const UInt i, const UInt j) const override
        { return i==j ? M_permeability[i] : 0; };

    //! Get the permeability at position (i) of the vector associated with the tensor, reading by row
    /*!
     * @param i position of the permeability into the vector
     * @return the permeability at i
     */
    Real operator()(const UInt i) const override
        { return i%4==0 ? M_permeability[i/4] : 0; };

    //! Set the permeability at position (i,j) of the 3x3 tensor
    /*!
     * @param perm permeability at (i,j)
     * @param i row
     * @param j column
     */
    void setPermeability(const Real perm, const UInt i, const UInt j) override
        { M_permeability[i] = i==j ? perm : M_permeability[i];};

    //! Get the permeability at position (i) of the vector associated with the tensor, reading by row
    /*!
     * @param perm permeability at i
     * @param i position of the permeability into the vector
     */
    void setPermeability(const Real perm, const UInt i) override
        { M_permeability[i/4] = i%4==0 ? perm : M_permeability[i/4];};

    //! Vector-tensor product
    /*!
     * @param vector vector
     * @return the product between vector^T and the permeability tensor
     */
    Point3D vectorTensorProduct(const Point3D & vector) const override;

    //! Tensor-vector product
    /*!
     * @param vector vector
     * @return the product between the permeability tensor and vector
     */
    Point3D tensorVectorProduct(const Point3D & vector) const override;

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

    //! Vector-tensor product
    /*!
     * @param vector vector
     * @param tensor permeability tensor
     * @return the product between vector^T and tensor
     */
    friend Point3D operator*(const Point3D & vector, const std::shared_ptr<PermeabilityDiagonal> & tensor);

    //! Vector-tensor product
    /*!
     * @param vector vector
     * @param tensor permeability tensor
     * @return the product between vector^T and tensor
     */
    friend Point3D operator*(const std::shared_ptr<PermeabilityDiagonal> & tensor, const Point3D & vector);

    //! Show the permeability tensor
    /*!
     * @param os output stream
     * @return output stream
     */
    std::ostream & showMe(std::ostream & os = std::cout) const override;

    //! Show the permeability tensor
    /*!
     * @param os output stream
     * @param perm shared pointer to PermeabilityBase
     * @return output stream
     */
    friend std::ostream & operator<<(std::ostream & os, const std::shared_ptr<PermeabilityDiagonal> & perm);

    //! Destructor
    virtual ~PermeabilityDiagonal() = default;

private:
    //! Size of the permeability
    static const UInt M_size = 3;

public:
    //! No copy-constructor
    PermeabilityDiagonal(const PermeabilityDiagonal &) = default;

    //! No assignment operator
    PermeabilityDiagonal & operator=(const PermeabilityDiagonal &) = default;
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
    Real operator()(const UInt i, const UInt j) const override
        { return j>=i ? M_permeability[3*i+j-i-(i>1)] : operator()(j,i); };

    //! Get the permeability at position (i) of the vector associated with the tensor, reading by row
    /*!
     * @param i position of the permeability into the vector
     * @return the permeability at i
     */
    Real operator()(const UInt i) const override
        { return operator()(i/3, i%3); };

    //! Set the permeability at position (i,j) of the 3x3 tensor
    /*!
     * @param perm permeability at (i,j)
     * @param i row
     * @param j column
     */
    void setPermeability(const Real perm, const UInt i, const UInt j) override
        { M_permeability[3*i+j-i-(i>1)] = j>=i ? perm : M_permeability[3*i+j-i-(i>1)]; };

    //! Get the permeability at position (i) of the vector associated with the tensor, reading by row
    /*!
     * @param perm permeability at i
     * @param i position of the permeability into the vector
     */
    void setPermeability(const Real perm, const UInt i) override
        { setPermeability(perm, i/3, i%3); };

    //! Vector-tensor product
    /*!
     * @param vector vector
     * @return the product between vector^T and the permeability tensor
     */
    Point3D vectorTensorProduct(const Point3D & vector) const override;

    //! Tensor-vector product
    /*!
     * @param vector vector
     * @return the product between the permeability tensor and vector
     */
    Point3D tensorVectorProduct(const Point3D & vector) const override;

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

    //! Vector-tensor product
    /*!
     * @param vector vector
     * @param tensor permeability tensor
     * @return the product between vector^T and tensor
     */
    friend Point3D operator*(const Point3D & vector, const std::shared_ptr<PermeabilitySymTensor> & tensor);

    //! Vector-tensor product
    /*!
     * @param vector vector
     * @param tensor permeability tensor
     * @return the product between vector^T and tensor
     */
    friend Point3D operator*(const std::shared_ptr<PermeabilitySymTensor> & tensor, const Point3D & vector);

    //! Show the permeability tensor
    /*!
     * @param os output stream
     * @return output stream
     */
    std::ostream & showMe(std::ostream & os = std::cout) const override;

    //! Show the permeability tensor
    /*!
     * @param os output stream
     * @param perm shared pointer to PermeabilityBase
     * @return output stream
     */
    friend std::ostream & operator<<(std::ostream & os, const std::shared_ptr<PermeabilitySymTensor> & perm);

    //! Destructor
    virtual ~PermeabilitySymTensor() = default;

private:
    //! Size of the permeability
    static const UInt M_size = 6;
public:
    //! No copy-constructor
    PermeabilitySymTensor(const PermeabilitySymTensor &) = default;

    //! No assignment operator
    PermeabilitySymTensor & operator=(const PermeabilitySymTensor &) = default;
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
    Real operator()(const UInt i, const UInt j) const override
        { return M_permeability[3*i+j];};

    //! Get the permeability at position (i) of the vector associated with the tensor, reading by row
    /*!
     * @param i position of the permeability into the vector
     * @return the permeability at i
     */
    Real operator()(const UInt i) const override
        { return M_permeability[i]; };

    //! Set the permeability at position (i,j) of the 3x3 tensor
    /*!
     * @param perm permeability at (i,j)
     * @param i row
     * @param j column
     */
    void setPermeability(const Real perm, const UInt i, const UInt j) override
        { M_permeability[3*i+j] = perm; };

    //! Get the permeability at position (i) of the vector associated with the tensor, reading by row
    /*!
     * @param perm permeability at i
     * @param i position of the permeability into the vector
     */
    void setPermeability(const Real perm, const UInt i) override
        { M_permeability[i] = perm; };

    //! Vector-tensor product
    /*!
     * @param vector vector
     * @return the product between vector^T and the permeability tensor
     */
    Point3D vectorTensorProduct(const Point3D & vector) const override;

    //! Tensor-vector product
    /*!
     * @param vector vector
     * @return the product between the permeability tensor and vector
     */
    Point3D tensorVectorProduct(const Point3D & vector) const override;

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

    //! Vector-tensor product
    /*!
     * @param vector vector
     * @param tensor permeability tensor
     * @return the product between vector^T and tensor
     */
    friend Point3D operator*(const Point3D & vector, const std::shared_ptr<PermeabilityFullTensor> & tensor);

    //! Vector-tensor product
    /*!
     * @param vector vector
     * @param tensor permeability tensor
     * @return the product between vector^T and tensor
     */
    friend Point3D operator*(const std::shared_ptr<PermeabilityFullTensor> & tensor, const Point3D & vector);

    //! Show the permeability tensor
    /*!
     * @param os output stream
     * @return output stream
     */
    std::ostream & showMe(std::ostream & os = std::cout) const override;

    //! Show the permeability tensor
    /*!
     * @param os output stream
     * @param perm shared pointer to PermeabilityBase
     * @return output stream
     */
    friend std::ostream & operator<<(std::ostream & os, const std::shared_ptr<PermeabilityFullTensor> & perm);

    //! Destructor
    virtual ~PermeabilityFullTensor() = default;

private:
    //! Size of the permeability
    static const UInt M_size = 9;

public:
    //! No copy-constructor
    PermeabilityFullTensor(const PermeabilityFullTensor &) = default;

    //! No assignment operator
    PermeabilityFullTensor & operator=(const PermeabilityFullTensor &) = default;
};

} // namespace FVCode3D

#endif /* PERMEABILITY_HPP_ */
