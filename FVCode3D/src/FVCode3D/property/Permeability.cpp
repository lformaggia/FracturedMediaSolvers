/*!
 *  @file Permeability.cpp
 *  @brief Classes that handle the permeability (definitions).
 */

#include <FVCode3D/property/Permeability.hpp>

namespace FVCode3D
{

Real PermeabilityBase::norm() const
{
    Real norm(0.);
    std::for_each( std::begin(M_permeability),
                   std::end(M_permeability),
                   [&](Real item)
                   {
                       norm += item*item;
                   }
                 );
    return std::sqrt(norm);
}

Point3D operator*(const Point3D & vector, const PermeabilityBase & tensor)
{
    return tensor.vectorTensorProduct(vector);
}

Point3D operator*(const PermeabilityBase & tensor, const Point3D & vector)
{
    return tensor.tensorVectorProduct(vector);
}

Point3D operator*(const Point3D & vector, const PermPtr_Type & tensor)
{
    return tensor->vectorTensorProduct(vector);
}

Point3D operator*(const PermPtr_Type & tensor, const Point3D & vector)
{
    return tensor->tensorVectorProduct(vector);
}

std::ostream & operator<<(std::ostream & os, const PermPtr_Type & tensor)
{
    return tensor->showMe(os);
}


Point3D PermeabilityScalar::vectorTensorProduct(const Point3D & vector) const
{
    return M_permeability[0] * vector;
}

Point3D PermeabilityScalar::tensorVectorProduct(const Point3D & vector) const
{
    return M_permeability[0] * vector;
}

Point3D operator*(const Point3D & vector, const PermeabilityScalar & tensor)
{
    return tensor.M_permeability[0] * vector;
}

Point3D operator*(const PermeabilityScalar & tensor, const Point3D & vector)
{
    return tensor.M_permeability[0] * vector;
}

Point3D operator*(const Point3D & vector, const std::shared_ptr<PermeabilityScalar> & tensor)
{
    return vector * *tensor;
}

Point3D operator*(const std::shared_ptr<PermeabilityScalar> & tensor, const Point3D & vector)
{
    return *tensor * vector;
}

std::ostream & PermeabilityScalar::showMe(std::ostream & os) const
{
    os << "| " << M_permeability[0] << " 0 0 |" << std::endl;
    os << "| 0 " << M_permeability[0] << " 0 |" << std::endl;
    os << "| 0 0 " << M_permeability[0] << " |" << std::endl;
    return os;
}

std::ostream & operator<<(std::ostream & os, const std::shared_ptr<PermeabilityScalar> & perm)
{
    os << "| " << perm->M_permeability[0] << " 0 0 |" << std::endl;
    os << "| 0 " << perm->M_permeability[0] << " 0 |" << std::endl;
    os << "| 0 0 " << perm->M_permeability[0] << " |" << std::endl;

    return os;
}


Point3D PermeabilityDiagonal::vectorTensorProduct(const Point3D & vector) const
{
    return Point3D( M_permeability[0] * vector.x(),
                    M_permeability[1] * vector.y(),
                    M_permeability[2] * vector.z());
}

Point3D PermeabilityDiagonal::tensorVectorProduct(const Point3D & vector) const
{
    return Point3D( M_permeability[0] * vector.x(),
                    M_permeability[1] * vector.y(),
                    M_permeability[2] * vector.z());
}

Point3D operator*(const Point3D & vector, const PermeabilityDiagonal & tensor)
{
    return Point3D( tensor.M_permeability[0] * vector.x(),
                    tensor.M_permeability[1] * vector.y(),
                    tensor.M_permeability[2] * vector.z());
}

Point3D operator*(const PermeabilityDiagonal & tensor, const Point3D & vector)
{
    return Point3D( tensor.M_permeability[0] * vector.x(),
                    tensor.M_permeability[1] * vector.y(),
                    tensor.M_permeability[2] * vector.z());
}

Point3D operator*(const Point3D & vector, const std::shared_ptr<PermeabilityDiagonal> & tensor)
{
    return vector * *tensor;
}

Point3D operator*(const std::shared_ptr<PermeabilityDiagonal> & tensor, const Point3D & vector)
{
    return *tensor * vector;
}

std::ostream & PermeabilityDiagonal::showMe(std::ostream & os) const
{
    os << "| " << M_permeability[0] << " 0 0 |" << std::endl;
    os << "| 0 " << M_permeability[1] << " 0 |" << std::endl;
    os << "| 0 0 " << M_permeability[2] << " |" << std::endl;
    return os;
}

std::ostream & operator<<(std::ostream & os, const std::shared_ptr<PermeabilityDiagonal> & perm)
{
    os << "| " << perm->M_permeability[0] << " 0 0 |" << std::endl;
    os << "| 0 " << perm->M_permeability[1] << " 0 |" << std::endl;
    os << "| 0 0 " << perm->M_permeability[2] << " |" << std::endl;

    return os;
}


Point3D PermeabilitySymTensor::vectorTensorProduct(const Point3D & vector) const
{
    Real x =    M_permeability[0] * vector.x() +
                M_permeability[1] * vector.y() +
                M_permeability[2] * vector.z();

    Real y =    M_permeability[1] * vector.x() +
                M_permeability[3] * vector.y() +
                M_permeability[4] * vector.z();

    Real z =    M_permeability[2] * vector.x() +
                M_permeability[4] * vector.y() +
                M_permeability[5] * vector.z();

    return Point3D(x, y, z);
}

Point3D PermeabilitySymTensor::tensorVectorProduct(const Point3D & vector) const
{
    Real x =    M_permeability[0] * vector.x() +
                M_permeability[1] * vector.y() +
                M_permeability[2] * vector.z();

    Real y =    M_permeability[1] * vector.x() +
                M_permeability[3] * vector.y() +
                M_permeability[4] * vector.z();

    Real z =    M_permeability[2] * vector.x() +
                M_permeability[4] * vector.y() +
                M_permeability[5] * vector.z();

    return Point3D(x, y, z);
}

Point3D operator*(const Point3D & vector, const PermeabilitySymTensor & tensor)
{
    Real x =    tensor.M_permeability[0] * vector.x() +
                tensor.M_permeability[1] * vector.y() +
                tensor.M_permeability[2] * vector.z();

    Real y =    tensor.M_permeability[1] * vector.x() +
                tensor.M_permeability[3] * vector.y() +
                tensor.M_permeability[4] * vector.z();

    Real z =    tensor.M_permeability[2] * vector.x() +
                tensor.M_permeability[4] * vector.y() +
                tensor.M_permeability[5] * vector.z();

    return Point3D(x, y, z);
}

Point3D operator*(const PermeabilitySymTensor & tensor, const Point3D & vector)
{
    Real x =    tensor.M_permeability[0] * vector.x() +
                tensor.M_permeability[1] * vector.y() +
                tensor.M_permeability[2] * vector.z();

    Real y =    tensor.M_permeability[1] * vector.x() +
                tensor.M_permeability[3] * vector.y() +
                tensor.M_permeability[4] * vector.z();

    Real z =    tensor.M_permeability[2] * vector.x() +
                tensor.M_permeability[4] * vector.y() +
                tensor.M_permeability[5] * vector.z();

    return Point3D(x, y, z);
}

Point3D operator*(const Point3D & vector, const std::shared_ptr<PermeabilitySymTensor> & tensor)
{
    return vector * *tensor;
}

Point3D operator*(const std::shared_ptr<PermeabilitySymTensor> & tensor, const Point3D & vector)
{
    return *tensor * vector;
}

std::ostream & PermeabilitySymTensor::showMe(std::ostream & os) const
{
    os << "| " << M_permeability[0] << " " << M_permeability[1] << " " << M_permeability[2] << " |" << std::endl;
    os << "| " << M_permeability[1] << " " << M_permeability[3] << " " << M_permeability[4] << " |" << std::endl;
    os << "| " << M_permeability[2] << " " << M_permeability[4] << " " << M_permeability[5] << " |" << std::endl;
    return os;
}

std::ostream & operator<<(std::ostream & os, const std::shared_ptr<PermeabilitySymTensor> & perm)
{
    os << "| " << perm->M_permeability[0] << " " << perm->M_permeability[1] << " " << perm->M_permeability[2] << " |" << std::endl;
    os << "| " << perm->M_permeability[1] << " " << perm->M_permeability[3] << " " << perm->M_permeability[4] << " |" << std::endl;
    os << "| " << perm->M_permeability[2] << " " << perm->M_permeability[4] << " " << perm->M_permeability[5] << " |" << std::endl;

    return os;
}


Point3D PermeabilityFullTensor::vectorTensorProduct(const Point3D & vector) const
{
    Real x =    M_permeability[0] * vector.x() +
                M_permeability[3] * vector.y() +
                M_permeability[6] * vector.z();

    Real y =    M_permeability[1] * vector.x() +
                M_permeability[4] * vector.y() +
                M_permeability[7] * vector.z();

    Real z =    M_permeability[2] * vector.x() +
                M_permeability[5] * vector.y() +
                M_permeability[8] * vector.z();

    return Point3D(x, y, z);
}

Point3D PermeabilityFullTensor::tensorVectorProduct(const Point3D & vector) const
{
    Real x =    M_permeability[0] * vector.x() +
                M_permeability[1] * vector.y() +
                M_permeability[2] * vector.z();

    Real y =    M_permeability[3] * vector.x() +
                M_permeability[4] * vector.y() +
                M_permeability[5] * vector.z();

    Real z =    M_permeability[6] * vector.x() +
                M_permeability[7] * vector.y() +
                M_permeability[8] * vector.z();

    return Point3D(x, y, z);
}

Point3D operator*(const Point3D & vector, const PermeabilityFullTensor & tensor)
{
    Real x =    tensor.M_permeability[0] * vector.x() +
                tensor.M_permeability[3] * vector.y() +
                tensor.M_permeability[6] * vector.z();

    Real y =    tensor.M_permeability[1] * vector.x() +
                tensor.M_permeability[4] * vector.y() +
                tensor.M_permeability[7] * vector.z();

    Real z =    tensor.M_permeability[2] * vector.x() +
                tensor.M_permeability[5] * vector.y() +
                tensor.M_permeability[8] * vector.z();

    return Point3D(x, y, z);
}

Point3D operator*(const PermeabilityFullTensor & tensor, const Point3D & vector)
{
    Real x =    tensor.M_permeability[0] * vector.x() +
                tensor.M_permeability[1] * vector.y() +
                tensor.M_permeability[2] * vector.z();

    Real y =    tensor.M_permeability[3] * vector.x() +
                tensor.M_permeability[4] * vector.y() +
                tensor.M_permeability[5] * vector.z();

    Real z =    tensor.M_permeability[6] * vector.x() +
                tensor.M_permeability[7] * vector.y() +
                tensor.M_permeability[8] * vector.z();

    return Point3D(x, y, z);
}

Point3D operator*(const Point3D & vector, const std::shared_ptr<PermeabilityFullTensor> & tensor)
{
    return vector * *tensor;
}

Point3D operator*(const std::shared_ptr<PermeabilityFullTensor> & tensor, const Point3D & vector)
{
    return *tensor * vector;
}

std::ostream &  PermeabilityFullTensor::showMe(std::ostream & os) const
{
    os << "| " << M_permeability[0] << " " << M_permeability[1] << " " << M_permeability[2] << " |" << std::endl;
    os << "| " << M_permeability[3] << " " << M_permeability[4] << " " << M_permeability[5] << " |" << std::endl;
    os << "| " << M_permeability[6] << " " << M_permeability[7] << " " << M_permeability[8] << " |" << std::endl;
    return os;
}

std::ostream & operator<<(std::ostream & os, const std::shared_ptr<PermeabilityFullTensor> & perm)
{
    os << "| " << perm->M_permeability[0] << " " << perm->M_permeability[1] << " " << perm->M_permeability[2] << " |" << std::endl;
    os << "| " << perm->M_permeability[3] << " " << perm->M_permeability[4] << " " << perm->M_permeability[5] << " |" << std::endl;
    os << "| " << perm->M_permeability[6] << " " << perm->M_permeability[7] << " " << perm->M_permeability[8] << " |" << std::endl;

    return os;
}

} // namespace FVCode3D
