/*!
 *	@file Permeability.cpp
 *	@brief Classes that handle the permeability (definitions).
 */

#include <FVCode3D/property/Permeability.hpp>

namespace FVCode3D
{

Point3D operator*(const Point3D & vector, const PermeabilityBase & tensor)
{
	return tensor.vectorTensorProduct(vector);
//	if(dynamic_cast<const PermeabilityScalar*>(&tensor))
//		return vector * *dynamic_cast<const PermeabilityScalar*>(&tensor);
//	else if (dynamic_cast<const PermeabilityDiagonal*>(&tensor))
//		return vector * *dynamic_cast<const PermeabilityDiagonal*>(&tensor);
//	else if (dynamic_cast<const PermeabilitySymTensor*>(&tensor))
//		return vector * *dynamic_cast<const PermeabilitySymTensor*>(&tensor);
//	else if (dynamic_cast<const PermeabilityFullTensor*>(&tensor))
//		return vector * *dynamic_cast<const PermeabilityFullTensor*>(&tensor);
//	else
//		return Point3D(nan(""), nan(""), nan(""));
}

Point3D operator*(const PermeabilityBase & tensor, const Point3D & vector)
{
	return tensor.tensorVectorProduct(vector);
//	if(dynamic_cast<const PermeabilityScalar*>(&tensor))
//		return *dynamic_cast<const PermeabilityScalar*>(&tensor) * vector;
//	else if (dynamic_cast<const PermeabilityDiagonal*>(&tensor))
//		return *dynamic_cast<const PermeabilityDiagonal*>(&tensor) * vector;
//	else if (dynamic_cast<const PermeabilitySymTensor*>(&tensor))
//		return *dynamic_cast<const PermeabilitySymTensor*>(&tensor) * vector;
//	else if (dynamic_cast<const PermeabilityFullTensor*>(&tensor))
//		return *dynamic_cast<const PermeabilityFullTensor*>(&tensor) * vector;
//	else
//		return Point3D(nan(""), nan(""), nan(""));
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



Point3D PermeabilityDiagonal::vectorTensorProduct(const Point3D & vector) const
{
	return Point3D(	M_permeability[0] * vector.x(),
					M_permeability[1] * vector.y(),
					M_permeability[2] * vector.z());
}

Point3D PermeabilityDiagonal::tensorVectorProduct(const Point3D & vector) const
{
	return Point3D(	M_permeability[0] * vector.x(),
					M_permeability[1] * vector.y(),
					M_permeability[2] * vector.z());
}

Point3D operator*(const Point3D & vector, const PermeabilityDiagonal & tensor)
{
	return Point3D(	tensor.M_permeability[0] * vector.x(),
					tensor.M_permeability[1] * vector.y(),
					tensor.M_permeability[2] * vector.z());
}

Point3D operator*(const PermeabilityDiagonal & tensor, const Point3D & vector)
{
	return Point3D(	tensor.M_permeability[0] * vector.x(),
					tensor.M_permeability[1] * vector.y(),
					tensor.M_permeability[2] * vector.z());
}



Point3D PermeabilitySymTensor::vectorTensorProduct(const Point3D & vector) const
{
	Real x = 	M_permeability[0] * vector.x() +
				M_permeability[1] * vector.y() +
				M_permeability[2] * vector.z();

	Real y = 	M_permeability[1] * vector.x() +
				M_permeability[3] * vector.y() +
				M_permeability[4] * vector.z();

	Real z = 	M_permeability[2] * vector.x() +
				M_permeability[4] * vector.y() +
				M_permeability[5] * vector.z();

	return Point3D(x, y, z);
}

Point3D PermeabilitySymTensor::tensorVectorProduct(const Point3D & vector) const
{
	Real x = 	M_permeability[0] * vector.x() +
				M_permeability[1] * vector.y() +
				M_permeability[2] * vector.z();

	Real y = 	M_permeability[1] * vector.x() +
				M_permeability[3] * vector.y() +
				M_permeability[4] * vector.z();

	Real z = 	M_permeability[2] * vector.x() +
				M_permeability[4] * vector.y() +
				M_permeability[5] * vector.z();

	return Point3D(x, y, z);
}

Point3D operator*(const Point3D & vector, const PermeabilitySymTensor & tensor)
{
	Real x = 	tensor.M_permeability[0] * vector.x() +
				tensor.M_permeability[1] * vector.y() +
				tensor.M_permeability[2] * vector.z();

	Real y = 	tensor.M_permeability[1] * vector.x() +
				tensor.M_permeability[3] * vector.y() +
				tensor.M_permeability[4] * vector.z();

	Real z = 	tensor.M_permeability[2] * vector.x() +
				tensor.M_permeability[4] * vector.y() +
				tensor.M_permeability[5] * vector.z();

	return Point3D(x, y, z);
}

Point3D operator*(const PermeabilitySymTensor & tensor, const Point3D & vector)
{
	Real x = 	tensor.M_permeability[0] * vector.x() +
				tensor.M_permeability[1] * vector.y() +
				tensor.M_permeability[2] * vector.z();

	Real y = 	tensor.M_permeability[1] * vector.x() +
				tensor.M_permeability[3] * vector.y() +
				tensor.M_permeability[4] * vector.z();

	Real z = 	tensor.M_permeability[2] * vector.x() +
				tensor.M_permeability[4] * vector.y() +
				tensor.M_permeability[5] * vector.z();

	return Point3D(x, y, z);
}



Point3D PermeabilityFullTensor::vectorTensorProduct(const Point3D & vector) const
{
	Real x = 	M_permeability[0] * vector.x() +
				M_permeability[3] * vector.y() +
				M_permeability[6] * vector.z();

	Real y = 	M_permeability[1] * vector.x() +
				M_permeability[4] * vector.y() +
				M_permeability[7] * vector.z();

	Real z = 	M_permeability[2] * vector.x() +
				M_permeability[5] * vector.y() +
				M_permeability[8] * vector.z();

	return Point3D(x, y, z);
}

Point3D PermeabilityFullTensor::tensorVectorProduct(const Point3D & vector) const
{
	Real x = 	M_permeability[0] * vector.x() +
				M_permeability[1] * vector.y() +
				M_permeability[2] * vector.z();

	Real y = 	M_permeability[3] * vector.x() +
				M_permeability[4] * vector.y() +
				M_permeability[5] * vector.z();

	Real z = 	M_permeability[6] * vector.x() +
				M_permeability[7] * vector.y() +
				M_permeability[8] * vector.z();

	return Point3D(x, y, z);
}

Point3D operator*(const Point3D & vector, const PermeabilityFullTensor & tensor)
{
	Real x = 	tensor.M_permeability[0] * vector.x() +
				tensor.M_permeability[3] * vector.y() +
				tensor.M_permeability[6] * vector.z();

	Real y = 	tensor.M_permeability[1] * vector.x() +
				tensor.M_permeability[4] * vector.y() +
				tensor.M_permeability[7] * vector.z();

	Real z = 	tensor.M_permeability[2] * vector.x() +
				tensor.M_permeability[5] * vector.y() +
				tensor.M_permeability[8] * vector.z();

	return Point3D(x, y, z);
}

Point3D operator*(const PermeabilityFullTensor & tensor, const Point3D & vector)
{
	Real x = 	tensor.M_permeability[0] * vector.x() +
				tensor.M_permeability[1] * vector.y() +
				tensor.M_permeability[2] * vector.z();

	Real y = 	tensor.M_permeability[3] * vector.x() +
				tensor.M_permeability[4] * vector.y() +
				tensor.M_permeability[5] * vector.z();

	Real z = 	tensor.M_permeability[6] * vector.x() +
				tensor.M_permeability[7] * vector.y() +
				tensor.M_permeability[8] * vector.z();

	return Point3D(x, y, z);
}

} // namespace FVCode3D
