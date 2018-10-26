#include <FVCode3D/core/TypeDefinition.hpp>
#include <cmath>
FVCode3D::Func SourceTerm = [](FVCode3D::Point3D p)
    {return 15.*( (
                    (p.x()-0.3)*(p.x()-0.3) +
                    (p.y()-0.3)*(p.y()-0.3) +
                    (p.z()-0.8)*(p.z()-0.8)
                          ) <=0.04
                        )
            -15.*( (
                      (p.x()-0.5)*(p.x()-0.5) +
                      (p.y()-0.5)*(p.y()-0.5) +
                      (p.z()-0.3)*(p.z()-0.3)
                    ) <=0.04
                 );
    };


