/*
 *  @file functions.hpp
 *  @brief This file contains the definitions of functions for BCs, source and sink.
 */

#ifndef FUNCTIONS_HPP_
#define FUNCTIONS_HPP_
#include <FVCode3D/core/TypeDefinition.hpp>
#include <cmath>

namespace FVCode3D
{

Func SourceDomain = [](Point3D p){return ((p.x() - 1.)*(p.x() - 1.) + (p.y() - 0.5)*(p.y() - 0.5) + (p.z() - 0.5)*(p.z()
- 0.5)) < 1;};
//Func SinkDomain = [](Point3D p){return (p.x()*p.x() + p.y()*p.y()) < 1;};

Func fZero = [](Point3D){ return 0.; };
Func fOne = [](Point3D){ return 1.; };
Func fTwo = [](Point3D){ return 2.; };
Func fFour = [](Point3D){ return 4.; };
Func fMinusTwo = [](Point3D){ return -2.; };
Func fMinusOne = [](Point3D){ return -1.; };
Func fOneZero = [](Point3D p){ return (2. - p.x()) / 2.; };
Func fTen = [](Point3D){ return 10.; };

Func fZ = [](Point3D p){ return p.z(); };

/* grid2 */
Func SourceGrid2 = [](Point3D p)
    {return 1*( (
                    (p.x()-3.15146e+07)*(p.x()-3.15146e+07) +
                    (p.y()-1.68972e+07)*(p.y()-1.68972e+07) +
                    (p.z()+13355.1    )*(p.z()+13355.1    )
                  ) <=1e6
                );
    };

Func SinkGrid2 = [](Point3D p)
    {return -1*( (
                     (p.x()-3.15147e+07)*(p.x()-3.15147e+07) +
                     (p.y()-1.6882e+07 )*(p.y()-1.6882e+07 ) +
                     (p.z()+13355.1    )*(p.z()+13355.1    )
                   ) <=1e6
                 );
    };

Func SSGrid2 = [](Point3D p)
    {return  1*( (
                   (p.x()-3.15146e+07)*(p.x()-3.15146e+07) +
                   (p.y()-1.68972e+07)*(p.y()-1.68972e+07) +
                   (p.z()+13355.1    )*(p.z()+13355.1    )
                 ) <=1e6
               )
            -1*( (
                   (p.x()-3.15147e+07)*(p.x()-3.15147e+07) +
                   (p.y()-1.6882e+07 )*(p.y()-1.6882e+07 ) +
                   (p.z()+13355.1    )*(p.z()+13355.1    )
                 ) <=1e6
               );
    };

/* grid3 parallel to fractures */
Func SourcePGrid3 = [](Point3D p)
    {return 100*( (
                    (p.x()-3.15389e+07)*(p.x()-3.15389e+07) +
                    (p.y()-1.68937e+07)*(p.y()-1.68937e+07) +
                    (p.z()+14048.1    )*(p.z()+14048.1    )
                  ) <=5e4
                );
    };

Func SinkPGrid3 = [](Point3D p)
    {return -100*( (
                     (p.x()-3.14986e+07)*(p.x()-3.14986e+07) +
                     (p.y()-1.68939e+07)*(p.y()-1.68939e+07) +
                     (p.z()+14041.5    )*(p.z()+14041.5    )
                   ) <=5e4
                 );
    };

/* grid3 orthogonal to fractures */
Func SourceOGrid3 = [](Point3D p)
    {return 1*( (
                    (p.x()-3.15185e7)*(p.x()-3.15185e7) +
                    (p.y()-1.6913e7 )*(p.y()-1.6913e7 ) +
                    (p.z()+1.80e4   )*(p.z()+1.80e4   )
                 ) <=1e6
               );
    };

Func SinkOGrid3 = [](Point3D p)
    {return -1*( (
                     (p.x()-3.15185e7)*(p.x()-3.15185e7) +
                     (p.y()-1.6874e7 )*(p.y()-1.6874e7 ) +
                     (p.z()+1.635e4  )*(p.z()+1.635e4  )
                  ) <=1e6
                );
    };

Func SSOGrid3 = [](Point3D p)
    {return 1*( (
                    (p.x()-3.15185e7)*(p.x()-3.15185e7) +
                    (p.y()-1.6913e7 )*(p.y()-1.6913e7 ) +
                    (p.z()+1.80e4   )*(p.z()+1.80e4   )
                 ) <=1e6
               )
            -1*( (
                    (p.x()-3.15185e7)*(p.x()-3.15185e7) +
                    (p.y()-1.6874e7 )*(p.y()-1.6874e7 ) +
                    (p.z()+1.635e4  )*(p.z()+1.635e4  )
                 ) <=1e6
               );
    };
    
Func SorgentSink = [](Point3D p)
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

/* TEST EDFM */
Func SSEDFM = [](Point3D p)
    { return 10. * ( p.x() <= 0.4 ); };
    
Func SSEDFMBC = [](Point3D p)
    { return 2. * ( (p.x() <= 0.2) && (p.x() >= -0.2 ) ); };

Func SS = SorgentSink;//SourceDomain; //fZero;//SSGrid2;

// TEST CASE: convergence order

Func SourceBulk = [](Point3D p)
	{ return (1-1e3)*cosh(1e-2/2.)*cos(p.x())*cos(p.z()); };
	
Func SourceFrac = [](Point3D p)
	{ return pow(1e3,2)*cos(p.x()) + 1e3*(1-1e3)*cosh(1e-2/2.)*cos(p.x())*cos(p.z()); };

} // namespace FVCode3D
#endif /* FUNCTIONS_HPP_ */
