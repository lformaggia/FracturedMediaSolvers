/*
 *	@file functions.hpp
 *	@brief This file contains the definitions of functions for BCs, source and sink.
 */

#ifndef FUNCTIONS_HPP_
#define FUNCTIONS_HPP_

//Func SourceDomain = [](Geometry::Point3D p){return (p.x()*p.x() + p.y()*p.y() + p.z()*p.z()) < 1;};
//Func SinkDomain = [](Geometry::Point3D p){return (p.x()*p.x() + p.y()*p.y()) < 1;};

Func fZero = [](Geometry::Point3D p){ return 0.; };
Func fOne = [](Geometry::Point3D p){ return 1.; };
Func fMinusOne = [](Geometry::Point3D p){ return -1.; };

/* grid2 */
Func SourceGrid2 = [](Geometry::Point3D p)
	{return 1*( (
					(p.x()-3.15146e+07)*(p.x()-3.15146e+07) +
					(p.y()-1.68972e+07)*(p.y()-1.68972e+07) +
					(p.z()+13355.1    )*(p.z()+13355.1    )
				  ) <=1e6
	 	 	 	);
	};

Func SinkGrid2 = [](Geometry::Point3D p)
	{return -1*( (
					 (p.x()-3.15147e+07)*(p.x()-3.15147e+07) +
					 (p.y()-1.6882e+07 )*(p.y()-1.6882e+07 ) +
					 (p.z()+13355.1    )*(p.z()+13355.1    )
				   ) <=1e6
	 	 	 	 );
	};

Func SSGrid2 = [](Geometry::Point3D p)
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
Func SourcePGrid3 = [](Geometry::Point3D p)
	{return 100*( (
					(p.x()-3.15389e+07)*(p.x()-3.15389e+07) +
					(p.y()-1.68937e+07)*(p.y()-1.68937e+07) +
					(p.z()+14048.1    )*(p.z()+14048.1    )
				  ) <=5e4
				);
	};

Func SinkPGrid3 = [](Geometry::Point3D p)
	{return -100*( (
					 (p.x()-3.14986e+07)*(p.x()-3.14986e+07) +
					 (p.y()-1.68939e+07)*(p.y()-1.68939e+07) +
					 (p.z()+14041.5    )*(p.z()+14041.5    )
				   ) <=5e4
				 );
	};

/* grid3 orthogonal to fractures */
Func SourceOGrid3 = [](Geometry::Point3D p)
	{return 1*( (
					(p.x()-3.15185e7)*(p.x()-3.15185e7) +
					(p.y()-1.6913e7 )*(p.y()-1.6913e7 ) +
					(p.z()+1.80e4   )*(p.z()+1.80e4   )
				 ) <=1e6
			   );
	};

Func SinkOGrid3 = [](Geometry::Point3D p)
	{return -1*( (
					 (p.x()-3.15185e7)*(p.x()-3.15185e7) +
		     	 	 (p.y()-1.6874e7 )*(p.y()-1.6874e7 ) +
		     	 	 (p.z()+1.635e4  )*(p.z()+1.635e4  )
		     	  ) <=1e6
				);
	};

Func SSOGrid3 = [](Geometry::Point3D p)
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

Func SS = fZero;//SSGrid2;

#endif /* FUNCTIONS_HPP_ */
