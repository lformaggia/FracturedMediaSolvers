#ifndef HELPERFUNCTION_HPP
#define HELPERFUNCTION_HPP
#include <iostream>
#include <string>
#include <cmath>
#include <cstddef>
#include "GetPot"
#include <string>
#include "./src_Turconi/TypeDefinition.hpp"


typedef	std::function<double (Geometry::Point2D point)> Func;

//! Helper function to read data from terminal or a source
void readParameters(const int argc, char** argv, double &L_max, Func &u, Func &derxu, Func &deryu, std::ostream & out=std::cout);

auto _f_sinsin = [](Geometry::Point2D p){return -(sin(p.x())*sin(p.y()));};
auto _f_sincos = [](Geometry::Point2D p){return -(sin(p.x())*cos(p.y()));};
auto _f_cossin = [](Geometry::Point2D p){return -(cos(p.x())*sin(p.y()));};
auto f_sinsin = [](Geometry::Point2D p){return (sin(p.x())*sin(p.y()));};
auto f_sincos = [](Geometry::Point2D p){return (sin(p.x())*cos(p.y()));};
auto f_cossin = [](Geometry::Point2D p){return (cos(p.x())*sin(p.y()));};
auto f_coscos = [](Geometry::Point2D p){return (cos(p.x())*cos(p.y()));};



void readParameters(const int argc, char** argv, double &L_max, Func &u, Func &derxu, Func &deryu, std::ostream & out)
{
  GetPot cl(argc, argv);
  //help case
  if(cl.search(2, "--help", "-h"))
  {
    out<<"Compute approximation of a PDE on the domain (0,pi)x(0,pi)."<<"\n";
    out<<"Choosen an exact solution u, between the choices," << std::endl;
	out<<"it is possible to solve -div(grad(x))+x=3u."<<"\n";
    out<<"There are Neumann BC on (0,0)-(0,pi) and on (pi,0)-(pi.pi)." << std::endl;
	out<<"Dirichlet on both other borders"<<"\n";
    out<<"The L2 norm of x-u is returned by the program"<<"\n";
    out<<"Possible arguments:"<<"\n";
    out<<"L_max=Value  (default 0.7) in the range (0.5,pi/2),\n is the maximum lenght of a facet of the mesh"<<"\n";
    out<<"u= exact solution (default fsinsin)"<<"\n";
    out<<"	Possible choices:"<<"\n";
    out<<"	f_sinsin: f(x)=sin(x)sin(y)"<<"\n";
    out<<"	f_coscos: f(x)=cos(x)cos(y)"<<"\n";
    out<<"	f_sincos: f(x)=sin(x)cos(y)"<<"\n";
    out<<"-h or --help : this help"<<"\n";
    std::exit(1);
  }
  //mesh biggest element
  L_max=cl("L_max", 0.5);
  //L_max=0.045+L_max*0.005;

  //define exact solution
  std::string w=cl("u","f_sinsin");
  if(w == "f_sinsin")
  {
	u=f_sinsin;
	derxu=f_cossin;
	deryu=f_sincos;
  }
  if(w == "f_coscos")
  {
	u=f_coscos;
	derxu=_f_sincos;
	deryu=_f_cossin;
  }
  if(w == "f_sincos")
  {
	u=f_sincos;
	derxu=f_coscos;
	deryu=_f_sinsin;
  }
}
#endif
