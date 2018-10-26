/*
 * readFunctionsFromGetpot.cpp
 *
 *  Created on: Oct 25, 2018
 *      Author: forma
 */
#include <FVCode3D/utility/readFunctionsFromGetpot.hpp>
#include <string>
#include <vector>
#include <FVCode3D/utility/functionLoader.hpp>
#include <iostream>
#include <stdexcept>
namespace FVCode3D
{
  const std::vector<std::pair<std::string,BorderLabel> >
    ReadFunctionsFromGetpot::bcList
  {
     {"Left",BorderLabel::Left},
     {"Right",BorderLabel::Right},
     {"Front",BorderLabel::Front},
     {"Back",BorderLabel::Back},
     {"Bottom",BorderLabel::Bottom},
     {"Top", BorderLabel::Top}
    };

  const std::vector<std::pair<std::string,BCType> >
   ReadFunctionsFromGetpot::bcType
   {
      {"Neumann",BCType::Neumann},
      {"Dirichlet",BCType::Dirichlet},
    };

  ReadFunctionsFromGetpot::ReadFunctionsFromGetpot(const GetPot &gp):M_gp(gp)
  {
    this->loadLib(M_gp);
  }

  void ReadFunctionsFromGetpot::setGetpotData(const GetPot & gp)
  {
    this->M_gp=gp;
    this->loadLib(this->M_gp);
  }

  void
  ReadFunctionsFromGetpot::loadLib(const std::string & libFileName)
  {
    try
    {
        this->M_functionLoader.load(libFileName);
        this->invalid=false;
    }
    catch (std::runtime_error & e)
    {
        std::cerr<<"Cannot read library "<<libFileName<<" with functions\n";
        std::cerr<<e.what();
        this->invalid=true;
    }
  }

  void
  ReadFunctionsFromGetpot::loadLib(const GetPot & gp)
  {
    std::string token{"functionLibrary"};
    std::string nullLib{"NULL"};
    // get library from getpot data structure
    std::string libFileName=gp(token,nullLib);
    // check if library is indicated
    if(libFileName==nullLib)
      {
        std::cerr<<"Cannot find library in the getpot file\n";
        this->invalid=true;
      }
    else
      // check if library is there!
    this->loadLib(libFileName);
  }


  std::vector<BoundaryConditions::BorderBC>
  ReadFunctionsFromGetpot::getBCForStandardDomain() const
  {
    GetPot const & gp=this->M_gp;
    std::vector<BoundaryConditions::BorderBC> bc;
    for (auto const & t : ReadFunctionsFromGetpot::bcList)
      {
        //BCType bctype;
        std::string nullFun{"NULL"};
        std::string token{std::string("bc/") + t.first + std::string("/type")};
        std::string funtoken{std::string("bc/") + t.first + std::string("/f")};
        auto val = gp(token,noType);
        if (val == noType)
          {
            std::cerr<<"BC type for "<<t.first<< "changed to Dirichlet\n";
            val=0;
          }
        auto theType = static_cast<BCType> (val);
        // read function name
        std::string funName= gp(funtoken,nullFun);

        if (funName==nullFun)
          {
          funName=std::string{"fZero"};
          std::cerr<<"BC fun for "<<t.first<< "changed to fZero\n";
          }

        auto fun = this->M_functionLoader.getFunction(funName);
        bc.emplace_back(t.second,theType,fun);
      }
    return bc;
  }
  
  Func ReadFunctionsFromGetpot::getSourceFunction() const
  {
    GetPot const & gp=this->M_gp;
    std::string token{"problem/sourceFunction"};
    std::string noVal{"NULL"};
    auto funName = gp(token,noVal);
    if (funName==noVal)
      {
        std::cerr<<"No function name for source found. I take fZero\n";
        funName="fZero";
      }
    return this->M_functionLoader.getFunction(funName);
  }

}
