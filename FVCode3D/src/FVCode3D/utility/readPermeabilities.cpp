#include "readPermeabilities.hpp"
#include<limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <algorithm>
//#include <iostream>// only for debugging
namespace
{
  Utility::PermeabilityData readPermeabilityDatum(int zone,std::string const & token, GetPot const & getpotFile)
  {
    constexpr int NZ = std::numeric_limits<int>::min();
    double constexpr INVALID = std::numeric_limits<double>::min();
    int type = getpotFile(token,0,NZ);
    if (type == NZ) throw std::runtime_error(std::string("Token ")+token+std::string(" not present in data file"));
    std::vector<double> values;
    switch(type)
      {
      case 0 :
        values.emplace_back(getpotFile(token,1,INVALID));
        break;
      case 1 :
        values.emplace_back(getpotFile(token,1,INVALID));
        values.emplace_back(getpotFile(token,2,INVALID));
        values.emplace_back(getpotFile(token,3,INVALID));
        break;
      case 2 :
        values.emplace_back(getpotFile(token,1,INVALID));
        values.emplace_back(getpotFile(token,2,INVALID));
        values.emplace_back(getpotFile(token,3,INVALID));
        values.emplace_back(getpotFile(token,4,INVALID));
        values.emplace_back(getpotFile(token,5,INVALID));
        values.emplace_back(getpotFile(token,6,INVALID));
        break;
      }
    for (auto & i:values)
      {
        if (i==INVALID)
          throw std::runtime_error(std::string("Invalid permebility value in data file. Token: ")+token);
      }
    return std::make_tuple(zone,std::make_tuple(type,values));
  }
}// End anonymous namespace
namespace Utility
{
  std::tuple<BulkPermeabilityData,FracturePermeabilityData>
  readPermeabilityData(GetPot const & getpotFile)
  {
    constexpr int NZ = std::numeric_limits<int>::min();
    BulkPermeabilityData bulk;
    FracturePermeabilityData fractures;
    // Read defaults
    bulk.emplace_back(readPermeabilityDatum(0,std::string("problem/bulk/permeability")        ,getpotFile));
    fractures.emplace_back(readPermeabilityDatum(0,std::string("problem/fracture/permeability"),getpotFile));
    // Now the zones
    int numBulkZones = getpotFile("problem/bulk/numzones",0);
    std::string head{"problem/bulk/zone/"};
    for (int z=1; z<=numBulkZones;++z)
      {
        int zonenumber = getpotFile("problem/bulk/zones",z-1,NZ);
        if (zonenumber==NZ) throw std::runtime_error
                              (
                               std::string("You have ")+std::to_string(numBulkZones)+
                               std::string(" zones in the bulk. You need to specify their value using zones=")
                               );
        std::string token=head+std::to_string(zonenumber)+std::string("/permeability");
        bulk.emplace_back(readPermeabilityDatum(zonenumber,token,getpotFile));
      }
    int numFractureZones = getpotFile("problem/fracture/numzones",0);
    head="problem/fracture/zone/";
    for (int z=1; z<=numFractureZones;++z)
      {
        int zonenumber = getpotFile("problem/fracture/zones",z-1,NZ);
        if (zonenumber==NZ) throw std::runtime_error
                              (
                               std::string("You have ")+std::to_string(numFractureZones)+
                               std::string(" zones in the fractures,  you need to specify their value using zones=")
                               );
        std::string token=head+std::to_string(zonenumber)+std::string("/permeability");
        fractures.emplace_back(readPermeabilityDatum(zonenumber,token,getpotFile));
      }
    return std::make_tuple(bulk,fractures);
  }

  //**
  Permeability
  getPermeability(PermeabilityDataList const & permeabilities,int zoneNumber)
  {
      static bool cached{false};
      static Utility::PermeabilityData cachedData=std::make_tuple(0,std::make_tuple(0,std::vector<double>()));
      if (!cached || std::get<0>(cachedData)!=zoneNumber)
      {
          auto it=std::find_if(
                  permeabilities.begin(), permeabilities.end(),
                  [&zoneNumber](const Utility::PermeabilityData & t){return std::get<0>(t)==zoneNumber;}
                  );
          if(it==permeabilities.end())
          {
              throw std::runtime_error(std::string("I cannot find zone Number ")
              +std::to_string(zoneNumber)+std::string(" in __FILE__:__LINE__"));
          }
          else
          {
              cachedData=*it;
              cached=true;
          }
      }
      return std::get<1>(cachedData);
  }


}
