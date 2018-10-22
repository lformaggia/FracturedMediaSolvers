
#include <FVCode3D/utility/readOtherData.hpp>
#include <FVCode3D/core/BasicType.hpp>
#include <limits>
#include <string>
namespace Utility
{
  
  fractureDataList readFractureData(GetPot const & gp)
  {
    fractureDataList l;
    constexpr int NZ = std::numeric_limits<int>::min();
    double constexpr INVALID = std::numeric_limits<double>::min();
    std::string head{"problem/fracture/"};
    double porosity = gp(head+std::string{"porosity"},INVALID);
    if (porosity==INVALID)
      throw std::runtime_error(head+std::string{"porosity"}+std::string(" not present in data file"));
    double aperture = gp(head+std::string{"aperture"},INVALID);
    if (aperture==INVALID)
      throw std::runtime_error(head+std::string{"aperture"}+std::string(" not present in data file"));
    l.emplace_back(std::make_tuple(0,porosity,aperture));
    // now the zones, if any
    int numFractureZones = gp("problem/fracture/numzones",0);
    head="problem/fracture/zone/";
    for (int z=1; z<=numFractureZones;++z)
      {
        int zonenumber = gp("problem/fracture/zones",z-1,NZ);
        if (zonenumber==NZ) throw std::runtime_error
                              (
                               std::string("You have ")+std::to_string(numFractureZones)+
                               std::string(" zones in the fractures,  you need to specify their value using zones=")
                               );
        std::string token=head+std::to_string(zonenumber)+std::string("/")+std::string{"porosity"};
        porosity = gp(token,INVALID);
        if (porosity==INVALID)
          throw std::runtime_error(token+std::string(" not present in data file"));
        token=head+std::to_string(zonenumber)+std::string("/")+std::string{"aperture"};
        aperture = gp(token,INVALID);
        if (aperture==INVALID)
          throw std::runtime_error(token+std::string(" not present in data file"));
        l.emplace_back(std::make_tuple(zonenumber,porosity,aperture));
      }
    return l;
  }

  bulkDataList readBulkData(GetPot const & gp)
  {
    bulkDataList l;
    constexpr int NZ = std::numeric_limits<int>::min();
    double constexpr INVALID = std::numeric_limits<double>::min();
    std::string head{"problem/bulk/"};
    std::string token=head + std::string{"porosity"};
    double porosity = gp(token,INVALID);
    if (porosity==INVALID)
      throw std::runtime_error(token+std::string(" not present in data file"));
    l.emplace_back(std::make_tuple(0,porosity));
    // now the zones, if any
    int numBulkZones = gp("problem/bulk/numzones",0);
    head="problem/bulk/zone/";
    for (int z=1; z<=numBulkZones;++z)
      {
        int zonenumber = gp("problem/bulk/zones",z-1,NZ);
        if (zonenumber==NZ) throw std::runtime_error
                              (
                               std::string("You have ")+std::to_string(numBulkZones)+
                               std::string(" zones in the bulk,  you need to specify their value using zones=")
                               );
        std::string token=head+std::to_string(zonenumber)+std::string("/")+std::string{"porosity"};
        porosity = gp(token,INVALID);
        if (porosity==INVALID)
          throw std::runtime_error(token+std::string(" not present in data file"));
        l.emplace_back(std::make_tuple(zonenumber,porosity));
      }
    return l;
  }

  std::tuple<double,double> getFractureData(fractureDataList const & list, int zoneNumber)
  {
      static bool cached{false};
      static std::tuple<int,double,double> cachedData{0,0.,0.};
      if (!cached || std::get<0>(cachedData)!=zoneNumber)
      {
          auto it=std::find_if
                  (
                  list.begin(), list.end(),
                  [&zoneNumber](std::tuple<int,double,double> const & t){return std::get<0>(t)==zoneNumber;}
                  );
          if(it==list.end())
          {
              throw std::runtime_error(std::string("I cannot find zone Number")
              +std::to_string(zoneNumber)+std::string(" __FILE__:__LINE__"));
          }
          else
          {
              cachedData=*it;
              cached=true;
          }
      }
      return std::make_tuple(std::get<1>(cachedData),std::get<2>(cachedData));
  }
  double getBulkData(bulkDataList const & list, int zoneNumber)
  {
      static bool cached{false};
      static std::tuple<int,double> cachedData{0,0.};
      if (!cached || std::get<0>(cachedData)!=zoneNumber)
      {
          auto it=std::find_if(
                  list.begin(), list.end(),
                  [&zoneNumber](std::tuple<int,double> const & t){return std::get<0>(t)==zoneNumber;}
          );
          if(it==list.end())
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
}// end namespace utility
