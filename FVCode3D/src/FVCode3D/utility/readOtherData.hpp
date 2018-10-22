#ifndef HH_READOTHERDATA_HH
#define HH_READOTHERDATA_HH
#include "GetPot"
#include <tuple>
#include <vector>
namespace Utility
{
  //! Data for zone (0 = default), porosity and aperture
  using fractureData = std::tuple<int,double,double>; // zone, porosity, aperture // @suppress("Type cannot be resolved") // @suppress("Symbol is not resolved")
  //! Data for zone (0 = default) and porosity
  using bulkData = std::tuple<int,double>; // zone, porosity // @suppress("Type cannot be resolved") // @suppress("Symbol is not resolved")
  //! List of data for fracture
  using fractureDataList=std::vector<fractureData>; // @suppress("Invalid template argument")
  //! List of data for bulk
  using bulkDataList=std::vector<bulkData>; // @suppress("Invalid template argument")
  //! Read porosity and aperture for all fracture zones, including default
  fractureDataList readFractureData(GetPot const &);
  //! Read porosity for all bulk zones, including default
  bulkDataList readBulkData(GetPot const &);
  //! Extracts fracture data from a fractureDataList
  //!
  //! @par list a list of data for fractures
  //! @par zoneNumber (a zone, 0 default).
  //! @return a tuple with porosity and aperture
  //!
  //! @throw a runtime_error is zone is not present
  //!
  std::tuple<double,double> getFractureData(fractureDataList const & list, int zoneNumber);
  //! Extracts matrix (bulk) data from a bulkDataList
  //!
  //! @par list a list of data for bulk
  //! @par zoneNumber (a zone, 0 default).
  //! @return the porosity
  //!
  //! @throw a runtime_error is zone is not present
  //!
  double getBulkData(bulkDataList const & list, int zoneNumber);


}
#endif
