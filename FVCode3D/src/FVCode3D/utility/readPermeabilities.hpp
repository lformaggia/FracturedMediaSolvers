#ifndef HH__READPERMEABILITY_HH
#define HH__READPERMEABILITY_HH
#include "GetPot"
#include <vector>
#include <tuple>
namespace Utility
{
  //! Permeability Data
/*!
 * In the format zonenymber type value(s)
 * The number of elements in values depends on type:
 * Scalar (0) -> 1 value
 * Diagonal (1) -> 3 values K11, K22, K33
 * SymTensor (2) -> 6 values K11, K12, K13, K22, K23, K33
 */
  using Permeability = std::tuple<int,std::vector<double>>;
  using PermeabilityData=std::tuple<int,Permeability>;
  //! To store peremabilities, generic name
  using PermeabilityDataList = std::vector<PermeabilityData>;
  //! To store bulk permeabilities
  using BulkPermeabilityData = PermeabilityDataList;
  //! To store fracture permeabilities
  using FracturePermeabilityData = PermeabilityDataList;
  //! Reads permeabilities from getpot object
  //! \par getpotFile A getpot object (not a file...)
  std::tuple<BulkPermeabilityData,FracturePermeabilityData>
  readPermeabilityData(GetPot const & getpotFile);

  //!Extracts a permeability record from the data structure
  //!
  //! @par permeabilities A PeremabilityData structure per zones
  //! @par zoneNumber a zone (0 default)
  //! @return A permeability structure with type and values
  //! @throw std::runtime_error if zonenumber not includes
  Permeability
  getPermeability(PermeabilityDataList const & permeabilities,int zoneNumber);
}
#endif
