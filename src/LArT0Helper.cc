/**
 *  @file   LArRecoND/src/LArT0Helper.cc
 *
 *  @brief  Implementation of the t0 (e.g. charge+light matched) helper class.
 *
 *  $Log: $
 */

#include "LArT0Helper.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include <algorithm>
#include <cmath>
#include <limits>

using namespace pandora;

namespace lar_nd_reco
{

//------------------------------------------------------------------------------------------------------------------------------------------

float LArT0Helper::GetMedianTime(const pandora::CaloHitList &inputCaloHitList)
{
  std::vector<float> times;

  // Get calo hit times and sort them lowest to highest
  for ( auto const& calohit : inputCaloHitList)
  {
    times.push_back(calohit->GetTime());
  }
  std::sort(times.begin(), times.end());

  if ( size(times) == 0 )
  {
    std::cout << "Empty calo hit list being used in t0 helper." << std::endl;
    return -1.f;
  }
  else if ( size(times) == 1 )
  {
    return times[0];
  }
  // If an even number of hits, average the values on either side of the middle
  else if ( size(times)%2 == 0 )
  {
    float mid_f = float(size(times))/2.0;
    unsigned int mid = (unsigned int)mid_f;
    return (times[mid]+times[mid-1])/2.;
  }
  // If an odd number of hits, grab the middle entry of the vector
  else
  {
    unsigned int mid = (unsigned int)(size(times)/2)+1;
    return times[mid];
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_nd_reco
