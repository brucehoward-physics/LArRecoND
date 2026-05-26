/**
 *  @file   LArRecoND/include/LArT0Helper.h
 *
 *  @brief  Header file for using t0 information from hits in ND-LAr
 *
 *  $Log: $
 */
#ifndef LAR_T0_HELPER_H
#define LAR_T0_HELPER_H 1

#include "Objects/CaloHit.h"

namespace lar_nd_reco
{

/**
 *  @brief  LArT0Helper class
 */
class LArT0Helper
{
public:
    /**
     *  @brief  Get the time for input calo hit list, here using median time
     *
     *  @param  inputHitList the list of input CaloHits
     *
     *  @return the median time
     */
    static float GetMedianTime(const pandora::CaloHitList &inputCaloHitList);

};

} // namespace lar_nd_reco

#endif // #ifndef LAR_T0_HELPER_H
