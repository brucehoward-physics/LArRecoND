/**
 *  @file   include/MasterThreeDAlgorithm.h
 *
 *  @brief  Header file for the master algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_MASTER_THREE_D_ALGORITHM_H
#define LAR_MASTER_THREE_D_ALGORITHM_H 1

#include "Pandora/AlgorithmTool.h"
#include "Pandora/ExternallyConfiguredAlgorithm.h"

#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"
#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include <unordered_map>

// Added by BH to do some time tracking:
// see https://en.cppreference.com/w/cpp/chrono/c/time
#include<ctime>

namespace lar_content
{

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  MasterThreeDAlgorithm class
 */
class MasterThreeDAlgorithm : public MasterAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    MasterThreeDAlgorithm() = default;

protected:
    pandora::StatusCode Run();

    /**
     *  @brief  Run the event slicing procedures, dividing available hits up into distinct 3D regions
     *
     *  @param  volumeIdToHitListMap the volume id to hit list map
     *  @param  sliceVector to receive the populated slice vector
     *
     *  @return whether slicing could be run
     */
    pandora::StatusCode RunSlicing(const VolumeIdToHitListMap &volumeIdToHitListMap, SliceVector &sliceVector) const;

    /**
     *  @brief  Create a pandora worker instance to handle a single LArTPC
     *
     *  @param  larTPC the lar tpc
     *  @param  gapList the gap list
     *  @param  settingsFile the pandora settings file
     *  @param  name the pandora instance name
     *
     *  @return the address of the pandora instance
     */
    const pandora::Pandora *CreateWorkerInstance(const pandora::LArTPC &larTPC, const pandora::DetectorGapList &gapList,
        const std::string &settingsFile, const std::string &name) const;

    /**
     *  @brief  Create a pandora worker instance to handle a number of LArTPCs
     *
     *  @param  larTPCMap the lar tpc map
     *  @param  gapList the gap list
     *  @param  settingsFile the pandora settings file
     *  @param  name the pandora instance name
     *
     *  @return the address of the pandora instance
     */
    const pandora::Pandora *CreateWorkerInstance(const pandora::LArTPCMap &larTPCMap, const pandora::DetectorGapList &gapList,
        const std::string &settingsFile, const std::string &name) const;

    /**
     *  @brief  Initialize pandora worker instances
     */
    pandora::StatusCode InitializeWorkerInstances();

    /**
     *  @brief  Get the mapping from lar tpc volume id to lists of all hits, and truncated hits
     *
     *  @param  volumeIdToHitListMap to receive the populated volume id to hit list map
     *
     *  @return status code
     */
    pandora::StatusCode GetVolumeIdToHitListMap(VolumeIdToHitListMap &volumeIdToHitListMap) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_MASTER_THREE_D_ALGORITHM_H
