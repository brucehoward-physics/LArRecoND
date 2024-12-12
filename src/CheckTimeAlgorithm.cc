/**
 *  @file   CheckTimeAlgorithm.cc
 *
 *  @brief  Dump time info
 
 *  $Log: $
 */

// BH: again with a shout-out to Andy C

#include "Pandora/AlgorithmHeaders.h"

#include "CheckTimeAlgorithm.h"

using namespace pandora;

namespace lar_content
{

  CheckTimeAlgorithm::CheckTimeAlgorithm()
  {
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  StatusCode CheckTimeAlgorithm::Run()
  {
    std::cout << "////// Check Time Algorithm //////" << std::endl;
    std::cout << "TIMER_TIMER_TIMER_alg " << std::time(nullptr) << m_auxiliaryInfo << std::endl;

    return STATUS_CODE_SUCCESS;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  StatusCode CheckTimeAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
  {
    // Auxiliary info to print
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "AuxInfo", m_auxiliaryInfo));

    return STATUS_CODE_SUCCESS;
  }

} // namespace lar_content