/**
 *  @file   src/PfoThreeDHitAssignmentAlgorithm.cc
 *
 *  @brief  
 *
 *  $Log: $
 */

#include "Api/PandoraApi.h"
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "PfoThreeDHitAssignmentAlgorithm.h"

using namespace pandora;

namespace lar_content
{

PfoThreeDHitAssignmentAlgorithm::PfoThreeDHitAssignmentAlgorithm() : m_inputCaloHitList3DName{""}
{
}

pandora::StatusCode PfoThreeDHitAssignmentAlgorithm::Run()
{
    if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    std::cout << "Running 3D hit assignment" << std::endl;

    const CaloHitList *pCaloHits3D{nullptr};
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputCaloHitList3DName, pCaloHits3D));

    CaloHitList availableHits;
    std::map<const CaloHit *, float> availableHitUPos, availableHitVPos, availableHitWPos;
    for (const CaloHit *pCaloHit : (*pCaloHits3D))
    {
        if (!PandoraContentApi::IsAvailable(*this, pCaloHit))
            continue;

        // !!!! TEMP
        std::cout << "3dHitAssignment hit " << pCaloHit->GetMCMatchPDG() << ": " << pCaloHit->GetMCMatchWeight() << std::endl;
        ///////// !!

        availableHits.emplace_back(pCaloHit);
        const CartesianVector pos3D = pCaloHit->GetPositionVector();
        availableHitUPos[pCaloHit] = PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->YZtoU(pos3D.GetY(), pos3D.GetZ());
        availableHitVPos[pCaloHit] = PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->YZtoV(pos3D.GetY(), pos3D.GetZ());
        availableHitWPos[pCaloHit] = PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->YZtoW(pos3D.GetY(), pos3D.GetZ());
    }

    // Maps to keep track of which 3D hit matches to a 2D hit in a given pfo
    std::map<const ParticleFlowObject *, std::string> pfoToClusterListName;
    std::map<const CaloHit *, const ParticleFlowObject *> availableHitToPfoU, availableHitToPfoV, availableHitToPfoW;
    for (unsigned int i = 0; i < m_inputPfoListNames.size(); ++i)
    {
        const std::string pfoListName(m_inputPfoListNames.at(i));
        const PfoList *pPfoList{nullptr};
        PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, pfoListName, pPfoList));
        if (!pPfoList || pPfoList->empty())
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "PfoThreeDHitAssignment: couldn't find pfo list " << pfoListName << std::endl;
            continue;
        }

        for (const ParticleFlowObject *pPfo : (*pPfoList))
        {
            pfoToClusterListName[pPfo] = m_outputClusterListNames.at(i);

            CaloHitList theseCaloHitsU, theseCaloHitsV, theseCaloHitsW;
            LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, theseCaloHitsU);
            LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, theseCaloHitsV);
            LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, theseCaloHitsW);

            for (const CaloHit *const pHit3D : availableHits)
            {
                const CartesianVector pos3D = pHit3D->GetPositionVector();
                for (const CaloHit *const pHit2D : theseCaloHitsU)
                {
                    const CartesianVector pos2D = pHit2D->GetPositionVector();
                    if (std::fabs(pos3D.GetX() - pos2D.GetX()) > std::numeric_limits<float>::epsilon())
                        continue;
                    if (std::fabs(availableHitUPos.at(pHit3D) - pos2D.GetZ()) > std::numeric_limits<float>::epsilon())
                        continue;
                    availableHitToPfoU[pHit3D] = pPfo;
                    break;
                }

                for (const CaloHit *const pHit2D : theseCaloHitsV)
                {
                    const CartesianVector pos2D = pHit2D->GetPositionVector();
                    if (std::fabs(pos3D.GetX() - pos2D.GetX()) > std::numeric_limits<float>::epsilon())
                        continue;
                    if (std::fabs(availableHitVPos.at(pHit3D) - pos2D.GetZ()) > std::numeric_limits<float>::epsilon())
                        continue;
                    availableHitToPfoV[pHit3D] = pPfo;
                    break;
                }

                for (const CaloHit *const pHit2D : theseCaloHitsW)
                {
                    const CartesianVector pos2D = pHit2D->GetPositionVector();
                    if (std::fabs(pos3D.GetX() - pos2D.GetX()) > std::numeric_limits<float>::epsilon())
                        continue;
                    if (std::fabs(availableHitWPos.at(pHit3D) - pos2D.GetZ()) > std::numeric_limits<float>::epsilon())
                        continue;
                    availableHitToPfoW[pHit3D] = pPfo;
                    break;
                }
            }
        }
    }

    CaloHitList threeDHitsMatchedToOnePfo;
    CaloHitList threeDHitsMatchedToMultiPfos;

    std::map<const CaloHit *, PfoVector> hits3DToPfos;
    for (const CaloHit *const pCaloHit3D : availableHits)
    {
        PfoVector matchedPfos;

        // For U we can always add the pfo as it is first
        if (availableHitToPfoU.count(pCaloHit3D))
        {
            if (std::find(matchedPfos.begin(), matchedPfos.end(), availableHitToPfoU.at(pCaloHit3D)) == matchedPfos.end())
                matchedPfos.emplace_back(availableHitToPfoU.at(pCaloHit3D));
        }

        if (availableHitToPfoV.count(pCaloHit3D))
        {
            if (std::find(matchedPfos.begin(), matchedPfos.end(), availableHitToPfoV.at(pCaloHit3D)) == matchedPfos.end())
                matchedPfos.emplace_back(availableHitToPfoV.at(pCaloHit3D));
        }

        if (availableHitToPfoW.count(pCaloHit3D))
        {
            if (std::find(matchedPfos.begin(), matchedPfos.end(), availableHitToPfoW.at(pCaloHit3D)) == matchedPfos.end())
                matchedPfos.emplace_back(availableHitToPfoW.at(pCaloHit3D));
        }

        if (!hits3DToPfos.count(pCaloHit3D))
            hits3DToPfos[pCaloHit3D] = matchedPfos;

        const unsigned int nPfos(matchedPfos.size());
        if (0 == nPfos)
            continue;
        else if (1 == nPfos)
            threeDHitsMatchedToOnePfo.emplace_back(pCaloHit3D);
        else
            threeDHitsMatchedToMultiPfos.emplace_back(pCaloHit3D);
    }

    std::map<const ParticleFlowObject *, CaloHitList> pfoToHits;

    // Deal with the unambiguous 3D hits first
    for (const CaloHit *const pCaloHit3D : threeDHitsMatchedToOnePfo)
    {
        const ParticleFlowObject *const pPfo = hits3DToPfos[pCaloHit3D][0];
        if (!pfoToHits.count(pPfo))
            pfoToHits[pPfo] = CaloHitList();
        pfoToHits[pPfo].emplace_back(pCaloHit3D);
    }

    // If a 3D hit has 2D hits in more than one pfo we have a decision to make
    for (const CaloHit *const pCaloHit3D : threeDHitsMatchedToMultiPfos)
    {
        if (!hits3DToPfos.count(pCaloHit3D))
            continue;

        const PfoVector &pfoList = hits3DToPfos[pCaloHit3D];

        int bestIndex{-999};
        int mostHits{0};
        int biggestPfoIndex{-999};
        float minDist = std::numeric_limits<float>::max();

        for (unsigned int iPfo = 0; iPfo < pfoList.size(); ++iPfo)
        {
            const ParticleFlowObject *const pPfo = pfoList[iPfo];

            const int n2DHits = LArPfoHelper::GetNumberOfTwoDHits(pPfo);
            if (mostHits < n2DHits)
            {
                mostHits = n2DHits;
                biggestPfoIndex = iPfo;
            }

            if (!pfoToHits.count(pPfo))
                continue;

            const CartesianVector hit3DPos = pCaloHit3D->GetPositionVector();
            for (const CaloHit *const pPfoHit : pfoToHits.at(pPfo))
            {
                const float dist = (hit3DPos - pPfoHit->GetPositionVector()).GetMagnitude();
                if (minDist > dist)
                {
                    minDist = dist;
                    bestIndex = iPfo;
                }
            }
        }

        // If no best index was found, use the pfo with the most 2D hits
        if (bestIndex < 0)
            bestIndex = biggestPfoIndex;

        if (!pfoToHits.count(pfoList[bestIndex]))
            pfoToHits[pfoList[bestIndex]] = CaloHitList();
        pfoToHits[pfoList[bestIndex]].emplace_back(pCaloHit3D);
    }

    // Assign the hits
    for (auto const &pfoPair : pfoToHits)
        AddHitsToPfo(pfoPair.first, pfoPair.second, pfoToClusterListName.at(pfoPair.first));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoThreeDHitAssignmentAlgorithm::AddHitsToPfo(const ParticleFlowObject *pPfo, const CaloHitList &hits, const std::string listName) const
{
    ClusterList clusters3D;
    LArPfoHelper::GetThreeDClusterList(pPfo, clusters3D);
    const unsigned int nClusters(clusters3D.size());

    if (0 == nClusters)
    {
        const ClusterList *pClusterList(nullptr);
        std::string clusterListName;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pClusterList, clusterListName));

        PandoraContentApi::Cluster::Parameters clusterParams;
        clusterParams.m_caloHitList.insert(clusterParams.m_caloHitList.end(), hits.begin(), hits.end());

        const Cluster *pCluster3D(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, clusterParams, pCluster3D));
        if (!pCluster3D || !pClusterList || pClusterList->empty())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, listName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pPfo, pCluster3D));
        std::cout << "Pfo " << pPfo << ": created cluster " << pCluster3D << " with " << hits.size() << " hits" << std::endl;
    }
    else if (1 == nClusters)
    {
        const Cluster *const pCluster3D = *(clusters3D.begin());
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pCluster3D, &hits));
        std::cout << "Pfo " << pPfo << ": added " << hits.size() << " hits to existing cluster " << pCluster3D << std::endl;
    }
    else
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode PfoThreeDHitAssignmentAlgorithm::ReadSettings(const pandora::TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitList3DName", m_inputCaloHitList3DName));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputPfoListNames", m_inputPfoListNames));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadVectorOfValues(xmlHandle, "OutputClusterListNames", m_outputClusterListNames));

    if (m_inputPfoListNames.size() != m_outputClusterListNames.size())
    {
        std::cout << "LArPfoThreeDHitAssignment: number of input pfo list names must be the same as the number of output cluster list names"
                  << std::endl;
        return STATUS_CODE_FAILURE;
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
