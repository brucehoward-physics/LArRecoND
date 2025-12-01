/**
 *  @file   LArRecoND/test/PandoraOuterface.cc
 *
 *  @brief  Implementation of the Post-Pandora high-level reco for DUNE ND
 *
 *  $Log: $
 */

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TSpline.h"
#include "TGraph.h"

#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoShape.h"
#include "TGeoVolume.h"

#include "Api/PandoraApi.h"
#include "Geometry/LArTPC.h"
#include "Helpers/XmlHelper.h"
#include "Managers/GeometryManager.h"
#include "Managers/PluginManager.h"
#include "Xml/tinyxml.h"

#include "larpandoracontent/LArContent.h"
#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"
#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"
#include "larpandoracontent/LArPlugins/LArPseudoLayerPlugin.h"
#include "larpandoracontent/LArPlugins/LArRotationalTransformationPlugin.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"

#ifdef LIBTORCH_DL
#include "larpandoradlcontent/LArDLContent.h"
#endif

#include "LArNDContent.h"
#include "LArNDGeomSimple.h"
#include "LArRay.h"
#include "PandoraOuterface.h"

#ifdef MONITORING
#include "TApplication.h"
#endif

#include <algorithm>
#include <cmath>
#include <getopt.h>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>
#include <map>
#include <functional>

using namespace pandora;
using namespace lar_nd_postreco;

int main(int argc, char *argv[])
{

  int errorNo(0);

  try
  {
    ParameterStruct pset;

    if (!ParseCommandLine(argc, argv, pset))
      return 1;

    if (!ReadSettings(pset))
      return 1;

    ProcessPostReco(pset);
  }
  catch (const StatusCodeException &statusCodeException)
  {
    std::cerr << "Pandora StatusCodeException: " << statusCodeException.ToString() << statusCodeException.GetBackTrace() << std::endl;
    errorNo = 1;
  }
  catch (...)
  {
    std::cerr << "Unknown exception: " << std::endl;
    errorNo = 1;
  }

  return errorNo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_nd_postreco
{

float LifetimeCorrectionFactor(const std::vector<float> &detAnodes, const float inputPos, const float lifetime, const float driftSpeed )
{
  float driftDist = std::numeric_limits<float>::max();
  for ( float wallX : detAnodes ) {
    float thisDist = fabs(wallX-inputPos);
    if ( thisDist < driftDist ) driftDist = thisDist;
  }
  float tDrift = driftDist / driftSpeed;
  return TMath::Exp(tDrift/lifetime);
}

float eVisWithRecombination(const ParameterStruct &parameters, const float inputQ)
{
  const float wIon = 23.6/1.0e6; // MeV/e-, a hard-coded for now value used later

  // MIP Recombination with the Q->E calculation as in FLOW file
  // see e.g. https://github.com/DUNE/ndlar_flow/blob/develop/src/proto_nd_flow/reco/charge/calib_prompt_hits.py#L289
  float dEdxMIP = 2.; // MeV/cm (using the value in above line)
  float recomb = 1.;
  if ( parameters.fBoxRecombination ) {
    float csi = parameters.fBoxBeta * dEdxMIP / (parameters.fEField * parameters.fDensity);
    recomb = TMath::Log(parameters.fBoxAlpha + csi)/csi;
  }
  else if ( parameters.fBirksRecombination ) {
    recomb = parameters.fBirksA / ( 1. + parameters.fBirksK * dEdxMIP / (parameters.fEField * parameters.fDensity) );
  }

  return inputQ * wIon / recomb;
}

float dEdxWithRecombination(const ParameterStruct &parameters, const float inputdQdx)
{
  float dEdx_val = 0.;
  const float wIon = 23.6/1.0e6; // MeV/e-, a hard-coded for now value used later

  if ( parameters.fShouldCorrectRecomb ) {
    if ( parameters.fFlowStyleRecombination ) {
      // MIP Recombination with the Q->E calculation as in FLOW file
      // see e.g. https://github.com/DUNE/ndlar_flow/blob/develop/src/proto_nd_flow/reco/charge/calib_prompt_hits.py#L289
      float dEdxMIP = 2.; // MeV/cm (using the value in above line)
      float recomb = 1.;
      if ( parameters.fBoxRecombination ) {
	float csi = parameters.fBoxBeta * dEdxMIP / (parameters.fEField * parameters.fDensity);
	recomb = TMath::Log(parameters.fBoxAlpha + csi)/csi;
      }
      else if ( parameters.fBirksRecombination ) {
	recomb = parameters.fBirksA / ( 1. + parameters.fBirksK * dEdxMIP / (parameters.fEField * parameters.fDensity) );
      }
      dEdx_val = inputdQdx / recomb * wIon;
    }
    else if ( parameters.fBoxRecombination ) {
      // Box style, LArSoft style, angular part turned off, we'll just use the box beta as-is
      // https://github.com/LArSoft/larreco/blob/develop/larreco/Calorimetry/CalorimetryAlg.cxx
      dEdx_val = (TMath::Exp(parameters.fBoxBeta * wIon * inputdQdx) - parameters.fBoxAlpha) / parameters.fBoxBeta;
    }
    else if ( parameters.fBirksRecombination ){
      // Birks style, LArSoft style
      // https://github.com/LArSoft/larreco/blob/develop/larreco/Calorimetry/CalorimetryAlg.cxx
      dEdx_val = inputdQdx / (parameters.fBirksA / wIon - parameters.fBirksK / parameters.fEField * inputdQdx);
    }
  }
  else {
    // Assume recomb = 1?
    float recomb = 1.;
    dEdx_val = inputdQdx / recomb * wIon;
  }

  if ( parameters.fApplyCalibrationFudgeFactor ) dEdx_val*=parameters.fCalibrationFudgeFactor;
  return dEdx_val;
}

float pFromRange_proton( const float inputRange )
{
  /// Set up the necessary pieces for proton momentum vs range spline: CSDA
  /// Result from LArSoft -> LArReco -> TrackMomentumCalculator (v09_26_02)
  float KE = 0.;
  if (inputRange > 0 && inputRange <= 80)
    KE = 29.9317 * std::pow(inputRange, 0.586304);
  else if (inputRange > 80 && inputRange <= 3.022E3)
    KE = 149.904 + (3.34146 * inputRange) + (-0.00318856 * inputRange * inputRange) +
      (4.34587E-6 * inputRange * inputRange * inputRange) +
      (-3.18146E-9 * inputRange * inputRange * inputRange * inputRange) +
      (1.17854E-12 * inputRange * inputRange * inputRange * inputRange * inputRange) +
      (-1.71763E-16 * inputRange * inputRange * inputRange * inputRange * inputRange * inputRange);
  else
    KE = -999;

  // convert KE to Momentum
  float massProton = 938.272;

  if ( KE < 0 )
    return 0.f;
  return std::sqrt((KE * KE) + (2 * massProton * KE))/1000.;
}

float KEFromRange_proton( const float inputRange )
{
  // Same as above but just KE
  float KE = 0.;
  if (inputRange > 0 && inputRange <= 80)
    KE = 29.9317 * std::pow(inputRange, 0.586304);
  else if (inputRange > 80 && inputRange <= 3.022E3)
    KE = 149.904 + (3.34146 * inputRange) + (-0.00318856 * inputRange * inputRange) +
      (4.34587E-6 * inputRange * inputRange * inputRange) +
      (-3.18146E-9 * inputRange * inputRange * inputRange * inputRange) +
      (1.17854E-12 * inputRange * inputRange * inputRange * inputRange * inputRange) +
      (-1.71763E-16 * inputRange * inputRange * inputRange * inputRange * inputRange * inputRange);
  else
    KE = -999;

  if ( KE < 0 )
    return 0.f;
  return KE/1000.;
}

void ProcessPostReco(const ParameterStruct &parameters)
{
  std::vector<float> posAnodes;
  std::vector<float> xBoundaries;
  std::vector<float> yBoundaries;
  std::vector<float> zBoundaries;
  if ( parameters.fDetector == 0 ) {
    // NDLAr anodes
    posAnodes.push_back(-50.);
    posAnodes.push_back(-150.);
    posAnodes.push_back(-250.);
    posAnodes.push_back(-350.);
    posAnodes.push_back(50.);
    posAnodes.push_back(150.);
    posAnodes.push_back(250.);
    posAnodes.push_back(350.);
    // NDLAr edges
    xBoundaries.push_back(-350.);
    xBoundaries.push_back(350.);
    zBoundaries.push_back(418.);
    zBoundaries.push_back(913.);
    yBoundaries.push_back(-216.);
    yBoundaries.push_back(82.);
  }
  else if ( parameters.fDetector == 1 ) {
    // 2x2 anodes
    posAnodes.push_back(3.0652);
    posAnodes.push_back(63.9273);
    posAnodes.push_back(-3.0652);
    posAnodes.push_back(-63.9273);
    // 2x2 edges
    xBoundaries.push_back(-63.9273);
    xBoundaries.push_back(63.9273);
    zBoundaries.push_back(-63.9273);
    zBoundaries.push_back(63.9273);
    yBoundaries.push_back(-63.9273);
    yBoundaries.push_back(63.9273);
  }

  /////////////////////////////////////////////////////////
  /// Set up the necessary pieces for muon momentum vs range spline: CSDA
  /// copied from LArSoft -> LArReco -> TrackMomentumCalculator
  ///   v09_26_02
  std::array<float, 29> Range_grampercm2{
    {9.833E-1, 1.786E0, 3.321E0, 6.598E0, 1.058E1, 3.084E1, 4.250E1, 6.732E1, 1.063E2, 1.725E2,
     2.385E2,  4.934E2, 6.163E2, 8.552E2, 1.202E3, 1.758E3, 2.297E3, 4.359E3, 5.354E3, 7.298E3,
     1.013E4,  1.469E4, 1.910E4, 3.558E4, 4.326E4, 5.768E4, 7.734E4, 1.060E5, 1.307E5}};
  for (float& value : Range_grampercm2) {
    value /= 1.396; // convert to cm
  }
  constexpr std::array<float, 29> KE_MeV{
    {10,    14,    20,    30,    40,     80,     100,    140,    200,   300,
     400,   800,   1000,  1400,  2000,   3000,   4000,   8000,   10000, 14000,
     20000, 30000, 40000, 80000, 100000, 140000, 200000, 300000, 400000}};
  TGraph const KEvsR{29, Range_grampercm2.data(), KE_MeV.data()};
  TSpline3 const KEvsR_spline3_muon{"KEvsRS", &KEvsR};

  /////////////////////////////////////////////////////////

  TFile *fileSource = TFile::Open(parameters.fileName.c_str(), "READ");
  if (!fileSource)
  {
      std::cout << "Error in ProcessPostReco(): can't open file " << parameters.fileName << std::endl;
      return;
  }

  TTree *recoTree = dynamic_cast<TTree *>(fileSource->Get("LArRecoND"));
  if (!recoTree)
  {
      std::cout << "Could not find the event tree LArRecoND" << std::endl;
      fileSource->Close();
      return;
  }

  std::unique_ptr<LArRecoNDFormat> pandoraIn = std::make_unique<LArRecoNDFormat>(recoTree);

  long nEntries = recoTree->GetEntries();

  std::cout << "Runninng ProcessEvents on " << nEntries << " entries with pixel pitch " 
	    << parameters.pixelPitch << " and track/shower separation score of " << parameters.trackScoreCut << std::endl;

  // Create the class where we'll store the output info
  //
  NDRecoOutputData fOut( parameters.outfileName, parameters.runTrackFit, parameters.runShowerFit, parameters.fShouldRunPID);

  // Loop events
  for ( long entryIdx = 0; entryIdx < nEntries; ++entryIdx ) {
    pandoraIn->GetEntry(entryIdx);

    // Fill up the branches of basic output
    fOut.FillBasicBranches(pandoraIn);

    // Track fit vectors of importance
    std::vector<float> trkStartX, trkStartY, trkStartZ, trkEndX, trkEndY, trkEndZ;
    std::vector<float> trkStartDirX, trkStartDirY, trkStartDirZ, trkEndDirX, trkEndDirY, trkEndDirZ;
    std::vector<float> trkLen, trk_KEFromLength_muon, trk_KEFromLength_proton, trk_pFromLength_muon, trk_pFromLength_proton;
    std::vector<bool> trkContained;

    // Track fit calo vectors of importance
    std::vector<float> trackFitTrackCaloE, trackFitVisE;

    // Track fit POINT values
    std::vector<int> trackFitSliceId, trackFitPfoId;
    std::vector<float> trackFitX, trackFitY, trackFitZ;
    std::vector<float> trackFitQ, trackFitRR, trackFitdx, trackFitdQdx, trackFitdEdx;

    // Per Particle PID
    std::vector<float> pid_muScore, pid_piScore, pid_kScore, pid_proScore;
    std::vector<int> pid_pdg, pid_ndf;

    //Shower fit vectors of importance
    std::vector<float> shwrCentroidX, shwrCentroidY, shwrCentroidZ, shwrStartX, shwrStartY, shwrStartZ;
    std::vector<float> shwrDirX, shwrDirY, shwrDirZ;
    std::vector<float> shwrLen;

    std::vector<int> shwrSliceId, shwrClusterId;

    std::vector<float> startTrkSlidingFitDirX,  startTrkSlidingFitDirY,  startTrkSlidingFitDirZ;    
    std::vector<float> shwrStartPointsX, shwrStartPointsY, shwrStartPointsZ;
    std::vector<double> shwrdEdx;
    std::vector<float> shwrTotalE;
    std::vector<int> shwrStartPointsRecoId;
    std::vector<float> minProjection, medianQ, chargePerHit, chargePerHitStartPoints, pitchValue;
    std::vector<float> shwrPCAX, shwrPCAY, shwrPCAZ;
    std::vector<float> shwrStartHitPositionX, shwrStartHitPositionY, shwrStartHitPositionZ;

    // Loop particles in the event
    unsigned int nParticles = pandoraIn->m_clusterID->size();
    for ( unsigned int particleIdx=0; particleIdx < nParticles; ++particleIdx ) {
      float trackScore = pandoraIn->m_trackScore->at(particleIdx);

      if ( !parameters.runTrackFit && !parameters.runShowerFit ) continue;

      int hitCounter(0);

      // Read in the vertex and point vector that will be the input to the track and shower fits
      int sliceID = pandoraIn->m_sliceID->at(particleIdx);
      int clusterID = pandoraIn->m_clusterID->at(particleIdx);
      CartesianVector vertexVector(pandoraIn->m_nuVtxX->at(particleIdx), pandoraIn->m_nuVtxY->at(particleIdx), pandoraIn->m_nuVtxZ->at(particleIdx));
      CaloHitList caloHitList;
      for ( unsigned int idxHits = 0; idxHits < pandoraIn->m_recoHitId->size(); ++idxHits ) {
	if ( pandoraIn->m_recoHitSliceId->at(idxHits)==sliceID && pandoraIn->m_recoHitClusterId->at(idxHits)==clusterID ) {
	  // Skip hit if it fails the threshold
	  if ( parameters.applyThreshold && pandoraIn->m_recoHitE->at(idxHits) < parameters.thresholdVal )
	    continue;
	  CartesianVector thisHit(pandoraIn->m_recoHitX->at(idxHits), pandoraIn->m_recoHitY->at(idxHits), pandoraIn->m_recoHitZ->at(idxHits));
	  lar_content::LArCaloHitParameters chParams;
	  chParams.m_positionVector = thisHit;
	  chParams.m_expectedDirection = pandora::CartesianVector(0.f, 0.f, 1.f);
	  chParams.m_cellNormalVector = pandora::CartesianVector(0.f, 0.f, 1.f);
	  chParams.m_cellGeometry = pandora::RECTANGULAR;
	  chParams.m_cellSize0 = parameters.pixelPitch;
	  chParams.m_cellSize1 = parameters.pixelPitch;
	  chParams.m_cellThickness = parameters.pixelPitch;
	  chParams.m_nCellRadiationLengths = 1.f;
	  chParams.m_nCellInteractionLengths = 1.f;
	  chParams.m_time = 0.f;
	  chParams.m_inputEnergy = pandoraIn->m_recoHitE->at(idxHits);
	  chParams.m_mipEquivalentEnergy = pandoraIn->m_recoHitE->at(idxHits);
	  chParams.m_electromagneticEnergy = pandoraIn->m_recoHitE->at(idxHits);
	  chParams.m_hadronicEnergy = pandoraIn->m_recoHitE->at(idxHits);
	  chParams.m_isDigital = false;
	  chParams.m_hitType = pandora::TPC_3D;
	  chParams.m_hitRegion = pandora::SINGLE_REGION;
	  chParams.m_layer = 0;
	  chParams.m_isInOuterSamplingLayer = false;
	  chParams.m_pParentAddress = (void*)(static_cast<uintptr_t>(++hitCounter));
	  chParams.m_larTPCVolumeId = 0;
	  chParams.m_daughterVolumeId = 0;
	  // push back the calo hit
	  lar_content::LArCaloHit *ch = new lar_content::LArCaloHit(chParams);
	  caloHitList.push_back( ch );
	}
      } // loop hits

      // Fit it as a track?
      if ( parameters.runTrackFit && (parameters.trackScoreCut < 0. || trackScore >= parameters.trackScoreCut) ) {
	// Run the track fit info:
	// TODO: Make the MinTrajectoryPoints(default=2) and SlidingFitHalfWindow(20) configurable
	int minTrajectoryPoints = 2;
	float slidingFitHalfWindow = 20;

	std::vector<float> trackVecDEDX;
	std::vector<float> trackVecRR;
	std::vector<float> trackVecDX;

	bool filledPID=false; // We'll check at the end and fill PID with bogus info if PID info not filled

	lar_content::LArTrackStateVector trackStateVector;
	std::vector<int> indexVector;
	bool trackStateSuccess=false;
	try {
	  lar_content::LArPfoHelper::GetSlidingFitTrajectory( &caloHitList, vertexVector, slidingFitHalfWindow, parameters.pixelPitch, trackStateVector, &indexVector, true );
	  trackStateSuccess=true;
	}
	catch (const pandora::StatusCodeException&) {
	  trackStateSuccess=false;
	}

	// If user has set the Voxelize Z function, then rerun the track fit, starting from the output of the first fit
	lar_content::LArTrackStateVector trackStateVector_v2;
	std::vector<int> indexVector_v2;
        bool trackStateSuccess_v2=false;
	if ( parameters.voxelizeZ && trackStateSuccess ) {
	  if( parameters.verbosity >= 1 ) {
	    std::cout << "    INFO: Since voxelization is turned on, we will take the output of the track fit and try to voxelize now." << std::endl;
	    std::cout << "    ----> Input track has " << trackStateVector.size() << " track points." << std::endl;
	  }
	  try
	  {
	    int hitCounter_v1p5(0);
	    int hitCounter_v2(0);

	    // Initial calohit vector
	    std::vector< lar_content::LArCaloHit* > caloHitVect_v1;
	    for (unsigned int idxPt=0; idxPt < trackStateVector.size(); ++idxPt ) {
	      const lar_content::LArTrackState& trackState = trackStateVector.at(idxPt);
	      lar_content::LArCaloHitParameters chParams;
	      chParams.m_positionVector = trackState.GetCaloHit()->GetPositionVector();
	      chParams.m_expectedDirection = pandora::CartesianVector(0.f, 0.f, 1.f);
	      chParams.m_cellNormalVector = pandora::CartesianVector(0.f, 0.f, 1.f);
	      chParams.m_cellGeometry = pandora::RECTANGULAR;
	      chParams.m_cellSize0 = parameters.pixelPitch;
	      chParams.m_cellSize1 = parameters.pixelPitch;
	      chParams.m_cellThickness = parameters.pixelPitch;
	      chParams.m_nCellRadiationLengths = 1.f;
	      chParams.m_nCellInteractionLengths = 1.f;
	      chParams.m_time = 0.f;
	      chParams.m_inputEnergy = trackState.GetCaloHit()->GetInputEnergy();
	      chParams.m_mipEquivalentEnergy = trackState.GetCaloHit()->GetInputEnergy();
	      chParams.m_electromagneticEnergy = trackState.GetCaloHit()->GetInputEnergy();
	      chParams.m_hadronicEnergy = trackState.GetCaloHit()->GetInputEnergy();
	      chParams.m_isDigital = false;
	      chParams.m_hitType = pandora::TPC_3D;
	      chParams.m_hitRegion = pandora::SINGLE_REGION;
	      chParams.m_layer = 0;
	      chParams.m_isInOuterSamplingLayer = false;
	      chParams.m_pParentAddress = (void*)(static_cast<uintptr_t>(++hitCounter_v1p5));
	      chParams.m_larTPCVolumeId = 0;
	      chParams.m_daughterVolumeId = 0;
	      lar_content::LArCaloHit *ch = new lar_content::LArCaloHit(chParams);
	      caloHitVect_v1.push_back( ch );
	    }
	    // Now let's construct the version that goes into the second pass track fit.
	    // 1. Loop through the vector and for each element, gather all the consecutive elements within epsilon of the z value
	    // 2. Within this subset, find the maximum Q hit, start here
	    //     a. Gather this hit and the ones within an x, y distance of the voxel setting
	    //     b. Make a new calohit that is the weighted mean of the (x, y, z) of these hits and the sum of the Q values
	    // 3. Repeat on the maximal Q value of the hits letf and continue repeating till all hits are swept up
	    // 4. Run track fit on this.
	    CaloHitList caloHitList_v2;
	    for ( unsigned int idxHit=0; idxHit < caloHitVect_v1.size(); ++idxHit ) {
	      if ( parameters.verbosity >= 2 ) std::cout << "      hit idx " << idxHit << " of " << caloHitVect_v1.size() << std::endl;
	      // Step 1
	      float thisZ = caloHitVect_v1.at(idxHit)->GetPositionVector().GetZ();
	      std::vector< lar_content::LArCaloHit* > caloHitVect_tmp;
	      caloHitVect_tmp.push_back( caloHitVect_v1.at(idxHit) );
	      bool stopLoop=false;
	      while ( !stopLoop && idxHit < caloHitVect_v1.size()-1 ) {
		if ( fabs(caloHitVect_v1.at(idxHit+1)->GetPositionVector().GetZ()-thisZ) < std::numeric_limits<float>::epsilon() ) {
		  caloHitVect_tmp.push_back( caloHitVect_v1.at(idxHit+1) );
		  idxHit+=1;
		}
		else stopLoop=true;
	      } // found all the hits that we need to check
	      if( parameters.verbosity >= 2 ) std::cout << "      --> At this stage of the voxelization, we have " << caloHitVect_tmp.size() << " hits to possibly merge." << std::endl;
	      // Step 2-3
	      if ( caloHitVect_tmp.size() == 1 ) {
		lar_content::LArCaloHitParameters chParams;
		chParams.m_positionVector = caloHitVect_tmp.at(0)->GetPositionVector();
		chParams.m_expectedDirection = pandora::CartesianVector(0.f, 0.f, 1.f);
		chParams.m_cellNormalVector = pandora::CartesianVector(0.f, 0.f, 1.f);
		chParams.m_cellGeometry = pandora::RECTANGULAR;
		chParams.m_cellSize0 = parameters.pixelPitch;
		chParams.m_cellSize1 = parameters.pixelPitch;
		chParams.m_cellThickness = parameters.pixelPitch;
		chParams.m_nCellRadiationLengths = 1.f;
		chParams.m_nCellInteractionLengths = 1.f;
		chParams.m_time = 0.f;
		chParams.m_inputEnergy = caloHitVect_tmp.at(0)->GetInputEnergy();
		chParams.m_mipEquivalentEnergy = caloHitVect_tmp.at(0)->GetInputEnergy();
		chParams.m_electromagneticEnergy = caloHitVect_tmp.at(0)->GetInputEnergy();
		chParams.m_hadronicEnergy = caloHitVect_tmp.at(0)->GetInputEnergy();
		chParams.m_isDigital = false;
		chParams.m_hitType = pandora::TPC_3D;
		chParams.m_hitRegion = pandora::SINGLE_REGION;
		chParams.m_layer = 0;
		chParams.m_isInOuterSamplingLayer = false;
		chParams.m_pParentAddress = (void*)(static_cast<uintptr_t>(++hitCounter_v2));
		chParams.m_larTPCVolumeId = 0;
		chParams.m_daughterVolumeId = 0;
		lar_content::LArCaloHit *ch = new lar_content::LArCaloHit(chParams);
		caloHitList_v2.push_back( ch );
	      }
	      else {
		while ( caloHitVect_tmp.size() > 0 ) {
		  float maxQ = 0.;
		  float maxQ_X = 0.;
		  float maxQ_Y = 0.;
		  for ( unsigned int idxHit_inner=0; idxHit_inner < caloHitVect_tmp.size(); ++idxHit_inner ) {
		    if ( caloHitVect_tmp.at(idxHit_inner)->GetInputEnergy() > maxQ ) {
		      maxQ = caloHitVect_tmp.at(idxHit_inner)->GetInputEnergy();
		      maxQ_X = caloHitVect_tmp.at(idxHit_inner)->GetPositionVector().GetX();
		      maxQ_Y = caloHitVect_tmp.at(idxHit_inner)->GetPositionVector().GetY();
		    }
		  }
		  if ( parameters.verbosity >= 2 ) std::cout << "      --> Max Hit X = " << maxQ_X << ", Y = " << maxQ_Y << ", Q = " << maxQ << std::endl;
		  std::vector<float> xs, ys, zs, qs;
		  std::vector<unsigned int> toDelete;
		  for ( unsigned int idxHit_inner=0; idxHit_inner < caloHitVect_tmp.size(); ++idxHit_inner ) {
		    float thisX_inner = caloHitVect_tmp.at(idxHit_inner)->GetPositionVector().GetX();
		    float thisY_inner = caloHitVect_tmp.at(idxHit_inner)->GetPositionVector().GetY();
		    if ( parameters.verbosity >= 2 ) std::cout << "      --> This Hit X = " << thisX_inner << ", Y = " << thisY_inner << ", Q = " << caloHitVect_tmp.at(idxHit_inner)->GetInputEnergy() <<std::endl;
                    if ( std::sqrt( std::pow(thisX_inner-maxQ_X,2) + std::pow(thisY_inner-maxQ_Y,2) ) < parameters.voxelZHW ) {
		      float thisZ_inner = caloHitVect_tmp.at(idxHit_inner)->GetPositionVector().GetZ();
		      float thisQ_inner = caloHitVect_tmp.at(idxHit_inner)->GetInputEnergy();
		      xs.push_back(thisX_inner);
		      ys.push_back(thisY_inner);
		      zs.push_back(thisZ_inner);
		      qs.push_back(thisQ_inner);
		      toDelete.push_back(idxHit_inner);
                    }
                  }
		  if ( parameters.verbosity >= 2 ) std::cout << "      --> Making a new hit from " << xs.size() << " hit(s) and deleting " << toDelete.size() << " hits." << std::endl;
		  std::vector< lar_content::LArCaloHit* > caloHitVect_tmp_prev = caloHitVect_tmp;
		  caloHitVect_tmp.clear();
		  for ( unsigned int idxCopy=0; idxCopy < caloHitVect_tmp_prev.size(); ++idxCopy ) {
		    bool skipCopy = false;
		    for ( unsigned int checkIdx=0; checkIdx < toDelete.size(); ++checkIdx ) {
		      if ( idxCopy == toDelete[checkIdx] ) {
			skipCopy = true;
			break;
		      }
		    }
		    if ( skipCopy ) continue;
		    caloHitVect_tmp.push_back( caloHitVect_tmp_prev.at(idxCopy) );
		  }
		  // Make new hit:
		  float newHitX(0.), newHitY(0.), newHitZ(0.), newHitQ(0.);
		  for ( unsigned int idxUse=0; idxUse<xs.size(); ++idxUse ) {
		    newHitX+=xs[idxUse]*qs[idxUse];
		    newHitY+=ys[idxUse]*qs[idxUse];
		    newHitZ+=zs[idxUse]*qs[idxUse];
		    newHitQ+=qs[idxUse];
		  }
		  if ( newHitQ > 0. ){
		    newHitX/=newHitQ;
		    newHitY/=newHitQ;
		    newHitZ/=newHitQ;
		  }
		  lar_content::LArCaloHitParameters chParams;
		  chParams.m_positionVector = {newHitX, newHitY, newHitZ};
		  chParams.m_expectedDirection = pandora::CartesianVector(0.f, 0.f, 1.f);
		  chParams.m_cellNormalVector = pandora::CartesianVector(0.f, 0.f, 1.f);
		  chParams.m_cellGeometry = pandora::RECTANGULAR;
		  chParams.m_cellSize0 = parameters.pixelPitch;
		  chParams.m_cellSize1 = parameters.pixelPitch;
		  chParams.m_cellThickness = parameters.pixelPitch;
		  chParams.m_nCellRadiationLengths = 1.f;
		  chParams.m_nCellInteractionLengths = 1.f;
		  chParams.m_time = 0.f;
		  chParams.m_inputEnergy = newHitQ;
		  chParams.m_mipEquivalentEnergy = newHitQ;
		  chParams.m_electromagneticEnergy = newHitQ;
		  chParams.m_hadronicEnergy = newHitQ;
		  chParams.m_isDigital = false;
		  chParams.m_hitType = pandora::TPC_3D;
		  chParams.m_hitRegion = pandora::SINGLE_REGION;
		  chParams.m_layer = 0;
		  chParams.m_isInOuterSamplingLayer = false;
		  chParams.m_pParentAddress = (void*)(static_cast<uintptr_t>(++hitCounter_v2));
		  chParams.m_larTPCVolumeId = 0;
		  chParams.m_daughterVolumeId = 0;
		  lar_content::LArCaloHit *ch = new lar_content::LArCaloHit(chParams);
		  caloHitList_v2.push_back( ch );
		  if( parameters.verbosity >= 2 )
		    std::cout << "      --> After this particular voxelization, we have " << caloHitVect_tmp.size() <<" hits remaining to possibly merge.\n"
			      << "          and caloHitList_v2 size is " << caloHitList_v2.size() << std::endl;
		}
	      } // Step 2-3
	    } // Steps 1-3
	    // Step 4
	    if( parameters.verbosity >= 1 ) std::cout << "    ----> DONE with the merging. Now running the new track fit." << std::endl;
	    lar_content::LArPfoHelper::GetSlidingFitTrajectory( &caloHitList_v2, vertexVector, slidingFitHalfWindow, parameters.pixelPitch, trackStateVector_v2, &indexVector_v2, true );

	    trackStateSuccess_v2=true;
	  }
	  catch(const pandora::StatusCodeException&) {
	    trackStateSuccess_v2=false;
	  }
	}

	lar_content::LArTrackStateVector trackStateVector_out = (parameters.voxelizeZ && trackStateSuccess_v2 ) ? trackStateVector_v2 : trackStateVector;
	if( parameters.verbosity >= 1 ) std::cout << "    INFO: The track state vector we are using for calorimetry analysis has " << trackStateVector_out.size() << " points." << std::endl;

	// Extract the track fit info
	if (!trackStateSuccess || trackStateVector.size() < minTrajectoryPoints) {
	  trkStartX.push_back(-9999.);
	  trkStartY.push_back(-9999.);
	  trkStartZ.push_back(-9999.);
	  trkStartDirX.push_back(1.);
	  trkStartDirY.push_back(0.);
	  trkStartDirZ.push_back(0.);
	  trkEndX.push_back(-9999.);
          trkEndY.push_back(-9999.);
          trkEndZ.push_back(-9999.);
          trkEndDirX.push_back(1.);
          trkEndDirY.push_back(0.);
          trkEndDirZ.push_back(0.);
	  trkLen.push_back(0.);
	  trkContained.push_back(false);
	  trackFitTrackCaloE.push_back(0.);
	  trackFitVisE.push_back(0.);
	  trk_KEFromLength_muon.push_back(0.);
	  trk_KEFromLength_proton.push_back(0.);
	  trk_pFromLength_muon.push_back(0.);
	  trk_pFromLength_proton.push_back(0.);
	}
	else {
	  const lar_content::LArTrackState& trackStateStart =
	    ( parameters.useVoxelizedStartStop && (trackStateSuccess_v2 && trackStateVector_out.size() >= minTrajectoryPoints) ) ? 
	    trackStateVector_out.front() : 
	    trackStateVector.front();
	  trkStartX.push_back(trackStateStart.GetPosition().GetX());
	  trkStartY.push_back(trackStateStart.GetPosition().GetY());
	  trkStartZ.push_back(trackStateStart.GetPosition().GetZ());
	  trkStartDirX.push_back(trackStateStart.GetDirection().GetX());
          trkStartDirY.push_back(trackStateStart.GetDirection().GetY());
          trkStartDirZ.push_back(trackStateStart.GetDirection().GetZ());
	  const lar_content::LArTrackState& trackStateEnd =
	    ( parameters.useVoxelizedStartStop && (trackStateSuccess_v2 && trackStateVector_out.size() >= minTrajectoryPoints) ) ?
            trackStateVector_out.back() :
            trackStateVector.back();
	  trkEndX.push_back(trackStateEnd.GetPosition().GetX());
          trkEndY.push_back(trackStateEnd.GetPosition().GetY());
          trkEndZ.push_back(trackStateEnd.GetPosition().GetZ());
          trkEndDirX.push_back(trackStateEnd.GetDirection().GetX());
          trkEndDirY.push_back(trackStateEnd.GetDirection().GetY());
          trkEndDirZ.push_back(trackStateEnd.GetDirection().GetZ());

	  // is the track contained?
	  bool thisTrackContained = true;
	  if      ( trackStateStart.GetPosition().GetX() < xBoundaries[0]+5. ) thisTrackContained = false;
	  else if ( trackStateStart.GetPosition().GetX() > xBoundaries[1]-5. ) thisTrackContained = false;
	  else if ( trackStateStart.GetPosition().GetY() < yBoundaries[0]+5. ) thisTrackContained = false;
          else if ( trackStateStart.GetPosition().GetY() > yBoundaries[1]-5. ) thisTrackContained = false;
	  else if ( trackStateStart.GetPosition().GetZ() < zBoundaries[0]+5. ) thisTrackContained = false;
          else if ( trackStateStart.GetPosition().GetZ() > zBoundaries[1]-5. ) thisTrackContained = false;
	  else if ( trackStateEnd.GetPosition().GetX() < xBoundaries[0]+5. ) thisTrackContained = false;
          else if ( trackStateEnd.GetPosition().GetX() > xBoundaries[1]-5. ) thisTrackContained = false;
          else if ( trackStateEnd.GetPosition().GetY() < yBoundaries[0]+5. ) thisTrackContained = false;
          else if ( trackStateEnd.GetPosition().GetY() > yBoundaries[1]-5. ) thisTrackContained = false;
          else if ( trackStateEnd.GetPosition().GetZ() < zBoundaries[0]+5. ) thisTrackContained = false;
          else if ( trackStateEnd.GetPosition().GetZ() > zBoundaries[1]-5. ) thisTrackContained = false;
	  trkContained.push_back(thisTrackContained);

	  float trklength = 0.;
	  // Get the length going point to point
	  if ( parameters.useVoxelizedStartStop && (trackStateSuccess_v2 && trackStateVector_out.size() >= minTrajectoryPoints) ) {
	    for (unsigned int idxPt=0; idxPt < trackStateVector_out.size()-1; ++idxPt) {
              const lar_content::LArTrackState& trackState = trackStateVector_out.at(idxPt);
              const lar_content::LArTrackState& trackStateNext = trackStateVector_out.at(idxPt+1);
              trklength+=std::sqrt( trackState.GetPosition().GetDistanceSquared( trackStateNext.GetPosition() ) );
            }
	  }
	  else {
	    for (unsigned int idxPt=0; idxPt < trackStateVector.size()-1; ++idxPt) {
	      const lar_content::LArTrackState& trackState = trackStateVector.at(idxPt);
	      const lar_content::LArTrackState& trackStateNext = trackStateVector.at(idxPt+1);
	      trklength+=std::sqrt( trackState.GetPosition().GetDistanceSquared( trackStateNext.GetPosition() ) );
	    }
	  }
	  trkLen.push_back(trklength);

	  // Track momentum from range:
	  trk_KEFromLength_proton.push_back( KEFromRange_proton(trklength) );
	  trk_pFromLength_proton.push_back( pFromRange_proton(trklength) );
	  float KEFromLength_muon = KEvsR_spline3_muon.Eval(trklength);
	  if (KEFromLength_muon > 0.) {
	    trk_KEFromLength_muon.push_back( KEFromLength_muon/1000. );
	    trk_pFromLength_muon.push_back( std::sqrt((KEFromLength_muon * KEFromLength_muon) + (2 * 105.7 * KEFromLength_muon))/1000. );
	  }
	  else {
	    trk_KEFromLength_muon.push_back( 0. );
	    trk_pFromLength_muon.push_back( 0. );
	  }

	  // Track calorimetry --> very rough first pass basically reimplemented from other test branch:
	  // ! Consider the first and last points, but here we only have one side of dx
	  // ! Does not do spacecharge, diffusion, etc. corrections at least yet
	  float summedTrkE = 0.;
	  float summedQinTrk = 0.;

	  // If we aren't going to do the track calorimetry, then say the track caloE = 0
          if ( !(trackStateVector_out.size() >= minTrajectoryPoints) ) {
            trackFitTrackCaloE.push_back(0.);
	    trackFitVisE.push_back(0.);
	  }

	  if ( trackStateVector_out.size() >= minTrajectoryPoints ) {
	    float lengthSoFar = 0.;
	    for (unsigned int idxPt=0; idxPt < trackStateVector_out.size(); ++idxPt ) {
	      const lar_content::LArTrackState& trackState = trackStateVector_out.at(idxPt);

	      // charge
	      float hitQ = trackState.GetCaloHit()->GetInputEnergy();
	      // residual range
	      if ( idxPt > 0 ) {
		const lar_content::LArTrackState& trackStatePrev = trackStateVector_out.at(idxPt-1);
		lengthSoFar+=std::sqrt( trackStatePrev.GetPosition().GetDistanceSquared( trackState.GetPosition() ) );
	      }
	      float hitRR = trklength - lengthSoFar;
	      // dx
	      float hitdx = 0.;
	      if ( idxPt == 0 ){
		if ( idxPt < trackStateVector_out.size()-1 ) {
                  const lar_content::LArTrackState& trackStateNext = trackStateVector_out.at(idxPt+1);
                  hitdx = std::sqrt( trackState.GetPosition().GetDistanceSquared( trackStateNext.GetPosition() ) )/2.;
                }
	      }
	      else {
		const lar_content::LArTrackState& trackStatePrev = trackStateVector_out.at(idxPt-1);
		// Middle Points
		if ( idxPt < trackStateVector_out.size()-1 ) {
                  const lar_content::LArTrackState& trackStateNext = trackStateVector_out.at(idxPt+1);
                  hitdx = std::sqrt( trackStatePrev.GetPosition().GetDistanceSquared( trackStateNext.GetPosition() ) )/2.;
                }
		// Last Point
		else {
                  hitdx = std::sqrt( trackStatePrev.GetPosition().GetDistanceSquared( trackState.GetPosition() ) )/2.;
                }
	      }
	      // dQdx
	      float hitdQdx = hitdx > 0. ? hitQ/hitdx : -5.f;

	      // Lifetime correction
	      if ( parameters.fShouldCorrectLifetime ) {
		hitdQdx*=( 1000. * LifetimeCorrectionFactor(posAnodes, trackState.GetPosition().GetX(), parameters.fElectronLifetime, parameters.fElectronDriftSpeed) ); // turn ke- to e- and do lifetime correction
		summedQinTrk += ( 1000. * hitQ * LifetimeCorrectionFactor(posAnodes, trackState.GetPosition().GetX(), parameters.fElectronLifetime, parameters.fElectronDriftSpeed) );
	      }
	      // Recombination correction
	      float hitdEdx = dEdxWithRecombination(parameters, hitdQdx);

	      // Outputs
	      trackFitSliceId.push_back(sliceID);
	      trackFitPfoId.push_back(clusterID);
	      trackFitX.push_back(trackState.GetPosition().GetX());
	      trackFitY.push_back(trackState.GetPosition().GetY());
	      trackFitZ.push_back(trackState.GetPosition().GetZ());
	      trackFitQ.push_back(hitQ);
	      trackFitRR.push_back(hitRR);
	      trackFitdx.push_back(hitdx);
	      trackFitdQdx.push_back(hitdQdx);
	      trackFitdEdx.push_back(hitdEdx);

	      trackVecDX.push_back(hitdx);
	      trackVecDEDX.push_back(hitdEdx);
	      trackVecRR.push_back(hitRR);

	      summedTrkE+=hitdEdx*hitdx; // sum up the energy along the track

	    } // loop points
	    // And now that we have dE/dx for all points, we can use the sum of that all to get the track calo E
	    trackFitTrackCaloE.push_back(summedTrkE);
	    // And calculate the total VisE for the track:
	    trackFitVisE.push_back( eVisWithRecombination(parameters, summedQinTrk) );

	    // Particle ID here
	    if ( parameters.fShouldRunPID ) {
	      if ( parameters.fPIDAlgChi2PID ) {
		// as in https://github.com/LArSoft/larana/blob/develop/larana/ParticleIdentification/Chi2PIDAlg.cxx#L90
		float chi2pro = 0.;
		float chi2ka = 0.;
		float chi2pi = 0.;
		float chi2mu = 0.;
		int nbins_dedx_range = parameters.templatesdEdxRR.at("proton")->GetNbinsX();
		int npts = 0;
		for ( unsigned int idxCaloPt=0; idxCaloPt < trackVecDEDX.size(); ++idxCaloPt ) {
		  if ( idxCaloPt == 0 || idxCaloPt == trackVecDEDX.size()-1 ) continue; // ignore 1st and last point
		  if ( trackVecDEDX[idxCaloPt] > 1000. ) continue; // ignore if too high dEdx
		  if ( trackVecDEDX[idxCaloPt] < parameters.fChi2RestrictDEDXLo ) continue; // also, optionally restrict unexpected low dE/dx. By default just requires it to be positive.
		  if ( parameters.fChi2RestrictDX && 
		       ( trackVecDX[idxCaloPt] < parameters.fChi2RestrictDXLo || 
			 (parameters.fChi2RestrictDXHi > 0. && trackVecDX[idxCaloPt] > parameters.fChi2RestrictDXHi) ) ) {
		    continue; // optionally skip this point if dx too small/large
		  }
		  int bin = parameters.templatesdEdxRR.at("proton")->FindBin(trackVecRR[idxCaloPt]);
		  if ( bin >= 1 && bin <= nbins_dedx_range ) {
		    // Content
		    float bincpro = parameters.templatesdEdxRR.at("proton")->GetBinContent(bin);
		    if ( bincpro < 1e-6 ) bincpro = ( parameters.templatesdEdxRR.at("proton")->GetBinContent(bin-1) + parameters.templatesdEdxRR.at("proton")->GetBinContent(bin+1) )/2.;
		    float bincka = parameters.templatesdEdxRR.at("kaon")->GetBinContent(bin);
                    if ( bincka < 1e-6 ) bincka = ( parameters.templatesdEdxRR.at("kaon")->GetBinContent(bin-1) + parameters.templatesdEdxRR.at("kaon")->GetBinContent(bin+1) )/2.;
		    float bincpi = parameters.templatesdEdxRR.at("pion")->GetBinContent(bin);
                    if ( bincpi < 1e-6 ) bincpi = ( parameters.templatesdEdxRR.at("pion")->GetBinContent(bin-1) + parameters.templatesdEdxRR.at("pion")->GetBinContent(bin+1) )/2.;
		    float bincmu = parameters.templatesdEdxRR.at("muon")->GetBinContent(bin);
                    if ( bincmu < 1e-6 ) bincmu = ( parameters.templatesdEdxRR.at("muon")->GetBinContent(bin-1) + parameters.templatesdEdxRR.at("muon")->GetBinContent(bin+1) )/2.;
		    // Error
		    float binepro = parameters.templatesdEdxRR.at("proton")->GetBinError(bin);
                    if ( binepro < 1e-6 ) binepro = ( parameters.templatesdEdxRR.at("proton")->GetBinError(bin-1) + parameters.templatesdEdxRR.at("proton")->GetBinError(bin+1) )/2.;
                    float bineka = parameters.templatesdEdxRR.at("kaon")->GetBinError(bin);
                    if ( bineka < 1e-6 ) bineka = ( parameters.templatesdEdxRR.at("kaon")->GetBinError(bin-1) + parameters.templatesdEdxRR.at("kaon")->GetBinError(bin+1) )/2.;
                    float binepi = parameters.templatesdEdxRR.at("pion")->GetBinError(bin);
                    if ( binepi < 1e-6 ) binepi = ( parameters.templatesdEdxRR.at("pion")->GetBinError(bin-1) + parameters.templatesdEdxRR.at("pion")->GetBinError(bin+1) )/2.;
                    float binemu = parameters.templatesdEdxRR.at("muon")->GetBinError(bin);
                    if ( binemu < 1e-6 ) binemu = ( parameters.templatesdEdxRR.at("muon")->GetBinError(bin-1) + parameters.templatesdEdxRR.at("muon")->GetBinError(bin+1) )/2.;
		    float errdedx = 0.04231 + 0.0001783 * trackVecDEDX[idxCaloPt] * trackVecDEDX[idxCaloPt];
		    errdedx *= trackVecDEDX[idxCaloPt];
		    float errdedx_square = errdedx*errdedx;
		    // chi2 values
		    float thisPointDEDX = trackVecDEDX[idxCaloPt];
		    if ( !parameters.fApplyCalibrationFudgeFactor && parameters.fApplyCalibrationFudgeFactor_PID ) thisPointDEDX*=parameters.fCalibrationFudgeFactor;
		    chi2pro += std::pow(thisPointDEDX - bincpro, 2) / (binepro*binepro + errdedx_square);
		    chi2ka += std::pow(thisPointDEDX - bincka, 2) / (bineka*bineka + errdedx_square);
		    chi2pi += std::pow(thisPointDEDX - bincpi, 2) / (binepi*binepi + errdedx_square);
		    chi2mu += std::pow(thisPointDEDX - bincmu, 2) / (binemu*binemu + errdedx_square);
		    npts+=1;
		  } // within bins
		} // loop calo points

		if ( npts > 0 ) {
		  int thisPDG = 0;
		  float thisChi2 = std::numeric_limits<float>::max();
		  if ( chi2pro/npts < thisChi2 ) {
		    thisPDG = 2212;
		    thisChi2 = chi2pro/npts;
		  }
		  if ( chi2ka/npts < thisChi2 ) {
                    thisPDG = 321;
                    thisChi2 = chi2ka/npts;
                  }
		  if ( chi2pi/npts < thisChi2 ) {
                    thisPDG = 211;
                    thisChi2 = chi2pi/npts;
                  }
		  if ( chi2mu/npts < thisChi2 ) {
                    thisPDG = 13;
                    thisChi2 = chi2mu/npts;
                  }
		  // prediction is minimum chi2/npts
		  filledPID=true;
		  pid_pdg.push_back(thisPDG);
		  pid_ndf.push_back(npts);
		  pid_muScore.push_back(chi2mu/npts);
		  pid_piScore.push_back(chi2pi/npts);
		  pid_kScore.push_back(chi2ka/npts);
		  pid_proScore.push_back(chi2pro/npts);
		}
	      } // use Chi2PID
	    } // getting PID
	    ///////////////////////

	  } // if trackstate has stuff needed to do dEdx
	} // if we have a track state

	if ( !filledPID ) {
	  // if PID isn't filled then we need to add in the defaults for this track
	  pid_pdg.push_back(0);
	  pid_ndf.push_back(0);
	  pid_muScore.push_back(-5.);
	  pid_piScore.push_back(-5.);
	  pid_kScore.push_back(-5.);
	  pid_proScore.push_back(-5.);
	}
      } // TRACK FIT

      if ( parameters.runShowerFit && (parameters.trackScoreCut < 0. || trackScore < parameters.trackScoreCut) ) {
 //         std::cout << "I would have fit this as a shower..." << std::endl;
	//Begin Defining Shower Direction Through a PCA      
          std::cout<< "NEW SHOWER ALERT" << std::endl;
        shwrSliceId.push_back(sliceID);
        shwrClusterId.push_back(clusterID);

        std::cout << "SHOWER EVENT: " << entryIdx << " SHOWER CLUSTER: " << clusterID << std::endl;
        CartesianVector centroid(0.f, 0.f, 0.f);
        lar_content::LArPcaHelper::EigenVectors eigenVecs;
        lar_content::LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
        lar_content::LArPcaHelper::RunPca(caloHitList, centroid, eigenValues, eigenVecs);
	
	//Define directions to be positive
	const CartesianVector axisDirection(eigenVecs.at(0).GetZ() > 0.f ? eigenVecs.at(0) : eigenVecs.at(0) * -1.f);

    shwrCentroidX.push_back(centroid.GetX());
    shwrCentroidY.push_back(centroid.GetY());
    shwrCentroidZ.push_back(centroid.GetZ());

    std::cout << "Entries: " << shwrDirX.size() << std::endl;

	//Define Shower Length
	float NSigma = 2;
	float primaryEigenValue = eigenValues.GetX();
	float showerLength = std::sqrt(primaryEigenValue) * 2 * NSigma;
	
    shwrLen.push_back(showerLength);

	//Define the shower start position
	//loop over the caloHitList 
	float projection;
	std::map<float, CartesianVector> projectionMap;

	CartesianVector hitPosition(0.f, 0.f, 0.f);

    //Find projections for each hit along the primary axis and save them into a map from least to greatest	

	for (const CaloHit *const pCaloHit3D : caloHitList){

		projection = axisDirection.GetDotProduct(pCaloHit3D->GetPositionVector() - centroid); 
        hitPosition = pCaloHit3D->GetPositionVector();
		projectionMap.insert({projection, hitPosition});	
	}
    
    std::cout << "NUMBER OF HITS: " << caloHitList.size() << std::endl;


    // constants for looping through projection
	//Define a proximity radius and proximity threshold
	//

    CartesianVector showerStartHitPos(0.f, 0.f, 0.f);
    float showerStartHitProjectionValue(9999);

    CartesianVector hit_i_pos(0.f, 0.f, 0.f), hit_j_pos(0.f, 0.f, 0.f);
    float hit_i_proj(9999);

    int hitProximityRadius = 4;
    int proximityHitsCounter;
    float hit_i_j_dist;
    int proximityHitsThreshold = 10;

    for(auto iMapEntry: projectionMap){
        proximityHitsCounter = 0;
       // std::cout << "******************New Potentail Start Point*****************" << std::endl;
        hit_i_pos = iMapEntry.second;
        hit_i_proj = iMapEntry.first;
        //std::cout << "Hit Position Start of Loop: " << hit_i_pos << std::endl;
        //std::cout << "Hit Projection Start of Loop: " << hit_i_proj << std::endl;
        
        for(auto jMapEntry : projectionMap){
           // std::cout << "#### NEW HIT ####" << std::endl; 
            hit_j_pos = jMapEntry.second;
            //std::cout << "Hit Position: " << hit_j_pos << std::endl;
            if(hit_j_pos == hit_i_pos){
                continue;
            }

            hit_i_j_dist = std::sqrt(hit_i_pos.GetDistanceSquared(hit_j_pos));
        
            if(hit_i_j_dist <= hitProximityRadius){
            proximityHitsCounter ++;
         //   std::cout << "Distance: " << hit_i_j_dist << std::endl;
            
            }
            //std::cout << "Hits COunter: " << proximityHitsCounter << std::endl;

            if(proximityHitsCounter > proximityHitsThreshold){
                std::cout << "Passed! Counter is at: " << proximityHitsCounter << std::endl;
                showerStartHitPos = hit_i_pos;
                showerStartHitProjectionValue = hit_i_proj;
                break;
            }    
        }   
        if(proximityHitsCounter > proximityHitsThreshold){
            break;
        }
    }
    
    
    if(showerStartHitProjectionValue == 9999){
        showerStartHitPos = projectionMap.begin()->second;
        showerStartHitProjectionValue = projectionMap.begin()->first;
    }
  

    CartesianVector showerStartPosition = centroid + axisDirection*showerStartHitProjectionValue;
    std::cout << "Projected Shower Start: " << showerStartPosition << std::endl;
    std::cout << "Hit Shower Start: " << showerStartHitPos << std::endl;
    

      //end shower start point calculation
 	minProjection.push_back(showerStartHitProjectionValue);       	

	float showerStartX, showerStartY, showerStartZ;

	int showerStartLength = 5;
	int showerStartWidth = 4;


	showerStartX = showerStartHitPos.GetX();
	showerStartY = showerStartHitPos.GetY();
	showerStartZ = showerStartHitPos.GetZ();

    shwrStartHitPositionX.push_back(showerStartHitPos.GetX());
    shwrStartHitPositionY.push_back(showerStartHitPos.GetY());
    shwrStartHitPositionZ.push_back(showerStartHitPos.GetZ());

    //Define shower direction as a vector passing through both the start point and the centroid

    CartesianVector showerDirection = (centroid - showerStartHitPos);
    showerDirection = showerDirection.GetUnitVector();

    shwrPCAX.push_back(axisDirection.GetX());
    shwrPCAY.push_back(axisDirection.GetY());
    shwrPCAZ.push_back(axisDirection.GetZ());
    
    shwrDirX.push_back(showerDirection.GetX());
    shwrDirY.push_back(showerDirection.GetY());
    shwrDirZ.push_back(showerDirection.GetZ());

	std::cout << "Start Position is: " << showerStartPosition << std::endl;
    
        shwrStartX.push_back(showerStartX);
        shwrStartY.push_back(showerStartY);
        shwrStartZ.push_back(showerStartZ);

	//Define dE/dx of the shower
	
	float distanceFromShowerStart;	
	
	float hitPCAOpeningAngle, hitPositionAlongAxis, hitPositionFromAxis;
	float totalQ = 0;

	CartesianVector showerStartCurrentHit(0.f, 0.f, 0.f);
	CartesianVector showerStartPCAProjection(0.f, 0.f, 0.f);
	CaloHitList showerStartCaloHitList;	
	CartesianVector hitProjectedPosition(0.f, 0.f, 0.f);
	
    int caloHitIndex = 0;

    showerStartCaloHitList.clear();

	for(const CaloHit *const pShowerStartCaloHit3D : caloHitList){
		showerStartPCAProjection = centroid + ( showerDirection * showerStartHitProjectionValue);
        showerStartCurrentHit = pShowerStartCaloHit3D->GetPositionVector();
		totalQ += pShowerStartCaloHit3D->GetInputEnergy();	
		chargePerHit.push_back(pShowerStartCaloHit3D->GetInputEnergy());
        if(showerStartPCAProjection == showerStartCurrentHit){
            continue;
			}
        else{
        CartesianVector input = showerStartCurrentHit - showerStartPCAProjection;
            if(input.GetMagnitude() < 0.001){
               continue;
            }
		hitPCAOpeningAngle = showerDirection.GetOpeningAngle(showerStartCurrentHit - showerStartPCAProjection);
		distanceFromShowerStart = std::sqrt(showerStartCurrentHit.GetDistanceSquared(showerStartPCAProjection));
		hitPositionAlongAxis = distanceFromShowerStart * (std::cos(hitPCAOpeningAngle));
        hitPositionFromAxis = distanceFromShowerStart * std::sin(hitPCAOpeningAngle);		
        }

        if(hitPositionAlongAxis > 0 && hitPositionAlongAxis < showerStartLength && hitPositionFromAxis < showerStartWidth){
			showerStartCaloHitList.push_back(pShowerStartCaloHit3D);
            shwrStartPointsRecoId.push_back(clusterID);
            
            chargePerHitStartPoints.push_back(pShowerStartCaloHit3D->GetInputEnergy());
	       //Add positions of hits to a branch to look at later
			shwrStartPointsX.push_back(showerStartCurrentHit.GetX());
            shwrStartPointsY.push_back(showerStartCurrentHit.GetY());
            shwrStartPointsZ.push_back(showerStartCurrentHit.GetZ());
	        caloHitIndex++;		
            } 

		else{
			continue;
			}	

	}
    


   	std::cout << "Reco Hit List Size " << showerStartCaloHitList.size() << std::endl; 
	shwrTotalE.push_back(totalQ);

	if(showerStartCaloHitList.size() < 2){
		std::cout<< "Reco Hit List Too Small. Skipping fit. Event is: " << entryIdx << std::endl;
        
		startTrkSlidingFitDirX.push_back(-9999.);
        startTrkSlidingFitDirY.push_back(-9999.);
        startTrkSlidingFitDirZ.push_back(-9999.);


        shwrdEdx.push_back(-9999.);
		pitchValue.push_back(-9999.);
        medianQ.push_back(-9999.);
        continue;
		}

	else{
//track fit of beginning points
	float showerSlidingFitHalfWindow = 20;
	std::vector<int> showerFitIndexVector;

	lar_content::LArTrackStateVector showerFitTrackStateVector;
	
	lar_content::LArPfoHelper::GetSlidingFitTrajectory( &showerStartCaloHitList, vertexVector, showerSlidingFitHalfWindow, parameters.pixelPitch, showerFitTrackStateVector, &showerFitIndexVector, true );

if(showerFitTrackStateVector.size() < 1){

std::cout << "TRACK STATE FAILED " << " Event is: " << entryIdx << "Slice: "<< sliceID << "Cluster: " << clusterID << std::endl;

        startTrkSlidingFitDirX.push_back(-9999.);
        startTrkSlidingFitDirY.push_back(-9999.);
        startTrkSlidingFitDirZ.push_back(-9999.);
	
	shwrdEdx.push_back(-9999.);
	pitchValue.push_back(-9999.);
    medianQ.push_back(-9999.);	
}

else{

	startTrkSlidingFitDirX.push_back(showerFitTrackStateVector.front().GetDirection().GetX());
	startTrkSlidingFitDirY.push_back(showerFitTrackStateVector.front().GetDirection().GetY());
	startTrkSlidingFitDirZ.push_back(showerFitTrackStateVector.front().GetDirection().GetZ());


	//dQ/dx
	double hitQ_Shower;
	float cosgamma, projectedShowerDir, pitch;	
	std::vector<double> Q_ShowerVector;

	for(unsigned int iTrackHit = 0; iTrackHit < showerFitTrackStateVector.size(); ++iTrackHit){
		const lar_content::LArTrackState& showerFitTrackState = showerFitTrackStateVector.at(iTrackHit);
		hitQ_Shower = showerFitTrackState.GetCaloHit()->GetInputEnergy();
     //   std::cout << "Q Hit Energy" << hitQ_Shower << std::endl;
		Q_ShowerVector.push_back(hitQ_Shower);
		}

	projectedShowerDir = std::sqrt(std::pow(showerFitTrackStateVector.front().GetDirection().GetY(), 2) +
            std::pow(showerFitTrackStateVector.front().GetDirection().GetZ(), 2));
    
    //std::cout << "Projected Shower Dir" << projectedShowerDir << std::endl;
 
if(projectedShowerDir == 0){
	std::cout << "Pitch is calculated to be inf. dEdx calculation will be skipped." << std::endl;
	shwrdEdx.push_back(-9999.);
	pitchValue.push_back(-9999.);
    medianQ.push_back(-9999.);
}
else{
	cosgamma = std::abs(showerFitTrackStateVector.front().GetDirection().GetZ() / projectedShowerDir);
		if(cosgamma == 0){
			std::cout << "Pitch is calculated to be inf. dEdx calculation will be skipped." << std::endl;
			shwrdEdx.push_back(-9999.);
		    pitchValue.push_back(-9999.);
            medianQ.push_back(-9999.);
        }	
		else{
			pitch = parameters.pixelPitch/cosgamma;
            if(pitch == 0){
				std::cout << "Pitch is calculated to be 0. dEdx calculation will be skipped." << std::endl;
				shwrdEdx.push_back(-9999.);
			    pitchValue.push_back(-9999.);
                medianQ.push_back(-9999.);
            }
			else{
				double Q = TMath::Median(Q_ShowerVector.size(), &Q_ShowerVector[0]);
                double dQdx =  Q / pitch;
                std::cout << "Uncalibrated dQdx = " << dQdx << std::endl;
                double dQdx_Calibrated = dQdx*(1000);
                std::cout << "Calibrated dQdx = " << dQdx_Calibrated << std::endl;
				double dEdx = dEdxWithRecombination(parameters, dQdx_Calibrated);
                std::cout << "Total dEdx = " << dEdx << std::endl;
                shwrdEdx.push_back(dEdx);
				pitchValue.push_back(pitch);
                medianQ.push_back(Q);
            }	
		}

    }
}

}

	 } // SHOWER FIT


}

    if ( parameters.runTrackFit ) {
      fOut.FillTrackBranches(trkStartX,trkStartY,trkStartZ,trkStartDirX,trkStartDirY,trkStartDirZ,trkEndX,trkEndY,trkEndZ,trkEndDirX,trkEndDirY,trkEndDirZ,trkLen,trkContained,
			     trk_KEFromLength_muon,trk_KEFromLength_proton,trk_pFromLength_muon,trk_pFromLength_proton);
      fOut.FillTrackCaloBranches(trackFitTrackCaloE,trackFitVisE,trackFitSliceId,trackFitPfoId,trackFitX,trackFitY,trackFitZ,trackFitQ,trackFitRR,trackFitdx,trackFitdQdx,trackFitdEdx);
      if ( parameters.fShouldRunPID ) fOut.FillTrackPID(pid_pdg,pid_ndf,pid_muScore,pid_piScore,pid_kScore,pid_proScore);
    }
   if ( parameters.runShowerFit){
     fOut.FillShowerBranches(shwrCentroidX,shwrCentroidY,shwrCentroidZ,shwrStartX, shwrStartY, shwrStartZ, shwrDirX, shwrDirY, shwrDirZ, shwrLen, 
             shwrSliceId, shwrClusterId, startTrkSlidingFitDirX, startTrkSlidingFitDirY, startTrkSlidingFitDirZ,
             shwrStartPointsX, shwrStartPointsY, shwrStartPointsZ, shwrdEdx, shwrTotalE, shwrStartPointsRecoId, minProjection, 
             medianQ,chargePerHit, chargePerHitStartPoints, pitchValue,
             shwrPCAX, shwrPCAY, shwrPCAZ,
             shwrStartHitPositionX, shwrStartHitPositionY, shwrStartHitPositionZ);
}


    // Write our branches to the output tree
    fOut.WriteToFile();
  } // loop entries

  // Close our output file
  fOut.CloseFile();

  return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ParseCommandLine(int argc, char *argv[], ParameterStruct &parameters)
{
    if (1 == argc)
        return PrintOptions();

    int cOpt(0);

    bool hasInputFile=false;
    bool hasXmlFile=false;

    while ((cOpt = getopt(argc, argv, "x:f:o:h")) != -1)
    {
      switch (cOpt)
      {
        case 'x':
	  parameters.xmlName = optarg;
	  hasXmlFile = true;
	  break;
        case 'f':
	  parameters.fileName = optarg;
	  hasInputFile = true;
	  break;
        case 'o':
	  parameters.outfileName = optarg;
	  break;
        case 'h':
        default:
	  return PrintOptions();
      }
    }

    bool passed=hasXmlFile && hasInputFile;
    if(!passed)
    {
        return PrintOptions();
    }
    return passed;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ReadSettings(ParameterStruct &parameters)
{
  TiXmlDocument xmlDocument(parameters.xmlName);

  if (!xmlDocument.LoadFile()) {
    std::cout << "XML document (" << parameters.xmlName << ") not loaded. Returning." << std::endl;
    return false;
  }

  const TiXmlHandle xmlDocumentHandle(&xmlDocument);
  const TiXmlHandle xmlHandle(TiXmlHandle(xmlDocumentHandle.FirstChildElement().Element()));

  try
  {
    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShouldRunTrackFit", parameters.runTrackFit) );
    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShouldRunShowerFit", parameters.runShowerFit) );
    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrackScoreCut", parameters.trackScoreCut) );
    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PixelPitch", parameters.pixelPitch) );

    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShouldApplyHitThreshold", parameters.applyThreshold) );
    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ChargeThreshold", parameters.thresholdVal) );

    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShouldVoxelizeZ", parameters.voxelizeZ) );
    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VoxelHalfWidthZ", parameters.voxelZHW) );
    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UseVoxelizedStartStop", parameters.useVoxelizedStartStop) );

    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Detector", parameters.fDetector) );

    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShouldCorrectLifetime", parameters.fShouldCorrectLifetime) );
    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ElectronLifetime", parameters.fElectronLifetime) );
    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ElectronDriftSpeed", parameters.fElectronDriftSpeed) );

    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShouldCorrectRecombination", parameters.fShouldCorrectRecomb) );
    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FlowStyleRecombination", parameters.fFlowStyleRecombination) );
    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "BoxRecombination", parameters.fBoxRecombination) );
    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "BirksRecombination", parameters.fBirksRecombination) );
    // Check if we have set to do BOTH Box & Birks if correcting lifetime, or if we have set to correct lifetime and chosen NEITHER
    if ( parameters.fShouldCorrectRecomb ) {
      if ( parameters.fBoxRecombination && parameters.fBirksRecombination ) {
	std::cout << "BOTH Box and Birks are true. Please set one to false." << std::endl;
	return false;
      }
      if ( !parameters.fBoxRecombination && !parameters.fBirksRecombination ) {
	std::cout << "BOTH Box and Birks are false. Please set one to true." << std::endl;
	return false;
      }
    }

    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Density", parameters.fDensity) );
    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EField", parameters.fEField) );
    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "BoxBeta", parameters.fBoxBeta) );
    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "BoxAlpha", parameters.fBoxAlpha) );
    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "BirksA", parameters.fBirksA) );
    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "BirksK", parameters.fBirksK) );
    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ApplyCalibrationFudgeFactor", parameters.fApplyCalibrationFudgeFactor) );
    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ApplyCalibrationFudgeFactorPID", parameters.fApplyCalibrationFudgeFactor_PID) );
    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CalibrationFudgeFactor", parameters.fCalibrationFudgeFactor) );

    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShouldRunPID", parameters.fShouldRunPID) );
    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PIDAlgChi2PID", parameters.fPIDAlgChi2PID) );
    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Chi2RestrictDX", parameters.fChi2RestrictDX) );
    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Chi2RestrictDXLo", parameters.fChi2RestrictDXLo) );
    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Chi2RestrictDXHi", parameters.fChi2RestrictDXHi) );
    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Chi2RestrictDEDXLo", parameters.fChi2RestrictDEDXLo) );
    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "dEdxResTempFile", parameters.fdEdxResTempFile) );

    if ( parameters.fShouldRunPID && parameters.fPIDAlgChi2PID ) {
      // Following scheme as in LArSoft larana/ParticleIdentification/Chi2PIDAlg.cxx
      TFile *tempFile = TFile::Open(parameters.fdEdxResTempFile.c_str());
      parameters.templatesdEdxRR["muon"] = (TProfile*)tempFile->Get("dedx_range_mu");
      parameters.templatesdEdxRR["pion"] = (TProfile*)tempFile->Get("dedx_range_pi");
      parameters.templatesdEdxRR["proton"] = (TProfile*)tempFile->Get("dedx_range_pro");
      parameters.templatesdEdxRR["kaon"] = (TProfile*)tempFile->Get("dedx_range_ka");
      try {
	if ( parameters.templatesdEdxRR.at("muon")->GetNbinsX() < 1 &&
	     parameters.templatesdEdxRR.at("pion")->GetNbinsX() < 1 &&
	     parameters.templatesdEdxRR.at("kaon")->GetNbinsX() < 1 &&
	     parameters.templatesdEdxRR.at("proton")->GetNbinsX() < 1 ) {
	  std::cout << "Not enough data in dEdx vs RR templates. Returning" << std::endl; 
	  return false;
	}
      }
      catch (StatusCodeException &statusCodeException) {
	std::cout << "Issue loading dEdx vs RR templates. Returning" << std::endl;
	return false;
      }
    }

    PANDORA_RETURN_RESULT_IF_AND_IF( pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Verbosity", parameters.verbosity) );
  }
  catch (StatusCodeException &statusCodeException)
  {
    std::cout << "Failed to initialized parameters in the XML file. Status code " << statusCodeException.GetStatusCode() << ". Returning." << std::endl;
    return false;
  }

  return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool PrintOptions()
{
    std::cout << std::endl
              << "./bin/PandoraOuterface -f [file name] -o [out name] -p [pixel pitch] -c [track score cut] -t -s" << std::endl
	      << "    -t = run track fit" << std::endl
	      << "    -s = run shower fit"
              << std::endl;

    return false;
}

} // namespace lar_nd_postreco
