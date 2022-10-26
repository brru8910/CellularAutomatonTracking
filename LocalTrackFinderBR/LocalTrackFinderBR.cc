/**
  \file
  Declaration of LocalTrackFinderBR

  \author B. Rumberger
  \version $Id$
  \date 31 Oct 2017
*/

#include "LocalTrackFinderBR.h"
#include <det/Detector.h>
#include <det/Target.h>
#include <det/TPC.h>
#include <det/TPCChamber.h>
#include <evt/Event.h>
#include <evt/ProcEvent.h>
#include <evt/IndexedObjectLinker.h>
#include <fwk/CentralConfig.h>
#include <modutils/MagneticTrackFitter.h>
#include <modutils/TrackAutopsy.h>
#include <utl/BranchIterator.h>
#include <utl/ErrorLogger.h>
#include <utl/GeometryUtilities.h>
#include <utl/Verbosity.h>

using namespace evt;
using namespace fwk;
using namespace std;
using namespace utl;

namespace LocalTrackFinderBR {
  
  VModule::EResultFlag LocalTrackFinderBR::Init() {
    CentralConfig& cc = CentralConfig::GetInstance();
    
    utl::Branch topBranch = cc.GetTopBranch("LocalTrackFinderBR");
    InitVerbosity(topBranch);
    
    //Get module-wide parameters.
    topBranch.GetChild("binDensityFactor").GetData(fBinDensityFactor);
    topBranch.GetChild("minKalmanFilterStepSize").GetData(fMinKalmanFilterStepSize);
    topBranch.GetChild("maxKalmanFilterStepSize").GetData(fMaxKalmanFilterStepSize);
      
    vector<string> detectorList;

    utl::Branch detectorSetupBranch = topBranch.GetChild("detectorSetupList");
    for (utl::BranchIterator setupIt = detectorSetupBranch.ChildrenBegin(),
	   setupEnd = detectorSetupBranch.ChildrenEnd(); setupIt != setupEnd; ++setupIt) {
      //Read in setup for each detector from .xml file. If the detector doesn't have a 
      //"detectorSetup" entry, it will not be treated.
      utl::Branch setup = *setupIt;
      modutils::DetectorSetup detectorSetup;
      setup.GetChild("xDisplacementTolerance").GetData(detectorSetup.fXDisplacementTolerance);
      setup.GetChild("yDisplacementTolerance").GetData(detectorSetup.fYDisplacementTolerance);
      setup.GetChild("maxAngleChangeXZ").GetData(detectorSetup.fMaxAngleChangeXZ);
      setup.GetChild("maxAngleChangeYZ").GetData(detectorSetup.fMaxAngleChangeYZ);
      setup.GetChild("xParameterWeight").GetData(detectorSetup.fXParameterWeight);
      setup.GetChild("yParameterWeight").GetData(detectorSetup.fYParameterWeight);
      setup.GetChild("aParameterWeight").GetData(detectorSetup.fAParameterWeight);
      setup.GetChild("bParameterWeight").GetData(detectorSetup.fBParameterWeight);
      setup.GetChild("qOverPParameterWeight").GetData(detectorSetup.fQOverPParameterWeight);
      setup.GetChild("minimumClustersForTrackSeed").GetData(detectorSetup.fMinimumClustersForTrackSeed);
      setup.GetChild("minimumClustersOnLocalTrack").GetData(detectorSetup.fMinimumClustersOnLocalTrack);
      setup.GetChild("maxMergeMetric").GetData(detectorSetup.fMaxMergeMetric);
      setup.GetChild("maxAverageResidual").GetData(detectorSetup.fMaxAverageResidual);

      //Get detector ID and add to list.
      string detectorName = (setup.GetAttributes())["detectorName"];
      //We only need to get one of the MTPCs, since the cluster planes span both.
      if (detectorName == "MTPCs")
	detectorName = "MTPCL";
      det::Const::EId globalId = det::Const::GetId(detectorName);
      if (globalId == det::Const::eUnknown) {
        ostringstream err;
        err << "Unknown detector " << detectorName;
        ERROR(err);
        return eFailure;
      }
      detectorSetup.fGlobalId = globalId;

      //Get tracking algorithm to use in this detector.
      string fitterName = setup.GetChild("fitterToUse").GetDataString();
      if (fitterName == "kalmanFilter")
	detectorSetup.fFitterType = modutils::eKalmanFilter;
      else if (fitterName == "straightTrackFitter")
	detectorSetup.fFitterType = modutils::eStraightTrackFitter;
      else if (fitterName == "parabolicTrackFitter")
	detectorSetup.fFitterType = modutils::eParabolicTrackFitter;
      else
	throw utl::DoesNotComputeException("Fitter type not correctly chosen! Check "
					   "LocalTrackFinderBR.xml file!");
      
      fDetectorSetups.push_back(detectorSetup);
    }
    
    ostringstream info;
    info << "Local tracker initialization.";

    if ( fVerbosity >= utl::Verbosity::eTerse )
      INFO(info);

    return eSuccess;
  }
  
  /// Process function for LocalTrackFinderBR
  VModule::EResultFlag LocalTrackFinderBR::Process(Event& event, const utl::AttributeMap& /*attr*/) {

    //Get necessary event and detector instances.
    evt::RecEvent& recEvent = event.GetRecEvent();
    const evt::EventHeader& eventHeader = event.GetEventHeader();
    const unsigned int runNumber = eventHeader.GetRunNumber();
    const unsigned int eventNumber = eventHeader.GetId();

    const det::Detector& detector = det::Detector::GetInstance();
    fDetectorCS = detector.GetDetectorCoordinateSystem();
    const det::TPC& detTPC = detector.GetTPC();
    //Initialize Kalman Fiter infrastructure.
    fKalmanFilter.Init(recEvent);
    fKalmanFilter.SetMinimumStepSize(fMinKalmanFilterStepSize);
    fKalmanFilter.SetMaximumStepSize(fMaxKalmanFilterStepSize);
    //Initialize least squares fitters..
    fStraightTrackFitter.Init(recEvent,fDetectorCS);
    fParabolicTrackFitter.Init(recEvent,fDetectorCS);

    if ( fVerbosity >= utl::Verbosity::eAnnoying ) { 
      ostringstream msg;
      msg << "Searching for local tracks. Run " << runNumber
	  << ", event " << eventNumber;
      INFO(msg);
    }

    //Initialize grid objects on all TPC detection planes. Clear last event's cluster positions.
    fClusterPositionsByIndex.clear();
    fClusterZPlanesByIndex.clear();
    fClusterErrorsByIndex.clear();
    InitializeDetectionPlaneGrids(event,detTPC);
    //Get proc event with detection plane grids initialized.
    evt::proc::Reconstruction& procReconstruction = event.GetProcEvent().GetReconstruction();
    
    for (std::vector<modutils::DetectorSetup>::iterator setupIt = fDetectorSetups.begin(),
	   detectorEnd = fDetectorSetups.end(); setupIt != detectorEnd; ++setupIt) {

      det::Const::EId detectorId = (*setupIt).fGlobalId;

      if (detectorId == det::Const::eUnknown) {
	std::ostringstream err;
	err << "Global detector ID= " << detectorId
	    << " has no corresponding detector.";
	ERROR(err.str());
	throw utl::NonExistentComponentException(err.str());
      }
      //TPC-specific. Change at some point to include vertex detector/future detectors.
      const det::TPCConst::EId& chamberId =
	det::TPCConst::GetId(det::Const::GetName(detectorId));
      //Make sure this chamber is included in the current detector configuration.
      if (!detTPC.HasChamber(chamberId))
	continue;
      const det::TPCChamber& chamber = detTPC.GetChamber(chamberId);
      const det::TPCChamber& detMTPCR = detTPC.GetChamber(det::TPCConst::eMTPCR);
      modutils::DetectorSetup setup = *setupIt;
      //Set detector z-plane ID range.
      setup.fZPlaneIdRange = chamber.GetZPlaneIdRange();

      //Initialize Cellular Automaton.
      fCellularAutomaton.Initialize(setup);
      //Initialize Tracklet Constructor.
      fTrackletConstructor.Initialize(setup);
      //Initialize Track Segment Fitter.
      fTrackSegmentFitter.Initialize(setup,
				     detMTPCR);
      fTrackSegmentFitter.SetKalmanFilter(fKalmanFilter);
      fTrackSegmentFitter.SetStraightTrackFitter(fStraightTrackFitter);
      fTrackSegmentFitter.SetParabolicTrackFitter(fParabolicTrackFitter);
      //Initialize Track Segment Extender.
      fTrackSegmentExtender.Initialize(setup,
				       fClusterPositionsByIndex,
				       fClusterZPlanesByIndex,
				       detMTPCR);
      fTrackSegmentExtender.SetKalmanFilter(fKalmanFilter);
      fTrackSegmentExtender.SetStraightTrackFitter(fStraightTrackFitter);
      fTrackSegmentExtender.SetParabolicTrackFitter(fParabolicTrackFitter);

      //Construct tracklets between adjacent detection planes first (separated by just one z-plane).
      //This is run iteratively, with the allowed search space increasing with each iteration.
      //This keeps combinatorics down.
      const unsigned int nIterations = 5;
      for (unsigned int i = 1; i <= nIterations; i++) {
      	fTrackletConstructor.ConstructTracklets(procReconstruction,
      						fClusterPositionsByIndex,
      						fCellularAutomaton,
      						i);
      }
      //Filter tracklets.
      fCellularAutomaton.ApplyTrackletPointingRule();
      //Calculate clusters represeting contiguous tracks.
      fCellularAutomaton.FindTrackCandidates();
      // Get lists of clusters corresponding to local tracks. These will be fit with the Kalman Filter.
      modutils::TrackCandidates& trackCandidates = fCellularAutomaton.GetTrackCandidates();
      // //Write tracklets. For testing.
      // WriteTracklets(event,trackCandidates,setup,chamber);
      //Fit entire track segment using desired track model and filter clusters based on result.
      fTrackSegmentFitter.FitTrackCandidates(procReconstruction,
      					     trackCandidates,
      					     fClusterPositionsByIndex,
      					     fClusterZPlanesByIndex,
      					     chamber);
      modutils::LocalTrackSegments& localTrackSegments = fTrackSegmentFitter.GetTrackSegments();
      modutils::ClusterList& clustersOnTracks = fTrackSegmentFitter.GetUsedClusters();

      fTrackSegmentExtender.MergeTrackSegments(procReconstruction,
      					       localTrackSegments,
      					       chamber);
      //Pick up lone clusters left out due to noise or detector inefficiencies.
      fTrackSegmentExtender.PerformClusterPickup(procReconstruction,
      						 localTrackSegments,
      						 clustersOnTracks,
      						 chamber);

      //Store tracks found in RecEvent.
      WriteTracks(event,localTrackSegments,setup);

    } //End detector setups loop.
 

    return eSuccess;
  }


  /// Finish function for LocalTrackFinderBR
  VModule::EResultFlag LocalTrackFinderBR::Finish() {
    ostringstream info;
    info << "Local tracking Finished." << "\n";
    INFO(info);
    
    return eSuccess;
  }

  //Function for adding local tracks to recEvent.
  void LocalTrackFinderBR::WriteTracks(evt::Event& event,
				       modutils::LocalTrackSegments& localTrackSegments,
				       const modutils::DetectorSetup& setup) {
    
    //Count number of tracks written.
    unsigned int nTracks = 0;
    
    //FIXME!! Get track const...this is very bad! Need a getter function that returns
    //trackConst from detector ID!!
    const det::TPCConst::EId& chamberId =
      det::TPCConst::GetId(det::Const::GetName(setup.fGlobalId));
    rec::TrackConst::ETPC detector = rec::TrackConst::eAll;
    if (chamberId == det::TPCConst::eVTPC1) detector = rec::TrackConst::eVTPC1;
    if (chamberId == det::TPCConst::eVTPC2) detector = rec::TrackConst::eVTPC2;
    if (chamberId == det::TPCConst::eGTPC)  detector = rec::TrackConst::eGTPC;
    if (chamberId == det::TPCConst::eMTPCL) detector = rec::TrackConst::eMTPC;
    if (chamberId == det::TPCConst::eMTPCR) detector = rec::TrackConst::eMTPC;
    if (chamberId == det::TPCConst::eFTPC1) detector = rec::TrackConst::eFTPC1;
    if (chamberId == det::TPCConst::eFTPC2) detector = rec::TrackConst::eFTPC2;
    if (chamberId == det::TPCConst::eFTPC3) detector = rec::TrackConst::eFTPC3;

    //Loop through saved local tracks.
    for (modutils::LocalTrackSegments::iterator trackSegmentsIt = localTrackSegments.begin(),
	   trackSegmentsEnd = localTrackSegments.end();
	 trackSegmentsIt != trackSegmentsEnd; ++trackSegmentsIt) {
      
      modutils::LocalTrackSegment& trackSegment = *trackSegmentsIt;
      //Sort clusters by z-position (just in case).
      trackSegment.SortClusters();
      
      //Get fit parameters.
      const modutils::KalmanTrackParameters& trackParameters =
	trackSegment.GetTrackParameters();

      const utl::Point& lastClusterPosition =
	fClusterPositionsByIndex[trackSegment.GetLastClusterIndex()];

      //Get rec event and proc event.
      evt::RecEvent& recEvent = event.GetRecEvent();
      evt::ProcEvent& procEvent = event.GetProcEvent();

      //Create track in recEvent.
      evt::rec::Track& track = recEvent.Make<rec::Track>();

      //Set "good" status.
      track.SetStatus(0);

      const utl::Point& momentumPoint = trackParameters.GetPosition();
      const utl::Vector& momentum = trackParameters.GetMomentum();
      const int charge = trackParameters.GetCharge();
      const unsigned int nClusters = trackSegment.GetNumberOfClusters();

      //Store track information.
      track.SetMomentumPoint(momentumPoint);
      track.SetMomentum(momentum);
      track.SetCharge(charge);
      track.SetChi2(trackParameters.GetChi2());
      track.SetNdf(trackParameters.GetNdf());
      track.SetCovarianceMatrix(trackParameters.GetCovarianceMatrix(utl::TrackParameters::eNA61));
      track.SetNumberOfClusters(detector,nClusters);
      track.SetNumberOfClusters(rec::TrackConst::eAll,nClusters);
      track.SetNumberOfFitClusters(detector,nClusters);
      track.SetNumberOfFitClusters(rec::TrackConst::eAll,nClusters);
      track.SetFirstPointOnTrack(momentumPoint);
      track.SetLastPointOnTrack(lastClusterPosition);

      //Loop through clusters and link to parent track.
      for (proc::ClusterVector::const_iterator clusterIt = trackSegment.GetClusters().begin(),
	     clusterEnd = trackSegment.GetClusters().end(); clusterIt != clusterEnd; clusterIt++) {
	evt::rec::Cluster& cluster = recEvent.Get(*clusterIt);
	try {
	  evt::IndexedObjectLinker::LinkDaughterToParent(cluster, track);
	}
	catch (exception& exc) {
	  ERROR("Couldn't link cluster to track! INVESTIGATE THIS!!");
	  continue;
	}
      }

      //Add tracks by chamber to ProcEvent detection plane grids.
      proc::Reconstruction& procReconstruction = procEvent.GetReconstruction();
      procReconstruction.AddTrackByDetector(track.GetIndex(),setup.fGlobalId);
      proc::DetectionPlaneGrid& grid =
	procReconstruction.GetDetectionPlaneGrid(trackSegment.GetFirstZPlaneId());
      grid.AddTrackToGrid(track.GetIndex(),track.GetMomentumPoint());

      nTracks++;
    }
    //Include some results info.
    if ( fVerbosity >= utl::Verbosity::eTerse ) {
      ostringstream info;
      const string detectorName = (setup.fGlobalId == det::Const::eMTPCL || 
				   setup.fGlobalId == det::Const::eMTPCR  ) ? 
	"MTPCs" : det::Const::GetName(setup.fGlobalId);
      info << "Found " << nTracks << " local tracks in " << detectorName << ".";
      INFO(info);
    }
  }
  
  //Temporary function for writing tracklets to RecEvent (before fitting).
  //This is to test the tracklet construction efficiency.
  void LocalTrackFinderBR::WriteTracklets(evt::Event& event,
					  modutils::TrackCandidates& trackCandidates,
					  const modutils::DetectorSetup& setup,
					  const det::TPCChamber& chamber) {
    //Loop through track candidates -- as of now, these are just connected clusters.
    for (modutils::TrackCandidates::iterator trackCandidatesIt = trackCandidates.begin(),
	   trackCandidatesEnd = trackCandidates.end();
	 trackCandidatesIt != trackCandidatesEnd; ++trackCandidatesIt) {

      //Get candidate.      
      modutils::TrackCandidate& candidate = *trackCandidatesIt;

      //Order clusters by z-position.
      candidate.SortClusters();

      unsigned int nClustersOnTracklet = candidate.GetTrackClusters().size();
      if ( nClustersOnTracklet < setup.fMinimumClustersOnLocalTrack )
	continue;
 
      //Get rec event and proc event.
      evt::RecEvent& recEvent = event.GetRecEvent();

      //Create track in recEvent.
      evt::rec::Track& track = recEvent.Make<rec::Track>();

      //Set first and last cluster positions.
      const evt::rec::Cluster& firstCluster = recEvent.Get(candidate.GetTrackClusters().front());
      const evt::rec::Cluster& lastCluster = recEvent.Get(candidate.GetTrackClusters().back());

      //Get first z-plane.
      const det::TPCPadrow& padrow =
	chamber.GetSector(firstCluster.GetSectorNumber()).GetPadrow(firstCluster.GetPadrowNumber());
      const utl::Plane& zPlane = chamber.GetZPlane(padrow.GetZPlaneId());
      
      //Set "good" status.
      track.SetStatus(0);

      modutils::StraightTrackFitResult fitResult;
      modutils::KalmanTrackParameters trackParameters;
      fStraightTrackFitter.FitStraightLine(candidate.GetTrackClusters(), fitResult);
      //Extrapolate fit to first cluster z-coordinate.
      //(Tracklets are not used for analysis, no need) to go crazy here.
      trackParameters.SetStraightFitParameters(fitResult,zPlane);
      track.SetMomentumPoint(trackParameters.GetPosition());
      track.SetMomentum(trackParameters.GetMomentum());
      track.SetFirstPointOnTrack(firstCluster.GetPosition());
      track.SetLastPointOnTrack(lastCluster.GetPosition());

      //Loop through clusters and link to parent track.
      for (proc::ClusterVector::const_iterator clusterIt = candidate.GetTrackClusters().begin(),
	     clusterEnd = candidate.GetTrackClusters().end(); clusterIt != clusterEnd; ++clusterIt) {
	evt::rec::Cluster& cluster = recEvent.Get(*clusterIt);
	try {
	  evt::IndexedObjectLinker::LinkDaughterToParent(cluster, track);
	}
	catch (exception& exc) {
	  ERROR("Couldn't link cluster to track! INVESTIGATE THIS!!");
	  continue;
	}
      }
      track.SetNumberOfClusters(rec::TrackConst::eAll,nClustersOnTracklet);
    }
  }

  //Function for initializing grid objects on all TPC detection planes.
  void LocalTrackFinderBR::InitializeDetectionPlaneGrids(evt::Event& event, const det::TPC& detTPC) {
    //Get RecEvent, procReconstruction, ZPlaneIdRange, and proc::Clusters.
    const RecEvent& recEvent = event.GetRecEvent();
    const proc::TPC& procTPC = event.GetProcEvent().GetTPC();
    proc::Reconstruction& procReconstruction = event.GetProcEvent().GetReconstruction();
    const det::TPCChamber::ZPlaneIdRange& zPlaneIdRange = detTPC.GetZPlaneIdRange();
    //Get the TPC cluster indices from ProcEvent ordered along detection planes.
    const proc::tpc::Clusters& clusters = procTPC.GetClusters();
    //Loop through all detection planes in chamber.
    for (unsigned int zPlaneId = zPlaneIdRange.GetBegin(), zPlaneEnd = zPlaneIdRange.GetEnd();
	 zPlaneId < zPlaneEnd; ++zPlaneId) {
      //Get z-plane's geometric plane z.
      const double zPlaneZ = detTPC.GetZPlane(zPlaneId).GetAnchor().GetZ();
      //Create container for quick access to cluster position by index.
      proc::ClusterPositionsList clusterPositionsByIndex;
      const proc::tpc::Clusters::ClusterIndices& clusterIndices =
	clusters.GetClusterIndicesOnTPCDetectionPlane(zPlaneId);      
      for ( proc::tpc::Clusters::ClusterIndices::const_iterator clusterIt = clusterIndices.begin(),
	      clusterEnd = clusterIndices.end(); clusterIt != clusterEnd; ++clusterIt ) {
	//Get the cluster. Use some protection.
	rec::Cluster cluster;
	try {
	  cluster = recEvent.Get(*clusterIt);
	}
	catch (exception& exc) {
	  ERROR("Problem getting cluster!! Did you fill the proc::tpc::Cluster containers correctly?!");
	  continue;
	}
	//Get cluster position.
	const utl::Point clusterPosition = cluster.GetPosition();
 	//Fill the container.
	modutils::ClusterIndexAndPositionPair positionAndIndexPair((int)(*clusterIt),clusterPosition);
	//Make a utl::Point out of the cluster errors and store.
	const utl::Point clusterErrors(cluster.GetPositionUncertainty(rec::ClusterConst::eX),
				       cluster.GetPositionUncertainty(rec::ClusterConst::eY),0);
	modutils::ClusterIndexAndPositionPair errorsAndIndexPair((int)(*clusterIt),clusterErrors);
	clusterPositionsByIndex.push_back(positionAndIndexPair);
	//Record z-plane ID by cluster index.
	modutils::ClusterIndexAndZPlanePair zPlaneAndIndexPair((int)(*clusterIt),(int)zPlaneId);
	//Also fill container for quick cluster position access throughout module.
	fClusterPositionsByIndex.insert(positionAndIndexPair);
	fClusterZPlanesByIndex.insert(zPlaneAndIndexPair);
	fClusterErrorsByIndex.insert(errorsAndIndexPair);
      }
      //Don't do anything if this grid has previously been initialized.
      if (procReconstruction.IsDetectionPlaneGridInitialized(zPlaneId)) continue;
      proc::DetectionPlaneGrid grid(zPlaneId,zPlaneZ,fBinDensityFactor);
      grid.SetClusterPositionsList(clusterPositionsByIndex);
      grid.Initialize();
      procReconstruction.AddDetectionPlaneGrid(zPlaneId,grid);
    }
  }

}

