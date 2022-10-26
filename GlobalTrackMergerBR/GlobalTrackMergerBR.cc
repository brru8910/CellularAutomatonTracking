/**
   \file
   Declaration of GlobalTrackMergerBR

   \author B. Rumberger
   \version $Id$
   \date 5 March 2018
*/

#include "GlobalTrackMergerBR.h"
#include <det/ConstantTranslation.h>
#include <det/Detector.h>
#include <det/TPCChamber.h>
#include <evt/Event.h>
#include <evt/IndexedObjectLinker.h>
#include <fwk/CentralConfig.h>
#include <utl/BranchIterator.h>
#include "../../Simulation/TrackingEfficiencyTestsBR/LocalTrackingEfficiencyTestsBR.h"

using namespace evt;
using namespace fwk;
using namespace std;
using namespace utl;

namespace GlobalTrackMergerBR {
  
  VModule::EResultFlag GlobalTrackMergerBR::Init() {
    CentralConfig& cc = CentralConfig::GetInstance();
    
    utl::Branch topBranch = cc.GetTopBranch("GlobalTrackMergerBR");
    InitVerbosity(topBranch);
    
    vector<string> detectorList;

    //Get detector-wide constants.
    topBranch.GetChild("minKalmanFilterStepSize").GetData(fMinKalmanFilterStepSize);
    topBranch.GetChild("maxKalmanFilterStepSize").GetData(fMaxKalmanFilterStepSize);
    topBranch.GetChild("maxMergeMetric").GetData(fMaxMergeMetric);
    topBranch.GetChild("maxAverageResidual").GetData(fMaxAverageResidual);
    utl::Branch detectorSetupBranch = topBranch.GetChild("detectorSetupList");
    for (utl::BranchIterator setupIt = detectorSetupBranch.ChildrenBegin(),
	   setupEnd = detectorSetupBranch.ChildrenEnd(); setupIt != setupEnd; ++setupIt) {
      //Read in setup for each detector from .xml file. If the detector doesn't have a 
      //"detectorSetup" entry, it will not be treated.
      utl::Branch setup = *setupIt;
      modutils::DetectorSetup detectorSetup;
      setup.GetChild("xParameterWeight").GetData(detectorSetup.fXParameterWeight);
      setup.GetChild("yParameterWeight").GetData(detectorSetup.fYParameterWeight);
      setup.GetChild("aParameterWeight").GetData(detectorSetup.fAParameterWeight);
      setup.GetChild("bParameterWeight").GetData(detectorSetup.fBParameterWeight);
      setup.GetChild("qOverPParameterWeight").GetData(detectorSetup.fQOverPParameterWeight);
      
      //Populate list of detectors to search for merges.
      utl::Branch mergeListBranch = setup.GetChild("detectorsToMerge");
      std::vector<string> mergeDetectorNames;
      mergeListBranch.GetData(mergeDetectorNames);
      //Temporary container for detector IDs to merge with.      
      std::vector<det::Const::EId> mergeDetectors;
      for (std::vector<string>::iterator mergeIt = mergeDetectorNames.begin(),
	     mergeEnd = mergeDetectorNames.end(); mergeIt != mergeEnd; ++mergeIt) {
	string mergeDetectorName = *mergeIt;
	//MTPCs use just one detection plane between them.
	if (mergeDetectorName == "MTPCs")
	  mergeDetectorName = "MTPCL";
	det::Const::EId mergeDetectorGlobalId = det::Const::GetId(mergeDetectorName);
	if (mergeDetectorGlobalId == det::Const::eUnknown) {
	  ostringstream err;
	  err << "Unknown detector " << mergeDetectorName;
	  ERROR(err);
	  return eFailure;
	}
	mergeDetectors.push_back(mergeDetectorGlobalId);
      }
      detectorSetup.fDetectorsToMerge = mergeDetectors;
  
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
      else
	throw utl::DoesNotComputeException("Fitter type not correctly chosen! "
					   "Check LocalTrackFinderBR.xml file!");
      
      //Insert into vector (ordered by desired detector merging order).
      fDetectorSetups.push_back(detectorSetup);
      //Insert into map in order to call detector setups by global ID later on.
      fDetectorSetupsByGlobalId
	.insert(pair<det::Const::EId,modutils::DetectorSetup>(globalId,detectorSetup));
    }
    
    ostringstream info;
    info << "Global track merger initialization.";

    if ( fVerbosity >= utl::Verbosity::eTerse )
      INFO(info);

    return eSuccess;
  }
  
  /// Process function for GlobalTrackMergerBR
  VModule::EResultFlag GlobalTrackMergerBR::Process(Event& event, const utl::AttributeMap& /*attr*/) {

    //Get necessary event and detector instances.
    evt::RecEvent& recEvent = event.GetRecEvent();
    proc::Reconstruction& procReconstruction = event.GetProcEvent().GetReconstruction();
    const evt::EventHeader& eventHeader = event.GetEventHeader();
    const unsigned int runNumber = eventHeader.GetRunNumber();
    const unsigned int eventNumber = eventHeader.GetId();
    const det::Detector& detector = det::Detector::GetInstance();
    fDetectorCS = detector.GetDetectorCoordinateSystem();
    const det::TPC& detTPC = detector.GetTPC();
    //Get MTPCR for plane retrieval.
    const det::TPCChamber& detMTPCR = detTPC.GetChamber(det::TPCConst::eMTPCR);
    //Initialize Kalman Fiter infrastructure.
    fKalmanFilter.Init(recEvent);
    fKalmanFilter.SetMinimumStepSize(fMinKalmanFilterStepSize);
    fKalmanFilter.SetMaximumStepSize(fMaxKalmanFilterStepSize);
    //Initialize Straight Track fitter.
    fStraightTrackFitter.Init(recEvent,fDetectorCS);

    if ( fVerbosity >= utl::Verbosity::eAnnoying ) {
      ostringstream msg;
      msg << "Merging local tracks into global tracks. Run " << runNumber
	  << ", event " << eventNumber;
      INFO(msg);
    }

    //Loop through designated detectors to perform track merging in.    
    for (std::vector<modutils::DetectorSetup>::iterator setupIt = fDetectorSetups.begin(),
	   detectorEnd = fDetectorSetups.end(); setupIt != detectorEnd; ++setupIt) {

      const det::Const::EId detectorId = (*setupIt).fGlobalId;
      modutils::DetectorSetup setup = *setupIt;

      //Statistics counter.
      unsigned int nMerges = 0;
      
      if (setup.fDetectorsToMerge.size() == 0)
	continue;

      if (detectorId == det::Const::eUnknown) {
	std::ostringstream err;
	err << "Global detector ID= " << detectorId
	    << " has no corresponding detector.";
	ERROR(err.str());
	throw utl::NonExistentComponentException(err.str());
      }
      
      //Get chamber and make sure it is included in the detector configuration.
      const det::TPCConst::EId& chamberId =
	det::TPCConst::GetId(det::Const::GetName(detectorId));
      if (!detTPC.HasChamber(chamberId))
	continue;
      const det::TPCChamber& chamber = detTPC.GetChamber(chamberId);

      //Loop through tracks by z-plane, so we have z-plane information about
      //the parent track during merging.
      const det::TPCChamber::ZPlaneIdRange parentZPlaneIdRange = chamber.GetZPlaneIdRange();

      ostringstream s;
      s << endl << ".===================================Extending " << det::TPCConst::GetName(chamberId) << " tracks.===================================";
      INFO(s);
      
      //Extend tracks to other desired active detectors.
      for (vector<det::Const::EId>::iterator otherDetectorIds = setup.fDetectorsToMerge.begin(),
	     otherDetectorIdsEnd = setup.fDetectorsToMerge.end();
	   otherDetectorIds != otherDetectorIdsEnd; ++otherDetectorIds) {
	
	//Get the other detector setup. 
	modutils::DetectorSetup otherDetectorSetup = fDetectorSetupsByGlobalId[*otherDetectorIds];
	//Don't do anything if this is the detector we're in.	  
	if (setup.fGlobalId == otherDetectorSetup.fGlobalId)
	  continue;
	
	//Get other chamber TPC instance.
	const det::Const::EId& otherDetectorId = otherDetectorSetup.fGlobalId;
	//TPC-specific. Change at some point to include vertex detector/future detectors.
	const det::TPCConst::EId& otherChamberId =
	  det::TPCConst::GetId(det::Const::GetName(otherDetectorSetup.fGlobalId));
	//Make sure this chamber is included in the current detector configuration.
	if (!detTPC.HasChamber(otherChamberId))
	  continue;
	const det::TPCChamber& otherChamber = detTPC.GetChamber(otherChamberId);
	const det::TPCChamber::ZPlaneIdRange zPlaneIdRange = otherChamber.GetZPlaneIdRange();
	
	//Loop through detection planes in original detector.
	for (unsigned int parentZPlaneId = parentZPlaneIdRange.GetBegin(),
	       parentZPlaneEnd = parentZPlaneIdRange.GetEnd();
	     parentZPlaneId < parentZPlaneEnd; ++parentZPlaneId) {
	  //Get this Z-Plane's detection plane grid.
	  proc::DetectionPlaneGrid& parentGrid =
	    procReconstruction.GetDetectionPlaneGrid(parentZPlaneId);
	  //Container for merged track's new indices and z-planes.
	  map<modutils::TrackIndex,modutils::ZPlaneId> newTracksAndZPlanes;
	  //Get active track indices.
	  const proc::TrackIndexVector& localTracks = parentGrid.GetActiveTracks();
	  //If there's no tracks in this range, move on.
	  if (localTracks.size() == 0)
	    continue;
	  
	  //Loop through local tracks in original TPC.
	  for (proc::TrackIndexVector::const_iterator trackIt = localTracks.begin(),
		 trackEnd = localTracks.end(); trackIt != trackEnd; trackIt++) {
	    
	    //Get the track.
	    rec::Track& track = recEvent.Get<rec::Track>(*trackIt);

	    ostringstream s1;
	    s1 << "Working with " << det::TPCConst::GetName(chamberId) << " track " << track.GetIndex() << ".";
	    INFO(s1);
	    
	    //Don't fit GTPC tracks that have not been merged with other tracks --
	    //we'll save this for the BeamlineTPCTrackFitter.
	    if (chamberId == det::TPCConst::eGTPC && track.GetMomentum().GetMag() == 1*utl::TeV)
	      continue;
	    //Enumerator for deciding which fitter type to use when re-fitting.
	    modutils::EFitter refitType = (setup.fFitterType == modutils::eStraightTrackFitter) ?
	      modutils::eStraightTrackFitter : modutils::eKalmanFilter;

	    //Container for holding merge candidates. We'll select the "best" one
	    //from this list -- this will merge one track per chamber.
	    modutils::MergeCandidates mergeCandidates;
	  
	    //Initialize track parameters for extrapolation using track properties.
	    modutils::KalmanTrackParameters trackParameters(track);

	    //Loop through detection planes in this detector.
	    for (unsigned int zPlaneId = zPlaneIdRange.GetBegin(),
		   zPlaneEnd = zPlaneIdRange.GetEnd(); zPlaneId < zPlaneEnd; ++zPlaneId) {

	      //Get this Z-Plane's detection plane grid and z-position.
	      proc::DetectionPlaneGrid& grid = procReconstruction.GetDetectionPlaneGrid(zPlaneId);
	      const double gridZ = grid.GetGridZ();

	      //If this grid has no tracks, skip it.
	      if (grid.GetNumberOfTracks() == 0)
		continue;

	      //Extrapolate track to this grid's z-position! Use
	      //designated fitter to perform the extrapolation.
	      if (setup.fFitterType == modutils::eKalmanFilter) 
		fKalmanFilter.ExtrapolateToGlobalZ(trackParameters,trackParameters,gridZ);
	      if (setup.fFitterType == modutils::eStraightTrackFitter)
		trackParameters.ExtrapolateStraightFitToZ(gridZ);
	      
	      //Get plane object corresponding to z-plane. If TPC Id
	      //is MTPCL and cluster X is negative, the cluster is in
	      //MTPCR.
	      const utl::Plane& zPlane = (otherChamberId == det::TPCConst::eMTPCL &&
					  trackParameters.GetPosition().GetX() < 0) ?
		detMTPCR.GetZPlane(zPlaneId) :
	        otherChamber.GetZPlane(zPlaneId);
	      if (setup.fFitterType == modutils::eKalmanFilter) 
		fKalmanFilter.ExtrapolateToPlane(trackParameters,
						 trackParameters,
						 zPlane);
	      if (setup.fFitterType == modutils::eStraightTrackFitter)
	        trackParameters.ExtrapolateStraightFitToPlane(zPlane);
	      
	      //Make sure extrapolation worked.	      
	      const utl::Point& position = trackParameters.GetPosition();
	      
	      if (position.GetX() != position.GetX() ||
	      	  position.GetY() != position.GetY() )
		continue;
	      
	      //Make sure we haven't left the chamber. MTPCL requires
	      //a check of MTPCR, since they are treated at the same
	      //time.
	      if (otherChamber.GetId() == det::TPCConst::eMTPCL) {
		if (!otherChamber.IsIn(position) &&
		    !detMTPCR.IsIn(position) )
		  continue;
	      }
	      else 
		if (!otherChamber.IsIn(position))
		  continue;
	      

	      //Get track momentum points from extrapolated track position.
	      const proc::TrackIndexVector& trackIndices =
		grid.GetTracksFromPositionAndRange(position.GetX(),
						   position.GetY(),
						   10*fMaxAverageResidual,
						   10*fMaxAverageResidual);

	      //If there's no tracks in this range, move on.
	      if (trackIndices.size() == 0)
		continue;
	      
	      //Loop through these tracks.
	      for (proc::TrackIndexVector::const_iterator pointedTrackIt = trackIndices.begin(),
		     pointedTrackEnd = trackIndices.end();
		   pointedTrackIt != pointedTrackEnd; pointedTrackIt++) {

		ostringstream s1;
		s1 << "Considering " << det::TPCConst::GetName(otherChamberId) << " track " << *pointedTrackIt << ".";
		INFO(s1);
		
		//Prevent re-merging of the same track, if possible.
		if (track.GetIndex() == *pointedTrackIt)
		  continue;

		//Get this track.
		rec::Track& pointedTrack = recEvent.Get<rec::Track>(*pointedTrackIt);

		//Make sure overlapping tracks aren't merged.
		if (pointedTrack.GetFirstPointOnTrack().GetZ() <
		    track.GetFirstPointOnTrack().GetZ() && 
		    pointedTrack.GetLastPointOnTrack().GetZ()  >
		    track.GetLastPointOnTrack().GetZ()    )
		  continue;
		if (pointedTrack.GetFirstPointOnTrack().GetZ() >
		    track.GetFirstPointOnTrack().GetZ() && 
		    pointedTrack.GetLastPointOnTrack().GetZ()  <
		    track.GetLastPointOnTrack().GetZ()    )
		  continue;
		
		//Get its position and momentum.
		const utl::Point& pointedTrackPosition = pointedTrack.GetMomentumPoint();
		const utl::Vector& pointedTrackMomentum = pointedTrack.GetMomentum();
		
		//Extrapolate track to other track's starting z-position.
		//(Grid z-positions can differ from clusters'!)
		if (setup.fFitterType == modutils::eKalmanFilter)
		  fKalmanFilter.ExtrapolateToGlobalZ(trackParameters,
						     trackParameters,
						     pointedTrackPosition.GetZ());
		if (setup.fFitterType == modutils::eStraightTrackFitter)
		  trackParameters.ExtrapolateStraightFitToZ(pointedTrackPosition.GetZ());
		
		const utl::Point& extrapolatedPosition = trackParameters.GetPosition();
		const utl::Vector& extrapolatedMomentum = trackParameters.GetMomentum();
		const utl::Vector& pointedTrackSlopes =
		  pointedTrackMomentum/pointedTrackMomentum.GetZ();

		//Calculate q/|p| for each track.
		const int trackCharge = track.GetCharge();
		const double trackQOverP = trackCharge/extrapolatedMomentum.GetMag();
		const int pointedTrackCharge = pointedTrack.GetCharge();
		const double pointedTrackQOverP = pointedTrackCharge/pointedTrackMomentum.GetMag();
	      
		const double differenceX =
		  extrapolatedPosition.GetX() - pointedTrackPosition.GetX();
		const double differenceY =
		  extrapolatedPosition.GetY() - pointedTrackPosition.GetY();
		const double differenceA = trackParameters.GetA() -  pointedTrackSlopes.GetX();
		const double differenceB = trackParameters.GetB() - pointedTrackSlopes.GetY();
		const double differenceQOverP = (trackQOverP - pointedTrackQOverP);

		modutils::MergeCandidate candidate(pointedTrack.GetIndex(),
						   differenceX, differenceY,
						   differenceA, differenceB,
						   differenceQOverP, otherDetectorSetup);

		ostringstream s2;
		s2 << det::TPCConst::GetName(chamberId) << " track " << track.GetIndex() << " merge metric: " << candidate.GetMetric() << "." << endl;
		s2 << "Parameters 1: (" << extrapolatedPosition.GetX() << "," << extrapolatedPosition.GetY() << "," << trackParameters.GetA() << "," << trackParameters.GetB() << ")" << endl
		   << "Parameters 2: (" << pointedTrackPosition.GetX() << "," << pointedTrackPosition.GetY() << "," << pointedTrackSlopes.GetX() << "," << pointedTrackSlopes.GetY() << ")" << endl;
		INFO(s2);
		
		candidate.SetZPlaneId(zPlaneId);
		mergeCandidates.insert(candidate);
	      } //End tracks loop.
	    } //End zPlaneId loop.

	    //See if we have any merge candidates. If we do, add the one with
	    //the best metric to the list of tracks to merge.
	    for (modutils::MergeCandidates::iterator candidateIt = mergeCandidates.begin(),
		   candidateEnd = mergeCandidates.end(); candidateIt != candidateEnd;
		 candidateIt++) {
	      const modutils::MergeCandidate& bestCandidate = *candidateIt;
	      rec::Track& bestTrack = recEvent.Get<rec::Track>(bestCandidate.GetRecTrackIndex());
	      if (bestCandidate.GetMetric() < fMaxMergeMetric) {
		//Sorted container of z-planes for tracks we've merged in order
		//to find the most-upstream one.
		const modutils::ZPlaneId& newZPlane =
		  (parentZPlaneId < bestCandidate.GetZPlaneId()) ?
		  parentZPlaneId: bestCandidate.GetZPlaneId();
		//If this track has magnetic field fit information,
		//use the Kalman Filter during refitting.
		if (bestCandidate.GetDetectorSetup().fFitterType == modutils::eKalmanFilter)
		  refitType = modutils::eKalmanFilter;
		//Persistent container for holding parent track index
		//while we erase the recEvent objects.
		const Index<rec::Track>& parentTrackIndex = track.GetIndex();
		//Determine new detector ID.
		const det::Const::EId& newDetectorId = (track.GetFirstPointOnTrack().GetZ() < 
							bestTrack.GetFirstPointOnTrack().GetZ()) ?
		  detectorId : otherDetectorId;
		//Try to merge the tracks.
		bool goodMerge = MergeTracks(track,newDetectorId,bestTrack,refitType,event,detTPC);
		//Only perform the administrative merging routines if the refit was good!
		if (goodMerge) {
		  //Remove merge track from the procEvent list of tracks in the detector.
		  procReconstruction.RemoveTrackByDetector(bestCandidate.GetRecTrackIndex(),
							   otherDetectorId);
		  proc::DetectionPlaneGrid& otherGrid =
		    procReconstruction.GetDetectionPlaneGrid(bestCandidate.GetZPlaneId());
		  otherGrid.AddTrackToDisabledList(bestCandidate.GetRecTrackIndex());
		  otherGrid.DisableTracks();
		  
		  //Remove this track from the procEvent list of tracks in the detector.
		  procReconstruction.RemoveTrackByDetector(parentTrackIndex,detectorId);
		  //Don't remove the track from the grid yet -- that
		  //would destroy the loop we're in!
		  parentGrid.AddTrackToDisabledList(parentTrackIndex);
		  const Index<rec::Track>& mergedIndex = recEvent.Back<rec::Track>().GetIndex();
		  //Store this track's index and z-plane for adding to
		  //the correct detection plane grid -- once we're
		  //done looping through this chamber's tracks!
		  newTracksAndZPlanes
		    .insert(pair<modutils::TrackIndex,modutils::ZPlaneId>(mergedIndex,newZPlane));
		  nMerges++;

		  break;
		} //End if (goodMerge).
	      } //End if (metic < fMaxMetric.)
	    } //End merge candidate loop.
	  } //End local tracks loop.
	  //Now disable merged tracks
	  parentGrid.DisableTracks();
	  //Add newly-merged tracks to grids.
	  for (auto iterator = newTracksAndZPlanes.begin(), itEnd = newTracksAndZPlanes.end();
	       iterator != itEnd; ++iterator) {
	    //Get most-upstream z-plane and add track to it.
	    proc::DetectionPlaneGrid& mostUpstreamGrid =
	      procReconstruction.GetDetectionPlaneGrid(iterator->second);
	    const rec::Track& mergedTrack = recEvent.Get<rec::Track>(iterator->first);
	    mostUpstreamGrid.AddTrackToGrid(mergedTrack.GetIndex(),mergedTrack.GetMomentumPoint());
	  }
	} //End parent z-plane loop.
      } //End other detectors loop.
      if (fVerbosity > utl::Verbosity::eTerse && nMerges > 0) {
	const string detectorName = (setup.fGlobalId == det::Const::eMTPCL || 
				     setup.fGlobalId == det::Const::eMTPCR  ) ? 
	  "MTPCs" : det::Const::GetName(setup.fGlobalId);
	ostringstream info;
	info << "Merged " << nMerges << " " << detectorName << " tracks.";
	INFO(info);
      }

    } //End detector setups loop.

    return eSuccess;
  }

  /// Finish function for GlobalTrackMergerBR
  VModule::EResultFlag GlobalTrackMergerBR::Finish() {

    ostringstream info;
    info << "Global tracking Finished." << "\n";
    INFO(info);

    return eSuccess;
  }

  //Function for merging a track with a list of compatible tracks. Returns index of the newly-merged track.
  bool GlobalTrackMergerBR::MergeTracks(rec::Track& parentTrack,
					const det::Const::EId& newDetectorId,
					rec::Track& trackToMerge,
					modutils::EFitter refitType,
					evt::Event& event,
					const det::TPC& detTPC) {

    //Create container for all clusters on new track. Reseve size.
    proc::ClusterVector clusters;
    double nNewClusters = parentTrack.GetNumberOfClusters() + trackToMerge.GetNumberOfClusters();
    clusters.reserve(nNewClusters);

    //Copy clusters from parent track.
    for (rec::ClusterIndexIterator clusterIt = parentTrack.ClustersBegin(),
	   clusterEnd = parentTrack.ClustersEnd(); clusterIt != clusterEnd; clusterIt++) {
      clusters.push_back(*clusterIt);
    }
    //Copy clusters from theeach merge track.
    for (rec::ClusterIndexIterator clusterIt = trackToMerge.ClustersBegin(),
	   clusterEnd = trackToMerge.ClustersEnd(); clusterIt != clusterEnd; clusterIt++) {
      clusters.push_back(*clusterIt);
    }
    

    //Calculate first and last clusters on new track.
    proc::ClusterIndex firstClusterIndex = *(clusters.begin());
    proc::ClusterIndex lastClusterIndex = firstClusterIndex;
    double firstClusterZ = (event.GetRecEvent().Get<rec::Cluster>(firstClusterIndex)).GetPosition().GetZ();
    double lastClusterZ = firstClusterZ;
    for (proc::ClusterVector::const_iterator clusterIt = clusters.begin(),
	   clusterEnd = clusters.end(); clusterIt != clusterEnd; clusterIt++) {
      const double clusterZ = (event.GetRecEvent().Get<rec::Cluster>(*clusterIt)).GetPosition().GetZ();
      if (clusterZ < firstClusterZ) {
	firstClusterIndex = *clusterIt;
	firstClusterZ = clusterZ;
      }
      if (clusterZ > lastClusterZ) {
	lastClusterIndex = *clusterIt;
	lastClusterZ = clusterZ;
      }
    }
    
    //Get recEvent.
    evt::RecEvent& recEvent = event.GetRecEvent();
    //Get first cluster and last on track.
    const rec::Cluster& firstCluster = recEvent.Get<rec::Cluster>(firstClusterIndex);
    const rec::Cluster& lastCluster = recEvent.Get<rec::Cluster>(lastClusterIndex);
    //Get first z-plane.
    const det::TPCConst::EId& chamberId =
      det::TPCConst::GetId(det::Const::GetName(newDetectorId));
    const det::TPCChamber& chamber = detTPC.GetChamber(chamberId);
    const det::TPCPadrow& padrow =
      chamber.GetSector(firstCluster.GetSectorNumber()).GetPadrow(firstCluster.GetPadrowNumber());
    const utl::Plane& zPlane = chamber.GetZPlane(padrow.GetZPlaneId());

    //Re-fit clusters using designated fitter.
    bool fitFlag = false;
    modutils::KalmanTrackParameters trackParameters;
    //Do a Kalman Filter fit if we should.
    if (refitType == modutils::eKalmanFilter) {
      fitFlag = fKalmanFilter.KalmanFit(clusters,trackParameters);
      fKalmanFilter.ExtrapolateToPlane(trackParameters,trackParameters,zPlane);
    }
    //Otherwise, use the Straight Track Fitter.
    else {
      modutils::StraightTrackFitResult fitResult;
      fitFlag = fStraightTrackFitter.FitStraightLine(clusters, fitResult);
      trackParameters.SetStraightFitParameters(fitResult,zPlane);
    }

    //Residual cut for nonsense refits.
    if (trackParameters.GetAverageXResidual() > fMaxAverageResidual ||
	trackParameters.GetAverageYResidual() > fMaxAverageResidual  )
      fitFlag = false;
    if (!fitFlag) {
      return fitFlag;
    }
    
    //Calculate/get track information.
    const utl::Point momentumPoint = trackParameters.GetPosition();
    const double trackA = trackParameters.GetA();
    const double trackB = trackParameters.GetB();
    const double largeMomentum = 100000*utl::GeV;
    const double pMag = (trackParameters.GetQOverP() == 0) ?
      largeMomentum : fabs(1/trackParameters.GetQOverP());
    const double pZ = pMag*(1/sqrt(1 + trackA*trackA + trackB*trackB));
    const double pX = trackA*pZ;
    const double pY = trackB*pZ;
    const int charge = round(trackParameters.GetQOverP()*pMag);
    const utl::Vector momentum(pX,pY,pZ);

    //Create track.
    rec::Track track;
    track.SetStatus(0);
    
    //Store track information.
    track.SetMomentumPoint(momentumPoint);
    track.SetMomentum(momentum);
    track.SetCharge(charge);
    track.SetChi2(trackParameters.GetChi2());
    track.SetNdf(trackParameters.GetNdf());
    track.SetCovarianceMatrix(trackParameters.GetCovarianceMatrix(utl::TrackParameters::eNA61));
    track.SetFirstPointOnTrack(momentumPoint);
    track.SetLastPointOnTrack(lastCluster.GetPosition());

    //Count number of fit clusters in each detector. 
    unsigned int nClustersAll = 0;
    unsigned int nFitClustersAll = 0;
    for (unsigned int iterator = 1, n = rec::TrackConst::ENumberOf::eNTPCs;
	 iterator < n; ++iterator) {
      //Cast the int as a TrackConst...
      rec::TrackConst::ETPC detector = (rec::TrackConst::ETPC)iterator;
      unsigned int nClusters = parentTrack.GetNumberOfClusters(detector);
      unsigned int nFitClusters = parentTrack.GetNumberOfFitClusters(detector);
      nClusters += trackToMerge.GetNumberOfClusters(detector);
      nFitClusters += trackToMerge.GetNumberOfFitClusters(detector);
      track.SetNumberOfClusters(detector,nClusters);
      track.SetNumberOfFitClusters(detector,nFitClusters);
  
      nClustersAll += nClusters; 
      nFitClustersAll += nFitClusters; 
    }
    track.SetNumberOfClusters(rec::TrackConst::eAll,nClustersAll);
    track.SetNumberOfFitClusters(rec::TrackConst::eAll, nFitClustersAll);
    
    //Add track to recEvent.
    evt::rec::Track& recTrack = recEvent.PushBack<evt::rec::Track>(track);

    //Remove old tracks from recEvent.
    parentTrack.Detach(event);
    event.GetRecEvent().Erase<rec::Track>(parentTrack.GetIndex());
    event.GetRecEvent().Erase<rec::Track>(trackToMerge.GetIndex());
    
    //Loop through clusters and link to parent track.
    for (proc::ClusterVector::const_iterator clusterIt = clusters.begin(),
	   clusterEnd = clusters.end(); clusterIt != clusterEnd; clusterIt++) {
      evt::rec::Cluster& cluster = recEvent.Get<rec::Cluster>(*clusterIt);
      evt::IndexedObjectLinker::LinkDaughterToParent(cluster, recTrack);
    }

    //Add this track to procEvent.
    proc::Reconstruction& procReconstruction = event.GetProcEvent().GetReconstruction();
    procReconstruction.AddTrackByDetector(recTrack.GetIndex(),newDetectorId);

    //Return success.
    return fitFlag;
  }
  
  //Function for removing tracks from ProcEvent after they've been merged.  
  void GlobalTrackMergerBR::RemoveTracksFromProcEvent(proc::Reconstruction& procReconstruction,
						      const rec::Track& track,
						      const det::Const::EId& detectorId,
						      const unsigned int& zPlaneId) {

    //Remove this track from the procEvent list of tracks in the detector.
    procReconstruction.RemoveTrackByDetector(track.GetIndex(),detectorId);
    //Remove this track from the corresponding grid structure in ProcEvent.
    //Get this Z-Plane's detection plane grid.
    proc::DetectionPlaneGrid& grid = procReconstruction.GetDetectionPlaneGrid(zPlaneId);
    grid.AddTrackToDisabledList(track.GetIndex());
    grid.DisableTracks();
  }  


}

