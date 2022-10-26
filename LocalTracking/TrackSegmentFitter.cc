#include "TrackSegmentFitter.h"

#include <utl/Point.h>

using namespace std;
using namespace evt;
using namespace evt::proc;


namespace modutils {
  
  bool TrackSegmentFitter::FitTrackCandidates(evt::proc::Reconstruction& procReconstruction,
					      vector<modutils::TrackCandidate>& trackCandidates,
					      const modutils::ClusterPositionsMap& clusterPositionsByIndex,
					      const modutils::ClusterZPlanesMap& clusterZPlanesByIndex,
					      const det::TPCChamber& chamber) {
    
    //Loop through all found potential track segments in chamber. Assumes track candidates have more than
    //three clusters -- check minimum number of clusters in parameters if this is giving an error!
    for (std::vector<modutils::TrackCandidate>::iterator trackCandidatesIt = trackCandidates.begin(),
	   trackCandidatesEnd = trackCandidates.end();
	 trackCandidatesIt != trackCandidatesEnd; ++trackCandidatesIt) {
      modutils::TrackCandidate& candidate = *trackCandidatesIt;
      
      //Make sure there are enough clusters on the candidate.
      const unsigned int nClustersOnTracklet = candidate.GetTrackClusters().size();
      if (nClustersOnTracklet < fDetectorSetup.fMinimumClustersForTrackSeed)
	continue;

      //Order clusters by z-position.
      candidate.SortClusters();

      //Record first z-plane.
      const modutils::ZPlaneId& zPlaneId =
	clusterZPlanesByIndex.at(*candidate.GetTrackClusters().begin());
      //Get plane object corresponding to z-plane. If TPC Id is MTPCL and
      //cluster X is negative, the cluster is in MTPCR.
      const utl::Plane& zPlane =
	(chamber.GetId() == det::TPCConst::eMTPCL &&
	 clusterPositionsByIndex.at(*candidate.GetTrackClusters().begin()).GetX() < 0) ?
	fMTPCR->GetZPlane(zPlaneId) :
	chamber.GetZPlane(zPlaneId);
      
      //Fit segments using Kalman Filter, with new track parameters.
      modutils::KalmanTrackParameters trackParameters;

      bool goodFit = false;
      if (fDetectorSetup.fFitterType == modutils::eKalmanFilter) {
	goodFit = fKalmanFilter.KalmanFit(candidate.GetTrackClusters(),trackParameters);
      }
      if (fDetectorSetup.fFitterType == modutils::eStraightTrackFitter) {
	modutils::StraightTrackFitResult fitResult;
	goodFit = fStraightTrackFitter.FitStraightLine(candidate.GetTrackClusters(), fitResult);
	trackParameters.SetStraightFitParameters(fitResult,zPlane);
      }
      if (fDetectorSetup.fFitterType == modutils::eParabolicTrackFitter) {
	modutils::ParabolicTrackFitResult fitResult;
	goodFit = fParabolicTrackFitter.FitTrack(candidate.GetTrackClusters(), fitResult);
	trackParameters.SetParabolicFitParameters(fitResult,zPlane);
      }
      if (!goodFit)
	continue;

      //Container for holding clusters compatible with fit.
      modutils::ClusterVector acceptedClusters = (trackParameters.GetTrimmedClusters().size() == 0) ?
        candidate.GetTrackClusters() : trackParameters.GetTrimmedClusters();
      
      //Record filtered first z-plane.
      const modutils::ZPlaneId& filteredZPlaneId = clusterZPlanesByIndex.at(*acceptedClusters.begin());
      //Get plane object corresponding to z-plane. If TPC Id is MTPCL and
      //cluster X is negative, the cluster is in MTPCR.
      const utl::Plane& filteredZPlane =
	(chamber.GetId() == det::TPCConst::eMTPCL &&
	 clusterPositionsByIndex.at(*acceptedClusters.begin()).GetX() < 0) ?
	fMTPCR->GetZPlane(filteredZPlaneId) :
	chamber.GetZPlane(filteredZPlaneId);

      //Find most upstream cluster.
      utl::Point firstClusterPosition = clusterPositionsByIndex.at(acceptedClusters.front());
      for (auto it = acceptedClusters.begin(), itEnd = acceptedClusters.end();
	   it != itEnd; ++it) {
	if (clusterPositionsByIndex.at(*it).GetZ() < firstClusterPosition.GetZ())
	  firstClusterPosition = clusterPositionsByIndex.at(*it);
      }
      
      //Re-fit using accepted clusters.
      if (fDetectorSetup.fFitterType == modutils::eKalmanFilter) {
      	fKalmanFilter.KalmanFit(acceptedClusters, trackParameters);
      }
      if (fDetectorSetup.fFitterType == modutils::eStraightTrackFitter) {
	modutils::StraightTrackFitResult fitResult;
	goodFit = fStraightTrackFitter.FitStraightLine(acceptedClusters, fitResult);
	trackParameters.SetStraightFitParameters(fitResult,filteredZPlane);
      }
      if (fDetectorSetup.fFitterType == modutils::eParabolicTrackFitter) {
	modutils::ParabolicTrackFitResult fitResult;
	goodFit = fParabolicTrackFitter.FitTrack(acceptedClusters, fitResult);
	trackParameters.SetParabolicFitParameters(fitResult,filteredZPlane);
      }

      //Average track residual check.
      const double maxResidual = fDetectorSetup.fMaxAverageResidual;
      const double averageTrackResidual =
	sqrt(trackParameters.GetAverageXResidual()*trackParameters.GetAverageXResidual() +
	     trackParameters.GetAverageYResidual()*trackParameters.GetAverageYResidual());
      if (averageTrackResidual > maxResidual)
	continue;
      
      //Container for holding local track properties.
      modutils::LocalTrackSegment localTrackSegment;
      //Make sure track candidate and local track segment indices match.
      localTrackSegment.SetIndex(candidate.GetIndex());

      for (proc::ClusterVector::const_iterator clusterIt = acceptedClusters.begin(), 
      	     clusterEnd = acceptedClusters.end(); clusterIt != clusterEnd; clusterIt++) {
      	//Get cluster information.
      	const utl::Point& clusterPosition = clusterPositionsByIndex.at(*clusterIt);
      	const modutils::ZPlaneId& zPlane = clusterZPlanesByIndex.at(*clusterIt);
      	//Add cluster to track segment.
      	localTrackSegment.AddCluster(*clusterIt,clusterPosition,zPlane);
      	//Keep track of which clusters are used on tracks.
      	fClustersOnTracks.insert(*clusterIt);
      }

      //Sort clusters we've added by z-position.
      localTrackSegment.SortClusters();
      
      if (fDetectorSetup.fFitterType == modutils::eKalmanFilter) {
	//Make sure track parameters are stored on z-plane.
	fKalmanFilter.ExtrapolateToPlane(trackParameters,
					 trackParameters,
					 filteredZPlane);
	localTrackSegment.SetTrackParameters(trackParameters);
      }
      if (fDetectorSetup.fFitterType == modutils::eStraightTrackFitter) {
	modutils::StraightTrackFitResult fitResult;
	fStraightTrackFitter.FitStraightLine(localTrackSegment.GetClusters(), fitResult);
	trackParameters.SetStraightFitParameters(fitResult,filteredZPlane);
	localTrackSegment.SetTrackParameters(trackParameters);
      }
      if (fDetectorSetup.fFitterType == modutils::eParabolicTrackFitter) {
	//FIXME: Extrapolate to z-plane!!!
	modutils::ParabolicTrackFitResult fitResult;
	fParabolicTrackFitter.FitTrack(localTrackSegment.GetClusters(), fitResult);
	trackParameters.SetParabolicFitParameters(fitResult,filteredZPlane);
	localTrackSegment.SetTrackParameters(trackParameters);
      }

      //Add track segment to the list.
      fLocalTrackSegments.push_back(localTrackSegment);
      
      //Add track segment to detection plane grid.
      proc::DetectionPlaneGrid& grid =
	procReconstruction.GetDetectionPlaneGrid(localTrackSegment.GetFirstZPlaneId());

      //Extrapolate to grid.
      utl::Point gridPosition;
      if (fDetectorSetup.fFitterType == modutils::eKalmanFilter) {
	fKalmanFilter.ExtrapolateToPlane(trackParameters,
					 trackParameters,
					 filteredZPlane);
	gridPosition = trackParameters.GetPosition();
      }
      if (fDetectorSetup.fFitterType == modutils::eStraightTrackFitter) {
	trackParameters.ExtrapolateStraightFitToPlane(filteredZPlane);
	gridPosition = trackParameters.GetPosition();
      }
      if (fDetectorSetup.fFitterType == modutils::eParabolicTrackFitter) {
	trackParameters.ExtrapolateParabolicFitToPlane(filteredZPlane);
	gridPosition = trackParameters.GetPosition();
      }
      grid.AddTrackSegmentToGrid(localTrackSegment.GetIndex(),gridPosition);

    } //End candidate loop.
    return true;
  }
  
}
