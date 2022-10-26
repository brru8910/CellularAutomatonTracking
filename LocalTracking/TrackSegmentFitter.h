/**
   \file
   Declaration of TrackSegmentFitter

   \author B. Rumberger
   \version $Id: TrackSegmentFitter.h 2019-10-01 brumberg $
   \date 01 Oct 2019
*/

#ifndef _TracklSegmentFitter_TracklSegmentFitter_h_
#define _TracklSegmentFitter_TracklSegmentFitter_h_

#include <modutils/LocalTrackingTypes.h>
#include <modutils/KalmanFilterWB.h>
#include <modutils/StraightTrackFitter.h>
#include <modutils/ParabolicTrackFitter.h>

#include <det/TPCChamber.h>

#include <evt/proc/DetectionPlaneGrid.h>
#include <evt/RecEvent.h>

namespace modutils {
 
  /**
     \class  TrackSegmentFitter
     \author B. Rumberger
     \brief  Fits track candidates found in TPCs. Models available are straight fit, 
     parabolic fit, and Kamlan Filter fit.
  */

  class TrackSegmentFitter {

  public:
    //Initialize with detector setup.
    void Initialize(const modutils::DetectorSetup& detectorSetup,
		    const det::TPCChamber& detMTPCR)
    {
      fDetectorSetup = detectorSetup;
      fLocalTrackSegments.clear();
      fClustersOnTracks.clear();
      fMTPCR = &detMTPCR;
    };
    bool FitTrackCandidates(evt::proc::Reconstruction& procReconstruction,
			    std::vector<modutils::TrackCandidate>& trackCandidates,
			    const modutils::ClusterPositionsMap& clusterPositionsByIndex,
			    const modutils::ClusterZPlanesMap& clusterZPlanesByIndex,
			    const det::TPCChamber& chamber);
    //Set fitter targets.
    void SetKalmanFilter(const modutils::KalmanFilter& kalmanFilter)
    { fKalmanFilter = kalmanFilter; };
    void SetStraightTrackFitter(const modutils::StraightTrackFitter& straightTrackFitter)
    { fStraightTrackFitter = straightTrackFitter; };
    void SetParabolicTrackFitter(const modutils::ParabolicTrackFitter& parabolicTrackFitter)
    { fParabolicTrackFitter = parabolicTrackFitter; };
    
    //Get vector of local tracks.
    modutils::LocalTrackSegments& GetTrackSegments()
      { return fLocalTrackSegments; };
    //Get clusters added to tracks.
    modutils::ClusterList& GetUsedClusters()
      { return fClustersOnTracks; };
    
      
  private:
    //Detector interface.
    modutils::DetectorSetup fDetectorSetup;
    //Kalman fitter for full local track construction and fitting.
    modutils::KalmanFilter fKalmanFilter;
    //Straight track fitter for local track construction in zero magnetic field regions.
    modutils::StraightTrackFitter fStraightTrackFitter;
    //Parabolic track fitter for fast local track construction in non-zero magnetic field regions.
    modutils::ParabolicTrackFitter fParabolicTrackFitter;

    //Storage container for good track segments.
    modutils::LocalTrackSegments fLocalTrackSegments;
    //Storage container for cluster indices that belong to a track.
    modutils::ClusterList fClustersOnTracks;

    //Container for holding MTPCR. Needed since MTPCL & MTPCR are treated with one set of z-planes.
    const det::TPCChamber* fMTPCR;
    
  };

}

#endif
