/**
   \file
   Declaration of TrackSegmentExtender
   
   \author B. Rumberger
   \version $Id: TrackSegmentExtender.h 2019-10-01 brumberg $
   \date 01 Oct 2019
*/

#ifndef _TrackSegmentExtender_TrackSegmentExtender_h_
#define _TrackSegmentExtender_TrackSegmentExtender_h_

#include <modutils/LocalTrackingTypes.h>
#include <modutils/KalmanFilterWB.h>
#include <modutils/StraightTrackFitter.h>
#include <modutils/ParabolicTrackFitter.h>

#include <det/TPCChamber.h>
#include <evt/proc/DetectionPlaneGrid.h>

#include <evt/Event.h>

#include "TFile.h"

namespace modutils {
  
  /**
     \class  TrackSegmentExtender
     \author B. Rumberger
     \brief  Merges compatible track segment within TPCs. Extends tracks to pick up compatible 
     lone clusters not incorporated into track fit.
  */
  
  class TrackSegmentExtender {
    
  public:
    //Initialize with detector setup and non-changing event structures.
    void Initialize(const modutils::DetectorSetup& detectorSetup,
		    const modutils::ClusterPositionsMap& clusterPositionsByIndex,
		    const modutils::ClusterZPlanesMap& clusterZPlanesByIndex,
		    const det::TPCChamber& detMTPCR)
    {
      fDetectorSetup = detectorSetup;
      fClusterPositionsByIndex = clusterPositionsByIndex;
      fClusterZPlanesByIndex = clusterZPlanesByIndex;
      fMTPCR = &detMTPCR;
    };
    
    //Merging function. Searches for and merges compatible local track segments.
    bool MergeTrackSegments(evt::proc::Reconstruction& procReconstruction,
			    modutils::LocalTrackSegments& localTrackSegments,
			    const det::TPCChamber& chamber);
    //Cluster pickup function. Searches for and adds compatible clusters to local tracks.
    bool PerformClusterPickup(evt::proc::Reconstruction& procReconstruction,
			      modutils::LocalTrackSegments& localTrackSegments,
			      modutils::ClusterList& clustersOnTracks,
			      const det::TPCChamber& chamber);
    //Copy track clusters/cluster info. Utility for track merging.
    void CopyClusters(modutils::LocalTrackSegment& trackToKeep,
		      const modutils::LocalTrackSegment& trackToRemove);
    
    //Set fitter targets.
    void SetKalmanFilter(const modutils::KalmanFilter& kalmanFilter)
    { fKalmanFilter = kalmanFilter; };
    void SetStraightTrackFitter(const modutils::StraightTrackFitter& straightTrackFitter)
    { fStraightTrackFitter = straightTrackFitter; };
    void SetParabolicTrackFitter(const modutils::ParabolicTrackFitter& parabolicTrackFitter)
    { fParabolicTrackFitter = parabolicTrackFitter; };
    
  private:
    //Detector interface.
    modutils::DetectorSetup fDetectorSetup;

    //Fitter interfaces.
    //Kalman fitter for full local track construction and fitting.
    modutils::KalmanFilter fKalmanFilter;
    //Straight track fitter for local track construction in zero magnetic field regions.
    modutils::StraightTrackFitter fStraightTrackFitter;
    //Parabolic track fitter for fast local track construction in non-zero magnetic field regions.
    modutils::ParabolicTrackFitter fParabolicTrackFitter;

    //Cluster positions indexed for fast lookup.
    modutils::ClusterPositionsMap fClusterPositionsByIndex;
    //Cluster z-planes indexed for fast lookup.
    modutils::ClusterZPlanesMap fClusterZPlanesByIndex;

    //Pointer to TPC.
    const det::TPCChamber* fMTPCR;
  };

}

#endif
