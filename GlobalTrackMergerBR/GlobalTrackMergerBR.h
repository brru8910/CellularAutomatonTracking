/**
  \file
  Declaration of GlobalTrackMergerBR

  \author B. Rumberger
  \version $Id$
  \date 5 March 2018
*/

#ifndef _GlobalTrackMergerBR_GlobalTrackMergerBR_h_
#define _GlobalTrackMergerBR_GlobalTrackMergerBR_h_

#include <fwk/VModule.h>

#include <modutils/LocalTrackingTypes.h>
#include <modutils/KalmanFilterWB.h>
#include <modutils/StraightTrackFitter.h>

#include <det/TPC.h>

#include <unordered_map>

namespace GlobalTrackMergerBR {

  /**
    \class GlobalTrackMergerBR
    \author B. Rumberger
    \brief Global Tracking module. Merges tracks between detectors using Kalman Filter extrapolation.

    \ingroup ReconstructionModules
  */

  class GlobalTrackMergerBR : public fwk::VModule {
    
    typedef std::pair<evt::Index<evt::rec::Track>,evt::Index<evt::rec::Track>> TrackIndexPair;

  public:
    fwk::VModule::EResultFlag Init();
    fwk::VModule::EResultFlag Process(evt::Event& event, const utl::AttributeMap& attr);
    fwk::VModule::EResultFlag Finish();
    
    //Function for merging a track with a list of compatible tracks.
    //Returns bool indicating whether or not fit is good.
    bool MergeTracks(evt::rec::Track& parentTrack,
		     const det::Const::EId& parentDetectorId,
		     evt::rec::Track& trackToMerge,
		     modutils::EFitter refitType,
		     evt::Event& event,
		     const det::TPC& detTPC);
    //Function for removing tracks from ProcEvent after they've been merged.
    void RemoveTracksFromProcEvent(evt::proc::Reconstruction& procReconstruction,
				   const evt::rec::Track& track,
				   const det::Const::EId& detectorId,
				   const unsigned int& zPlaneId);

    //Keep administrative functions defined in here.
    void MarkTracksForMerging(const evt::Index<evt::rec::Track> trackIndex1,
			      const evt::Index<evt::rec::Track> trackIndex2)
    {
      TrackIndexPair pair(trackIndex1,trackIndex2);
      fTrackMergePairs.push_back(pair);
    }

  private:
    //Kalman filter max and min step sizes.
    double fMinKalmanFilterStepSize;
    double fMaxKalmanFilterStepSize;
    //Max cap of merging metric across detectors.
    double fMaxMergeMetric;
    //Max average residual upon merging tracks.
    double fMaxAverageResidual;
    //Detector setup vector for initialization.
    std::vector<modutils::DetectorSetup> fDetectorSetups;
    //Detector setup map pulling setups by detector ID.
    std::map<det::Const::EId,modutils::DetectorSetup> fDetectorSetupsByGlobalId;
    //Kalman fitter for local track construction.
    modutils::KalmanFilter fKalmanFilter;
    //Straight track fitter for local track construction in zero magnetic field regions.
    modutils::StraightTrackFitter fStraightTrackFitter;
    //Detector's coordinate system.
    utl::CoordinateSystemPtr fDetectorCS;
    //List of tracks marked for merging.
    std::vector<TrackIndexPair> fTrackMergePairs;

    REGISTER_MODULE("GlobalTrackMergerBR",GlobalTrackMergerBR, "$Id$");
  };

  /**
    \class GlobalTrackMergerBR::MergeCandidate
    \author B. Rumberger
    \brief Class for holding and sorting merge candidate tracks. Sorts according to a metric 
    defined in the less operator.

    \ingroup ReconstructionModules
  */
  
  class MergeCandidate {
    
  public:
    //Default constructor using 
  MergeCandidate(const evt::proc::TrackIndex& index, const double& compatibilityMetric, 
		 const double& positionDifferenceX, const double& positionDifferenceY,
		 const double& momentumDifferenceX, const double& momentumDifferenceY,
		 const unsigned int& zPlaneId) :
    fTrackIndex(index),
      fCompatibilityMetric(compatibilityMetric),
      fPositionDifferenceX(positionDifferenceX),
      fPositionDifferenceY(positionDifferenceY),
      fMomentumDifferenceX(momentumDifferenceX),
      fMomentumDifferenceY(momentumDifferenceY),
      fZPlaneId(zPlaneId)
      {}
    
    ~MergeCandidate() {}
    
    const evt::proc::TrackIndex& GetTrackIndex() const
    { return fTrackIndex; }
    const double& GetCompatibilityMetric() const
    { return fCompatibilityMetric; }
    const double& GetPositionDifferenceX() const
    { return fPositionDifferenceX; }
    const double& GetPositionDifferenceY() const
    { return fPositionDifferenceY; }
    const double& GetMomentumDifferenceX() const
    { return fMomentumDifferenceX; }
    const double& GetMomentumDifferenceY() const
    { return fMomentumDifferenceY; }
    const unsigned int& GetZPlaneId() const
    { return fZPlaneId; }
    
    bool operator<(const MergeCandidate& rhs) const
    { return (fCompatibilityMetric < rhs.fCompatibilityMetric); }
    
  private:
    evt::proc::TrackIndex fTrackIndex;
    double fCompatibilityMetric;
    double fPositionDifferenceX;
    double fPositionDifferenceY;
    double fMomentumDifferenceX;
    double fMomentumDifferenceY;
    unsigned int fZPlaneId;
    
  };

  typedef std::set<MergeCandidate> MergeCandidates;
}


#endif
