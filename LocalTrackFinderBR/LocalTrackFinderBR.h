/**
  \file
  Declaration of LocalTrackFinderBR

  \author B. Rumberger
  \version $Id$
  \date 31 Oct 2017
*/

#ifndef _LocalTrackFinderBR_LocalTrackFinderBR_h_
#define _LocalTrackFinderBR_LocalTrackFinderBR_h_

#include <fwk/VModule.h>
#include <utl/Branch.h>
#include <det/MagneticFieldTracker.h>
#include <det/TPCConst.h>
#include <det/TPCChamber.h>

#include <modutils/CellularAutomaton.h>
#include <modutils/TrackletConstructor.h>
#include <modutils/TrackSegmentFitter.h>
#include <modutils/TrackSegmentExtender.h>
#include <modutils/LocalTrackingTypes.h>
#include <modutils/KalmanFilterWB.h>
#include <modutils/StraightTrackFitter.h>
#include <modutils/ParabolicTrackFitter.h>

#include <unordered_map>

namespace LocalTrackFinderBR {

  /**
    \class LocalTrackFinderBR
    \author B. Rumberger
    \brief Local track finder using Cellular Automaton algorithm. Can be used in or out of 
    magnetic fields and in any chamber.

    \ingroup ReconstructionModules
  */

  class LocalTrackFinderBR : public fwk::VModule {
    
  public:
    fwk::VModule::EResultFlag Init();
    fwk::VModule::EResultFlag Process(evt::Event& event, const utl::AttributeMap& attr);
    fwk::VModule::EResultFlag Finish();

    //Function for writing tracks to RecEvent.
    void WriteTracks(evt::Event& recEvent,
		     modutils::LocalTrackSegments& localTrackSegments,
		     const modutils::DetectorSetup& setup);
    //Temporary function for writing tracklets to RecEvent (before fitting).
    //This is to test the tracklet construction efficiency.
    void WriteTracklets(evt::Event& event,
			modutils::TrackCandidates& trackCandidates,
			const modutils::DetectorSetup& setup,
			const det::TPCChamber& chamber);
    //Function for initializing grid objects on TPC detection planes.
    void InitializeDetectionPlaneGrids(evt::Event& event,
				       const det::TPC& detTPC);

  private:
    //Factor for defining bin density in DetectionPlaneGrid objects.
    double fBinDensityFactor;
    //Kalman filter max and min step sizes.
    double fMinKalmanFilterStepSize;
    double fMaxKalmanFilterStepSize;
    //Detector setup vector to pass to Cellular Automaton for detector configuration.
    std::vector<modutils::DetectorSetup> fDetectorSetups;
    //Detector's coordinate system (for initializing fitters).
    utl::CoordinateSystemPtr fDetectorCS;

    ///Sub-modules.
    //Instance of CA track finder.
    modutils::CellularAutomaton fCellularAutomaton;
    //Instance of tracklet constructor.
    modutils::TrackletConstructor fTrackletConstructor;
    //Instance of track segment fitter.
    modutils::TrackSegmentFitter fTrackSegmentFitter;
    //Instance of track segment extender (for merging and cluster pickup).
    modutils::TrackSegmentExtender fTrackSegmentExtender;

    ///Fitter interfaces.
    //Kalman fitter for local track construction.
    modutils::KalmanFilter fKalmanFilter;
    //Straight track fitter for local track construction in zero magnetic field regions.
    modutils::StraightTrackFitter fStraightTrackFitter;
    //Parabolic track fitter for local track construction in non-zero magnetic field regions.
    modutils::ParabolicTrackFitter fParabolicTrackFitter;

    ///Cluster position and z-plane information.
    //Unordered map of cluster indices and corresponding positions.
    modutils::ClusterPositionsMap fClusterPositionsByIndex;
    //Unordered map of cluster indices and corresponding errors.
    modutils::ClusterPositionsMap fClusterErrorsByIndex;
    //Unordered map of cluster z-planes by index.
    modutils::ClusterZPlanesMap fClusterZPlanesByIndex;

    REGISTER_MODULE("LocalTrackFinderBR",LocalTrackFinderBR, "$Id$");
  };

}


#endif
