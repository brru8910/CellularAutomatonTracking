/**
   \file                                                                                                                         
   Declaration of TrackletConstructor                                                                                                  
   
   \author B. Rumberger                                                                                                              
   \version $Id: TrackletConstructor.h 1 2019-9-30 12:00:00Z brumberg $
   \date 30 Sept 2019                                                                                                                 
*/

#ifndef _TrackletConstructor_TrackletConstructor_h_
#define _TrackletConstructor_TrackletConstructor_h_

#include <modutils/CellularAutomaton.h>
#include <modutils/LocalTrackingTypes.h>

#include <evt/proc/DetectionPlaneGrid.h>
#include <utl/Point.h>

namespace modutils {
  

  /**
     \class TrackletConstructor                                                                                                                                                                                     
     \author B. Rumberger                                                                                                                    
     \brief Constructs links, aka tracklets, between clusters based on a pattern recognition model.
    */
    
  class TrackletConstructor {
    
  public:
    //Initialize TrackletConstructor with detector setup.
    void Initialize(const modutils::DetectorSetup& detectorSetup)
    { fDetectorSetup = detectorSetup;
      fTrackletIndex = 1; };
    //Function for forming tracklets only between the best candidates. Forward and backward tracklets
    //are considered, and then the straightest pair of forward/backward tracklets are formed. All
    //other possible connections are ignored. Search area can be adjusted with the first parameter
    //for use in multiple iterations.
    bool ConstructTracklets(evt::proc::Reconstruction& procReconstruction,
			    const modutils::ClusterPositionsMap& clusterPositionsByIndex,
			    modutils::CellularAutomaton& cellularAutomaton,
			    double areaSearchScaleFactor);
    
    
  private:
    modutils::DetectorSetup fDetectorSetup;
    unsigned int fTrackletIndex;
    
  };

}

#endif
