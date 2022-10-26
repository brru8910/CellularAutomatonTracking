/**
   \file                                                                                                                         
   Declaration of CellularAutomaton                                                                                                  
   
   \author B. Rumberger                                                                                                              
   \version $Id: CellularAutomaton.h 1 2017-10-30 12:00:00Z brumberg $                                                               
   \date 30 Oct 2017                                                                                                                 
*/

#ifndef _CellularAutomaton_CellularAutomaton_h_
#define _CellularAutomaton_CellularAutomaton_h_

#include <modutils/LocalTrackingTypes.h>
#include <evt/RecEvent.h>
#include <det/TPCConst.h>
#include <utl/Point.h>

//ROOT clutter to get rid of later.
#include <TFile.h>
#include <TNtuple.h>

namespace modutils {
  

  /**
     \class CellularAutomaton                                                                                                                                                                                     
     \author B. Rumberger                                                                                                                    
     \brief Finds local tracks in subdetectors using an algorithm based on the Cellular Automaton 
     tracking scheme.
    */
    
  class CellularAutomaton {
    
    typedef std::vector<std::map<modutils::ClusterIndex,modutils::Node> >::iterator NodePlaneIterator;
    typedef std::vector<std::map<modutils::TrackletIndex,modutils::Tracklet> >::iterator TrackletPlaneIterator;
    typedef std::map<modutils::ClusterIndex,modutils::Node>::iterator NodeIterator;
    typedef std::map<modutils::ClusterIndex,modutils::Node>::reverse_iterator NodeReverseIterator;
    typedef std::map<modutils::TrackletIndex,modutils::Tracklet>::iterator TrackletIterator;
    typedef std::map<modutils::TrackletIndex,modutils::Tracklet>::reverse_iterator TrackletReverseIterator;
    typedef std::vector<modutils::TrackCandidate>::iterator TrackCandidateIterator;
    typedef std::vector<modutils::ClusterIndex>::const_iterator ClusterIndexIterator;    

  public:
    //Initialize local tracker with .xml parameters. Possibly resize vectors with number of z-planes here!
    bool Initialize(const modutils::DetectorSetup& detectorSetup);
    
    //Add a forward-going tracklet to an existing node.
    void AddForwardTrackletToNode(const modutils::ClusterIndex& index,
				  const modutils::Tracklet& tracklet);
      
    //Add a backward-going tracklet to an existing node.
    void AddBackwardTrackletToNode(const modutils::ClusterIndex& index,
				   const modutils::Tracklet& tracklet);

    //Organize all clusters into lists of clusters connected by tracklets.
    void FindTrackCandidates();

    //Find conected clusters recursively from initial cluster.
    void FindConnectedClusters(modutils::TrackCandidate& trackCandidate,
			       const modutils::Node&);
      
    //Return lists of connected clusters constituting tracks.
    modutils::TrackCandidates& GetTrackCandidates()
      { return fTrackCandidates; }

    //Master function for finding local tracks in a detector.
    bool PlotTracklets(TFile* file);

    //Remove tracklets that do not fit local track criteria in YZ-plane.
    void ApplyYZSlopeRule();

    //Remove tracklets that do not fit local track criteria in XZ-plane.
    void ApplXZSlopeRule();

    //Function for filtering spurious tracklets.
    void ApplyTrackletPointingRule();

    //Function for filtering track candidates by rejecting tracklets considered outliers.
    void FilterTrackCandidates();

    //Function for finding most-upstream clusters on each track segment.
    //Creates list of these for joining segments and fitting. 
    void FindFirstClusterOnTrackSegments();

    //Sandbox function for other cuts to try.
    void ApplyCut(utl::Point targetPosition);
    
    //TODO: Move utility routines to here! (all that follow)
    
    //Function to remove a node (once it is assigned to a track) and all tracklets linked to that node.
    void RemoveNode(modutils::Node);

    //Function for purging tracklets from event.
    void RemoveTracklets(modutils::TrackletIndices trackletIndicesToDelete);

    //Add tracklet for consideration.
    void AddTracklet(const modutils::Tracklet& tracklet);
    
    //Add node to z-plane.
    void AddNode(modutils::Node& node);

    //Search for an existing tracklet by index.
    bool SearchTracklet(const modutils::TrackletIndex& index);

    //Search for an existing node on a z-plane by corresponding cluster index.
    bool SearchNode(const modutils::ClusterIndex& index);

    //Get a node from internal container.
    modutils::Node& GetNode(const ClusterIndex& index);

    //Set a node as used.
    void SetNodeUsed(const ClusterIndex& index, const bool isUsed);

    //Assign a starting node to a track segment.
    std::pair<modutils::ClusterIndex,modutils::TrackSegmentIndex>
      MakeNodeAndTrackSegmentPair(const ClusterIndex& clusterIndex, const TrackSegmentIndex& trackIndex) {
      return std::make_pair(clusterIndex,trackIndex);
    };

    //Search track segment first cluster list. Return true if cluster with index clusterIndex
    //appears in the list.
    bool SearchStartCluster(const ClusterIndex& clusterIndex) {

      std::map<modutils::ClusterIndex,modutils::TrackSegmentIndex>::iterator it =
	fTrackSegmentFirstClusters.find(clusterIndex);

      /* msg << endl << "Cluster " << clusterIndex << " corresponds to track " << it->second; */
      // if (it != fTrackSegmentFirstClusters.end()) INFO(msg);

      return (it != fTrackSegmentFirstClusters.end());
    }

    //Get the track segment associated with a start cluster.
    //Should be used after checking for existence with with SearchStartCluster().
    const modutils::TrackSegmentIndex&
      GetTrackSegmentIndexFromStartCluster(const ClusterIndex& clusterIndex) {
      std::map<modutils::ClusterIndex,modutils::TrackSegmentIndex>::iterator it =
	fTrackSegmentFirstClusters.find(clusterIndex);
      if (it == fTrackSegmentFirstClusters.end()) {
	std::ostringstream err;
	err << "Couldn't find a Start Cluster with index " << clusterIndex
	    << "! Node marked as existing in SearchStartCluster!";
	ERROR(err);
      }
      return it->second;
    }
    


  private:
    det::TPCConst::EId fChamberId;
    det::Const::EId fDetectorId;
    utl::CoordinateSystemPtr fDetectorCS;
    double fXZAngleTolerance;
    double fYZAngleTolerance;
    double fXDisplacementTolerance;
    double fYDisplacementTolerance;
    unsigned int fNZPlanesInChamber;
    unsigned int fNPlanesToSkip;
    unsigned int fNTrackSegments;
    unsigned int fMinimumClustersOnLocalTrack;

    //Vector of track candidates containing clusters to be fit to a track.
    std::vector<modutils::TrackCandidate> fTrackCandidates;
    //Vector of clusters at the upstream end of valid track segments and the track segment index.
    std::map<modutils::ClusterIndex,modutils::TrackSegmentIndex> fTrackSegmentFirstClusters;
    //Vector of used cluster indices for finding contiguous clusters.
    std::vector<bool> fUsedClusters;
    //Map of tracklets indexed by tracklet index.
    std::map<modutils::TrackletIndex,modutils::Tracklet> fTracklets;
    //List of nodes indexed by node index.
    std::map<modutils::ClusterIndex,modutils::Node> fNodes;
    //Store detector setup after initialization.
    modutils::DetectorSetup fDetectorSetup;

    //ROOT crap.
    TFile* file;
    TFile* file2;
    TFile* file3;
    TNtuple* tracklets;
    TNtuple* nodes;
    
  };

}

#endif
