/**
  \file
  Declaration of types necessary for local tracking.

  \author B. Rumberger
  \version $Id: LocalTrackingTypes.h 1 2017-10-30 12:00:00Z brumberg $
  \date 30 Oct 2017
*/

#ifndef _LocalTrackingTypes_LocalTrackingTypes_h_
#define _LocalTrackingTypes_LocalTrackingTypes_h_

#include <evt/Event.h>
#include <evt/IndexedObject.h>
#include <evt/IndexList.h>
#include <modutils/KalmanTrackParameters.h>
#include <utl/Range.h>
#include <utl/TrackLocalParameters.h>

#include <unordered_map>

namespace modutils {

  //Enum for using different fitter types.
  enum EFitter {
    eStraightTrackFitter = 0,
    eKalmanFilter = 1,
    eParabolicTrackFitter = 2
  };
  
  typedef unsigned int ZPlaneId;
  typedef utl::Range<unsigned int> ZPlaneIdRange;

  typedef unsigned int NodeIndex;

  typedef evt::Index<evt::rec::Cluster> ClusterIndex;
  typedef evt::Index<evt::rec::Track> TrackIndex;
  typedef unsigned int TrackletIndex;
  typedef int TrackSegmentIndex;
  typedef std::vector<TrackSegmentIndex> TrackSegmentIndexVector;

  typedef std::set<TrackletIndex> TrackletIndices;
  typedef std::set<ClusterIndex> ClusterList;
  typedef std::vector<ClusterIndex> ClusterVector;

  typedef std::unordered_map<unsigned int, utl::Point> ClusterPositionsMap;
  typedef std::pair<ClusterIndex, utl::Point> ClusterIndexAndPositionPair;

  typedef std::unordered_map<unsigned int, unsigned int> ClusterZPlanesMap;
  typedef std::pair<ClusterIndex, unsigned int> ClusterIndexAndZPlanePair;
  
  //Detector setup structure for passing parameters to the CA class.
  struct DetectorSetup {
    det::Const::EId fGlobalId;
    double fXDisplacementTolerance;
    double fYDisplacementTolerance;
    double fMaxAngleChangeXZ;
    double fMaxAngleChangeYZ;
    double fPadLength;
    unsigned int fNPlanesToSkip;    
    unsigned int fMinimumClustersOnLocalTrack;
    unsigned int fMinimumClustersForTrackSeed;
    modutils::EFitter fFitterType;

    //ZPlane ranges for local tracking.
    ZPlaneIdRange fZPlaneIdRange;
    
    //Track merging parameters. Larger weights give more significance to the associated term
    //when deciding if two track segments are compatible. See MergeCandidate class for details.
    double fXParameterWeight;
    double fYParameterWeight;
    double fAParameterWeight;
    double fBParameterWeight;
    double fQOverPParameterWeight;

    //Max merge metric and average residuals for local tracking algorithms.
    double fMaxMergeMetric;
    double fMaxAverageResidual;
    
    //List of detectors in which merging should be performed. Used in Global Tracking.
    std::vector<det::Const::EId> fDetectorsToMerge;

  };

  //A node is a an object that contains links to other nodes (in the form of tracklets). 
  //Each cluster has a corresponding node.
  class Node {
    
  public:
    
    //Initalize node using corresponding cluster index?
  Node() :
    fIndex(0),
      fZPlaneId(0),
      fTPCId(det::TPCConst::eUnknown),
      fIsUsed(false),
      fFirstOnTrackSegment(false)
	{}

    ~Node() {}
    
    void SetIndex(evt::Index<evt::rec::Cluster> clusterIndex) { fIndex = clusterIndex; }
    const evt::Index<evt::rec::Cluster>& GetIndex() const { return fIndex; }

    void SetZPlaneId(modutils::ZPlaneId zPlaneId) { fZPlaneId = zPlaneId; }
    const modutils::ZPlaneId& GetZPlaneId() const {return fZPlaneId; }

    void SetTPCId(det::TPCConst::EId tpcId) { fTPCId = tpcId; }
    const det::TPCConst::EId& GetTPCId() const {return fTPCId; }

    void AddForwardTrackletIndex(TrackletIndex trackletIndex)
    { fForwardTrackletIndices.insert(trackletIndex); }
    void AddBackwardTrackletIndex(TrackletIndex trackletIndex)
    { fBackwardTrackletIndices.insert(trackletIndex); }

    void RemoveForwardTrackletIndex(TrackletIndex trackletIndex) { 
      if (fForwardTrackletIndices.size() == 0) return;
      TrackletIndices::iterator indexIt = fForwardTrackletIndices.find(trackletIndex);
      if (indexIt == fForwardTrackletIndices.end()) return;
      fForwardTrackletIndices.erase(indexIt); 
    }
    void RemoveBackwardTrackletIndex(TrackletIndex trackletIndex) { 
      if (fBackwardTrackletIndices.size() == 0) return;
      TrackletIndices::iterator indexIt = fBackwardTrackletIndices.find(trackletIndex);
      if (indexIt == fBackwardTrackletIndices.end()) return;
      fBackwardTrackletIndices.erase(trackletIndex); 
    }
    
    const TrackletIndices& GetForwardTrackletIndices() const { return fForwardTrackletIndices; }
    const TrackletIndices& GetBackwardTrackletIndices() const { return fBackwardTrackletIndices; }

    void KeepTracklets(TrackletIndices forwardIndicesToKeep, TrackletIndices backwardIndicesToKeep) {
      fForwardTrackletIndices = forwardIndicesToKeep;
      fBackwardTrackletIndices = backwardIndicesToKeep;
    }

    utl::Point GetPosition() const { return fPosition; }
    void SetPosition(utl::Point position) { fPosition = position; }

    bool operator<(const Node& rhs) const {
      return (fPosition.GetZ() < rhs.fPosition.GetZ());
    }

    bool IsUsed() const { return fIsUsed; }

    void SetUsed(bool isUsed) { fIsUsed = isUsed; }

    void SetFirstNodeOnTrackSegment(bool isFirstOnTrackSegment)
    { fFirstOnTrackSegment = isFirstOnTrackSegment; }

    bool IsFirstNodeOnTrackSegment() { return fFirstOnTrackSegment; }

  private:
    evt::Index<evt::rec::Cluster> fIndex;
    //Global z-plane Id.
    modutils::ZPlaneId fZPlaneId;
    det::TPCConst::EId fTPCId;
    TrackletIndices fForwardTrackletIndices;
    TrackletIndices fBackwardTrackletIndices;
    double fMostCommonSlopeXZ;
    double fMostCommonSlopeYZ;
    utl::Point fPosition;
    bool fIsUsed;
    bool fFirstOnTrackSegment;
  };

  typedef std::set<Node> NodeList;

  //A tracklet is a potential local track segment. Tracklets are initially created between
  //every cluster combination on adjacent
  //padrows. These potential track segments are filtered by the CA algorithm and thrown
  //away if they do not meet local track
  //criteria. The most important tracklet characteristics are XZ/YZ angles, position, and linked clusters.  
  class Tracklet {

  public:
  Tracklet() : 
    fFirstZPlaneId(0), 
      fLastZPlaneId(0),
      fFirstClusterIndex(0),
      fLastClusterIndex(0),
      fDetectorId(det::TPCConst::eUnknown),
      fNumberOfNeighbors(0),
      fPosition(utl::Point(0,0,0)),
      fAngleTangentXZ(0),
      fAngleTangentYZ(0),
      fAngleTangentXY(0)
	{}
    
  Tracklet(const ClusterIndex& firstNodeIndex,
	   const utl::Point& firstNodePosition,
	   const ClusterIndex& lastNodeIndex,
	   const utl::Point& lastNodePosition,
	   TrackletIndex index,
	   const det::TPCConst::EId& detectorId) : 
    fFirstClusterIndex(firstNodeIndex),
      fLastClusterIndex(lastNodeIndex),
      fDetectorId(detectorId),
      fIndex(index) {

	//Calculate average position and store.
	const double x = (firstNodePosition.GetX() + lastNodePosition.GetX())/2;
	const double y = (firstNodePosition.GetY() + lastNodePosition.GetY())/2;
	const double z = (firstNodePosition.GetZ() + lastNodePosition.GetZ())/2;
	fPosition.Set(x,y,z);

	//Now we mess around with the STL container holding the clusters.
	//This is the current bottleneck in the constructor! 
	fConstituentClusters.reserve(2);
	//TODO: MAKE THESE AN ENUM!!
	fConstituentClusters.push_back(firstNodeIndex);
	fConstituentClusters.push_back(lastNodeIndex);

	//Calculate fAngleTangents.
	const double fDisplacementZ = lastNodePosition.GetZ() - firstNodePosition.GetZ();
	fAngleTangentXZ = (lastNodePosition.GetX() - firstNodePosition.GetX())/fDisplacementZ;
	fAngleTangentYZ = (lastNodePosition.GetY() - firstNodePosition.GetY())/fDisplacementZ;
	fAngleTangentXY = fAngleTangentXZ / fAngleTangentYZ;
      }
    
    ~Tracklet() {}
    
    const utl::Point& GetPosition() const { return fPosition;}
    double GetAngleXZ() const {return std::atan(fAngleTangentXZ);}
    double GetAngleYZ() const {return std::atan(fAngleTangentYZ);}
    double GetAngleXY() const {return std::atan(fAngleTangentXY);}
    const utl::Vector GetAngleTangents() const {return utl::Vector(fAngleTangentXZ,fAngleTangentYZ,1);}
    const TrackletIndex& GetIndex() const {return fIndex;} 
    void AddCluster(ClusterIndex index) {
      //Since fConstituentClusters is a std::set and rec::Clusters are compared by z-position,
      //this container will remain sorted.
      fConstituentClusters.push_back(index);
    }
    const std::vector<ClusterIndex>& GetConstituentClusters() const {return fConstituentClusters;}
    const ClusterIndex& GetFirstClusterIndex() const {return fFirstClusterIndex;}
    const ClusterIndex& GetLastClusterIndex() const {return fLastClusterIndex;}
    //Work in progress. Not sure what other fields/methods are required.

  private:
    ZPlaneId fFirstZPlaneId;
    ZPlaneId fLastZPlaneId;
    ClusterIndex fFirstClusterIndex;
    ClusterIndex fLastClusterIndex;
    det::TPCConst::EId fDetectorId;
    unsigned int fNumberOfNeighbors;
    utl::Point fPosition;
    TrackletIndex fIndex;
    double fAngleTangentXZ;
    double fAngleTangentYZ;
    double fAngleTangentXY;
    std::vector<ClusterIndex> fConstituentClusters;

  };

  typedef std::vector<Tracklet> TrackletVector;

  //Class for collecting clusters connected by tracklets and deciding whether or not to add them
  //to a potential track.
  class TrackCandidate {

  public:
    
  TrackCandidate(const Node& node) :
    fNClusters(1),
      fIndex(-1),
      fAnglesSorted(false)
	{ 
	  fTrackClusters.push_back(node.GetIndex());
	  fTrackNodes.insert(node); 
	  fNodeIndicesByZPosition.insert(std::make_pair(node.GetPosition().GetZ(),node.GetIndex()));
	}
    
    ~TrackCandidate() {}
    
    void AddNodeToTrack(const Node& node, const Tracklet& tracklet) {
      fNClusters++;
      fTrackNodes.insert(node);
      fTrackTracklets.push_back(tracklet);
      fTrackClusters.push_back(node.GetIndex());
      fNodeIndicesByZPosition.insert(std::make_pair(node.GetPosition().GetZ(),node.GetIndex()));
      //Update angle containers.
      fXZAngles.push_back(tracklet.GetAngleXZ());
      fYZAngles.push_back(tracklet.GetAngleYZ());
    }

    void CalculateAngleStatistics() {

      //Sort angle containers, if we haven't done so.
      if (fAnglesSorted == false) {
	std::sort(fXZAngles.begin(),fXZAngles.end());
	std::sort(fYZAngles.begin(),fYZAngles.end());
	fAnglesSorted = true;
      }

      unsigned int nAngles = fXZAngles.size();

      //XZ statistics calculation.
      fXZAngleMedian = (nAngles % 2 == 0) ?
      	(fXZAngles[nAngles/2 - 1] + fXZAngles[nAngles/2])/2 : fXZAngles[floor(nAngles/2)];

      fXZAngleLowerQuartile = (nAngles % 4 == 0) ?
      	(fXZAngles[nAngles/4 - 1] + fXZAngles[nAngles/4])/2 : fXZAngles[floor(nAngles/4)];

      fXZAngleUpperQuartile = (3*nAngles % 4 == 0) ?
      	(fXZAngles[3*nAngles/4 - 1] + fXZAngles[3*nAngles/4])/2 : fXZAngles[floor(3*nAngles/4)];

      //YZ statistics calculation.
      fYZAngleMedian = (nAngles % 2 == 0) ?
      	(fYZAngles[nAngles/2 - 1] + fYZAngles[nAngles/2])/2 : fYZAngles[floor(nAngles/2)];

      fYZAngleLowerQuartile = (nAngles % 4 == 0) ?
      	(fYZAngles[nAngles/4 - 1] + fYZAngles[nAngles/4])/2 : fYZAngles[floor((nAngles+1)/4)];

      fYZAngleUpperQuartile = (3*nAngles % 4 == 0) ?
      	(fYZAngles[3*nAngles/4 - 1] + fYZAngles[3*nAngles/4])/2 : fYZAngles[floor(3*nAngles/4)];
      
    }
    //Stats getters.
    double GetXZAngleMedian() { return fXZAngleMedian;}
    double GetXZAngleLowerQuartile() { return fXZAngleLowerQuartile;}
    double GetXZAngleUpperQuartile() { return fXZAngleUpperQuartile;}
    double GetXZAngleInterQuartileRange() { return fabs(fXZAngleUpperQuartile - fXZAngleLowerQuartile);}
    double GetYZAngleMedian() { return fYZAngleMedian;}
    double GetYZAngleLowerQuartile() { return fYZAngleLowerQuartile;}
    double GetYZAngleUpperQuartile() { return fYZAngleUpperQuartile;}
    double GetYZAngleInterQuartileRange() { return fabs(fYZAngleUpperQuartile - fYZAngleLowerQuartile);}

    /* double GetXZAngleRange() { */
    /*   if (fAnglesSorted == false) { */
    /* 	std::sort(fXZAngles.begin(),fXZAngles.end()); */
    /* 	std::sort(fYZAngles.begin(),fYZAngles.end()); */
    /* 	fAnglesSorted = true; */
    /*   } */

    /*   double maxXZ =  */
    /* } */

    //Index handling.
    void SetIndex(const TrackSegmentIndex& index) { fIndex = index; }
    const TrackSegmentIndex& GetIndex() const { return fIndex; }

    const ClusterVector& GetTrackClusters() const { return fTrackClusters; }
    const TrackletVector& GetTracklets() const { return fTrackTracklets; }
    const NodeList& GetTrackNodes() const { return fTrackNodes; }
    const unsigned int& GetNClusters() const { return fNClusters; }

    void SortClusters() {
      fTrackClusters.clear();
      for (std::map<double,int>::iterator indexIt = fNodeIndicesByZPosition.begin(),
	     indexEnd = fNodeIndicesByZPosition.end(); indexIt != indexEnd; indexIt++) {
	fTrackClusters.push_back(indexIt->second);
      }
    }
    
  private:
    unsigned int fNClusters;
    TrackSegmentIndex fIndex;
    ClusterVector fTrackClusters;
    TrackletVector fTrackTracklets;
    NodeList fTrackNodes;
    //Map of node indices and corresponding z-positions.
    std::map<double,int> fNodeIndicesByZPosition;
    std::vector<double> fXZAngles;
    std::vector<double> fYZAngles;
    bool fAnglesSorted;
    //Statistics containers.
    double fXZAngleMedian;
    double fXZAngleLowerQuartile;
    double fXZAngleUpperQuartile;
    double fYZAngleMedian;
    double fYZAngleLowerQuartile;
    double fYZAngleUpperQuartile;
  };

  
  typedef std::vector<modutils::TrackCandidate> TrackCandidates;

  
  //Class for storing and sorting potential tracklets in minimal tracklet CA mode.
  class TrackletCandidate {
    
  public:
  TrackletCandidate(const double xzAngleDifference,
		    const double yzAngleDifference,
		    const ClusterIndex& clusterIndexPreviousLayer,
		    const ClusterIndex& centralClusterIndex,
		    const ClusterIndex& clusterIndexNextLayer,
		    const ZPlaneId& zPlaneId) :
    fXZAngleDifference(xzAngleDifference),
      fYZAngleDifference(yzAngleDifference),
      fClusterIndexPreviousLayer(clusterIndexPreviousLayer),
      fCentralClusterIndex(centralClusterIndex),
      fClusterIndexNextLayer(clusterIndexNextLayer),
      fZPlaneId(zPlaneId)
    //Sort on equal combinaion of XZ and YZ slopes.
      { 
	fConstructionMetric = sqrt(fYZAngleDifference*fYZAngleDifference +
				   fXZAngleDifference*fXZAngleDifference); }
    
    const ClusterIndex& GetClusterIndexNextLayer() const {return fClusterIndexNextLayer;}
    const ClusterIndex& GetCentralClusterIndex() const {return fCentralClusterIndex;}
    const ClusterIndex& GetClusterIndexPreviousLayer() const {return fClusterIndexPreviousLayer;}
    const double& GetConstructionMetric() const {return fConstructionMetric;}
    unsigned int GetZPlaneId() const {return fZPlaneId;}
    
    bool operator<(const TrackletCandidate& rhs) const {
      return (fConstructionMetric < rhs.fConstructionMetric);
    }

  private:
    double fXZAngleDifference;
    double fYZAngleDifference;
    ClusterIndex fClusterIndexPreviousLayer;
    ClusterIndex fCentralClusterIndex;
    ClusterIndex fClusterIndexNextLayer;
    ZPlaneId fZPlaneId;
    double fConstructionMetric;
  };

  //Class for storing track objects when fitting and extending TrackCandidates.
  //Routines for storing track parameters at start/end of tracks, adding clusters,
  //and copying/removing clusters.
  class LocalTrackSegment {
    
  public:
  LocalTrackSegment() :
      fIndex(-1),
      fIsSorted(false)
    {}

      //Add cluster to local track segment. Must include the cluster's position and
      //detection plane for sorting purposes.
      void AddCluster(const ClusterIndex& index,
		      const utl::Point& clusterPosition,
		      const modutils::ZPlaneId& zPlaneId) {
	fTrackClusters.push_back(index);
	fClusterIndicesByZPosition.insert(std::make_pair(clusterPosition.GetZ(),index));
	fClusterZPlanesByZPosition.insert(std::make_pair(clusterPosition.GetZ(),zPlaneId));
	fIsSorted = false;
      }
      //Sort vector of clusters by global z-coordinate upon demand.
      void SortClusters() {
	if (fIsSorted == true) return;
	fTrackClusters.clear();
	for (std::map<double,int>::iterator indexIt = fClusterIndicesByZPosition.begin(), 
	       indexEnd = fClusterIndicesByZPosition.end(); indexIt != indexEnd; indexIt++) {
	  fTrackClusters.push_back(indexIt->second);
	  fIsSorted = true;
	}
      }
    
      //Setters.
      void SetIndex(const TrackSegmentIndex& index) { fIndex = index; }
      void SetTrackParameters(const modutils::KalmanTrackParameters& parameters)
      { fTrackParameters = parameters; } 
      
      //Getters.
      const TrackSegmentIndex& GetIndex() const { return fIndex; }
      const std::vector<evt::Index<evt::rec::Cluster> >& GetClusters() const
      { return fTrackClusters; }
      std::vector<evt::Index<evt::rec::Cluster> >& GetClusters()
	{ return fTrackClusters; }
      const modutils::KalmanTrackParameters& GetTrackParameters() const
      { return fTrackParameters; } 
      //Get first/last cluster information according to maps kept sorted by cluster z-position.
      ClusterIndex GetFirstClusterIndex() const
      { return fClusterIndicesByZPosition.begin()->second; }
      ClusterIndex GetLastClusterIndex() const
      { return fClusterIndicesByZPosition.rbegin()->second; }
      double GetFirstClusterZPosition() const
      { return fClusterIndicesByZPosition.begin()->first; }
      double GetLastClusterZPosition() const
      { return fClusterIndicesByZPosition.rbegin()->first; }
      ZPlaneId GetFirstZPlaneId() const
      { return fClusterZPlanesByZPosition.begin()->second; }
      ZPlaneId GetLastZPlaneId() const
      { return fClusterZPlanesByZPosition.rbegin()->second; }
      unsigned int GetNumberOfClusters() const
      { return fTrackClusters.size(); }
      
      //Sort on number of clusters, from largest to smallest.
      bool operator<(const LocalTrackSegment& rhs) const {
	return GetNumberOfClusters() > rhs.GetNumberOfClusters();
      }
      
  private:
      TrackSegmentIndex fIndex;
      //fTrackClusters needs to be a vector of cluster indices in order to play nice with Kalman filters
      //and RecEvent.
      std::vector<evt::Index<evt::rec::Cluster> > fTrackClusters;
      //Track parameters at most-upstream point.
      modutils::KalmanTrackParameters fTrackParameters;
      //Map of cluster indices and corresponding z-positions.
      std::map<double,int> fClusterIndicesByZPosition;
      //Map of cluster indices and corresponding z-planes.
      std::map<double,modutils::ZPlaneId> fClusterZPlanesByZPosition;
      //Flag for sorting clusters by z-position.
      bool fIsSorted;
  };
  
  //Structure for finding local track segments by index (they must be ordered according to position
  //in their STL container!)
  class FindLocalTrackByIndex {

    TrackSegmentIndex index;    

  public:
  FindLocalTrackByIndex(const TrackSegmentIndex& index) : 
    index(index)
    {}

    bool operator()(const LocalTrackSegment& localTrackSegment) const
    { return localTrackSegment.GetIndex() == index; }
  };

  //Class for storing and sorting potential tracklets in minimal tracklet CA mode.
  class MergeCandidate {
    
  public:
  MergeCandidate(const LocalTrackSegment& segment,
		 const double& xDifference,
		 const double& yDifference,
		 const double& aDifference,
		 const double& bDifference,
		 const double& qOverPDifference,
		 const modutils::DetectorSetup& setup) :
    fLocalTrackSegment(segment),
      fXDifference(xDifference),
      fYDifference(yDifference),
      fADifference(aDifference),
      fBDifference(bDifference), 
      fQOverPDifference(qOverPDifference),
      fDetectorSetup(setup)
    //Sort on combination of track parameters. Play with this!
    //Currently involves detector-dependent weights, attempted to bring each term to O(1).
      { 
	if (fDetectorSetup.fQOverPParameterWeight != 0) {
	  fConstructionMetric = sqrt(fDetectorSetup.fXParameterWeight*fXDifference*fXDifference + 
				     fDetectorSetup.fYParameterWeight*fYDifference*fYDifference + 
				     fDetectorSetup.fAParameterWeight*fADifference*fADifference + 
				     fDetectorSetup.fBParameterWeight*fBDifference*fBDifference +
				     fDetectorSetup.fQOverPParameterWeight*qOverPDifference*qOverPDifference); 
	}
	else {
	  fConstructionMetric = sqrt(fDetectorSetup.fXParameterWeight*fXDifference*fXDifference + 
				     fDetectorSetup.fYParameterWeight*fYDifference*fYDifference + 
				     fDetectorSetup.fAParameterWeight*fADifference*fADifference + 
				     fDetectorSetup.fBParameterWeight*fBDifference*fBDifference);
	}
	mergeIndex = fLocalTrackSegment.GetIndex();
	fZPlaneId = fLocalTrackSegment.GetFirstZPlaneId();
      }
    
  MergeCandidate(const LocalTrackSegment& segment,
		 const double xDifference,
		 const double xUncertainty,
		 const double yDifference,
		 const double yUncertainty,
		 const double aDifference,
		 const double aUncertainty,
		 const double bDifference,
		 const double bUncertainty,
		 const double qOverPDifference,
		 const double qOverPUncertainty) :
    fLocalTrackSegment(segment),
      fXDifference(xDifference),
      fYDifference(yDifference),
      fADifference(aDifference),
      fBDifference(bDifference), 
      fQOverPDifference(qOverPDifference)
    //Sort on combination of track parameters. Play with this!
    //Currently involves detector-dependent weights, attempted to bring each term to O(1).
      { 
	if (fDetectorSetup.fQOverPParameterWeight != 0) {
	  fConstructionMetric = sqrt(fXDifference*fXDifference/(xUncertainty*xUncertainty) + 
				     fYDifference*fYDifference/(yUncertainty*yUncertainty) + 
				     fADifference*fADifference/(aUncertainty*aUncertainty) + 
				     fBDifference*fBDifference/(bUncertainty*bUncertainty) +
				     qOverPDifference*qOverPDifference/(qOverPUncertainty*qOverPUncertainty)); 
	}
	else {
	  fConstructionMetric = sqrt(fXDifference*fXDifference/(xUncertainty*xUncertainty) + 
				     fYDifference*fYDifference/(yUncertainty*yUncertainty) + 
				     fADifference*fADifference/(aUncertainty*aUncertainty) + 
				     fBDifference*fBDifference/(bUncertainty*bUncertainty) ); 
	}
	mergeIndex = fLocalTrackSegment.GetIndex();
	fZPlaneId = fLocalTrackSegment.GetFirstZPlaneId();
      }

  MergeCandidate(const evt::Index<evt::rec::Track>& index,
		 const double& xDifference,
		 const double& yDifference,
		 const double& aDifference,
		 const double& bDifference,
		 const double& qOverPDifference,
		 const modutils::DetectorSetup& setup) :
    fRecTrackIndex(index),
      fXDifference(xDifference),
      fYDifference(yDifference),
      fADifference(aDifference),
      fBDifference(bDifference), 
      fQOverPDifference(qOverPDifference),
      fDetectorSetup(setup)
    //Sort on combination of track parameters. Play with this!
    //Currently involves detector-dependent weights, attempted to bring each term to O(1).
      { 
	if (fDetectorSetup.fQOverPParameterWeight != 0) {
	  fConstructionMetric = sqrt(fDetectorSetup.fXParameterWeight*fXDifference*fXDifference + 
				     fDetectorSetup.fYParameterWeight*fYDifference*fYDifference + 
				     fDetectorSetup.fAParameterWeight*fADifference*fADifference + 
				     fDetectorSetup.fBParameterWeight*fBDifference*fBDifference +
				     fDetectorSetup.fQOverPParameterWeight*fQOverPDifference*fQOverPDifference); 
	}
	else {
	  fConstructionMetric = sqrt(fDetectorSetup.fXParameterWeight*fXDifference*fXDifference + 
				     fDetectorSetup.fYParameterWeight*fYDifference*fYDifference + 
				     fDetectorSetup.fAParameterWeight*fADifference*fADifference + 
				     fDetectorSetup.fBParameterWeight*fBDifference*fBDifference);
	}
	mergeIndex = index;
	fZPlaneId = 0;
      }
    
    const TrackIndex& GetIndex() const
    { return mergeIndex; }
    const evt::Index<evt::rec::Track>& GetRecTrackIndex() const
    { return fRecTrackIndex; }
    const LocalTrackSegment& GetLocalTrackSegment() const
    { return fLocalTrackSegment; }
    const modutils::DetectorSetup& GetDetectorSetup() const
    { return fDetectorSetup; }
    const double& GetMetric() const
    { return fConstructionMetric; }
    const double& GetXDifference() const
    { return fXDifference; }
    const double& GetYDifference() const
    { return fYDifference; }
    const double& GetADifference() const
    { return fADifference; }
    const double& GetBDifference() const
    { return fBDifference; }
    const double& GetQOverPDifference() const
    { return fQOverPDifference; }
    const ZPlaneId& GetZPlaneId() const
    { return fZPlaneId; }
    std::string PrintMetric() const {
      std::ostringstream info;
      info << "Metric: " << fConstructionMetric 
	   << ". X part: " << fDetectorSetup.fXParameterWeight*fXDifference*fXDifference
	   << ". Y part: " << fDetectorSetup.fYParameterWeight*fYDifference*fYDifference
	   << ". A part: " << fDetectorSetup.fAParameterWeight*fADifference*fADifference
	   << ". B part: " << fDetectorSetup.fBParameterWeight*fBDifference*fBDifference
	   << ". QOverP part: " << fDetectorSetup.fQOverPParameterWeight*fQOverPDifference*fQOverPDifference;
      return info.str();
    }
    
    void SetZPlaneId(const ZPlaneId& id)
    { fZPlaneId = id; }
    void SetChi2PerNDF(const double chi2PerNDF)
    { fChi2PerNDF = chi2PerNDF; }
    bool operator<(const MergeCandidate& rhs) const {
      return (fConstructionMetric < rhs.fConstructionMetric);
    }

  private:
    LocalTrackSegment fLocalTrackSegment;
    evt::Index<evt::rec::Track> fRecTrackIndex;
    KalmanTrackParameters fTrackParameters;
    double fXDifference;
    double fYDifference;
    double fADifference;
    double fBDifference;
    double fQOverPDifference;
    modutils::DetectorSetup fDetectorSetup;
    double fConstructionMetric;
    double fChi2PerNDF;
    TrackIndex mergeIndex;
    ZPlaneId fZPlaneId;
  };

  typedef std::set<TrackletCandidate> TrackletCandidates;
  //This cannot be a set, since set iterators are inherently const.
  typedef std::list<LocalTrackSegment> LocalTrackSegments;
  //Set for sorting merge candidates.
  typedef std::set<MergeCandidate> MergeCandidates;

}

#endif
