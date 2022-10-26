#include "CellularAutomaton.h"
#include "utl/Line.h"

using namespace std;

namespace modutils {
  
  bool CellularAutomaton::Initialize(const modutils::DetectorSetup& detectorSetup) {
    //Initialization function. Maybe geometry setup, etc. Get xml parameters
    //(detection planes that can be omitted, etc)
    fTracklets.clear();
    fNodes.clear();
    fDetectorId = detectorSetup.fGlobalId;
    fXDisplacementTolerance = detectorSetup.fXDisplacementTolerance;
    fYDisplacementTolerance = detectorSetup.fYDisplacementTolerance;
    fNPlanesToSkip = detectorSetup.fNPlanesToSkip;
    fMinimumClustersOnLocalTrack = detectorSetup.fMinimumClustersOnLocalTrack;

    //Save entire detector setup.
    fDetectorSetup = detectorSetup;

    //Remove tracks found in previous chambers and used cluster list.
    fTrackSegmentFirstClusters.clear();

    //For debugging.
    // cout << "Detector setup (in CA): " << endl;
    // cout << "Global ID: " << detectorSetup.fGlobalId << endl;
    // cout << "DisplacementTolerance: " << detectorSetup.fDisplacementTolerance << endl;
    // cout << "NPlanesToSkip: " << detectorSetup.fNPlanesToSkip << endl;

    return true;
  }
  
  bool CellularAutomaton::PlotTracklets(TFile* file) {
    //Inspect tracklets! Loop over nodes.
    file->cd();
    unsigned int detectorId = (unsigned int)fDetectorId;
    TString nodesString;
    TString trackletsString;

    unsigned int nNodes = 0;
    unsigned int nTracklets = 0;

    nodesString.Form("nodesTPC%i",detectorId);
    trackletsString.Form("trackletsTPC%i",detectorId);
    nodes = new TNtuple(nodesString,nodesString,"detectorId:index:x:y:z:isFirstNode");
    tracklets = new TNtuple(trackletsString,trackletsString,"detectorId:angleXZ:angleYZ:angleXY:firstClusterX:firstClusterY:firstClusterZ:lastClusterX:lastClusterY:lastClusterZ:index:direction");

    for (NodeIterator nodeIt = fNodes.begin(), nodeEnd = fNodes.end(); nodeIt != nodeEnd; nodeIt++) {
      nNodes++;	
      modutils::Node node = nodeIt->second;
      nodes->Fill(fDetectorId,node.GetIndex(),node.GetPosition().GetX(),
		  node.GetPosition().GetY(),node.GetPosition().GetZ(),node.IsFirstNodeOnTrackSegment());
      modutils::TrackletIndices forwardTrackletIndices = node.GetForwardTrackletIndices();
      modutils::TrackletIndices backwardTrackletIndices = node.GetBackwardTrackletIndices();
      for (modutils::TrackletIndices::iterator forwardIt = forwardTrackletIndices.begin(),
	     forwardEnd = forwardTrackletIndices.end(); forwardIt != forwardEnd; forwardIt++) {
	modutils::TrackletIndex forwardTrackletIndex = *forwardIt;
	if (fTracklets.find(forwardTrackletIndex) == fTracklets.end()) {
	  ostringstream msg;
	  msg << "(From PlotTracklets) Didn't get forward index correctly!! Tracklet index: "
	      << forwardTrackletIndex << ", node index: " << node.GetIndex();
	  WARNING(msg);
	  continue;
	}
	modutils::Tracklet tracklet = fTracklets[forwardTrackletIndex];
	nTracklets++;
	if (tracklet.GetConstituentClusters().size() == 0)
	  continue;
	std::vector<modutils::ClusterIndex>::const_iterator clIt = tracklet.GetConstituentClusters().begin();
	utl::Point firstClusterPosition = GetNode(*clIt).GetPosition();
	clIt++;
	utl::Point lastClusterPosition = GetNode(*clIt).GetPosition();
	tracklets->Fill(fDetectorId,tracklet.GetAngleXZ(),tracklet.GetAngleYZ(),
			tracklet.GetAngleXY(),firstClusterPosition.GetX(),firstClusterPosition.GetY(),
			firstClusterPosition.GetZ(),lastClusterPosition.GetX(),lastClusterPosition.GetY(),
			lastClusterPosition.GetZ(),tracklet.GetIndex(),1);
      }
      for (modutils::TrackletIndices::iterator backwardIt = backwardTrackletIndices.begin(),
	     backwardEnd = backwardTrackletIndices.end(); backwardIt != backwardEnd; backwardIt++) {
	modutils::TrackletIndex backwardTrackletIndex = *backwardIt;
	if (fTracklets.find(*backwardIt) == fTracklets.end()) {
	  ostringstream msg;
	  msg << "(From PlotTracklets) Didn't get backward index correctly!! Tracklet index: "
	      << backwardTrackletIndex << ", node index: " << node.GetIndex();
	  WARNING(msg);
	  continue;
	}
	modutils::Tracklet tracklet = fTracklets[*backwardIt];
	nTracklets++;
	if (tracklet.GetConstituentClusters().size() == 0)
	  continue;
	std::vector<modutils::ClusterIndex>::const_iterator clIt = tracklet.GetConstituentClusters().begin();
	utl::Point firstClusterPosition = GetNode(*clIt).GetPosition();
	clIt++;
	utl::Point lastClusterPosition = GetNode(*clIt).GetPosition();
	tracklets->Fill(fDetectorId,tracklet.GetAngleXZ(),tracklet.GetAngleYZ(),
			tracklet.GetAngleXY(),firstClusterPosition.GetX(),firstClusterPosition.GetY(),
			firstClusterPosition.GetZ(),lastClusterPosition.GetX(),lastClusterPosition.GetY(),
			lastClusterPosition.GetZ(),tracklet.GetIndex(),1);
	
      }
    }
    
    file->Write();
    
    return true;
  }
  
  void CellularAutomaton::ApplyYZSlopeRule() {
    
    //Currently, the philosophy for filtering tracklets is to keep ones that fit some criteria.
    //In the future, I want to discard tracks that DON'T fit criteria, so that the filtering 
    //methods can be run iteatively and remove bad tracklets after ones that might have been 
    //making them pass the filtering disappear. 
    modutils::TrackletIndices forwardIndicesToKeep;
    modutils::TrackletIndices backwardIndicesToKeep;

    //Infrastructure for creating subtractive tracklet filtering algorithms.
    modutils::TrackletIndices forwardIndicesToDelete;
    modutils::TrackletIndices backwardIndicesToDelete;
    
    //Loop through all nodes.
    for (NodeIterator nodeIt = fNodes.begin(), nodeEnd = fNodes.end(); nodeIt != nodeEnd; ++nodeIt) {
      modutils::TrackletIndices forwardTrackletIndices = nodeIt->second.GetForwardTrackletIndices();
      modutils::TrackletIndices backwardTrackletIndices = nodeIt->second.GetBackwardTrackletIndices();
      forwardIndicesToKeep.clear();
      backwardIndicesToKeep.clear();
      forwardIndicesToDelete.clear();
      backwardIndicesToDelete.clear();
      
      //Problem here for subtractive algorithm:
      //Don't throw away tracklets if node has zero forward-going or backward-going tracklets --
      //this will kill all tracks starting from upstream/downstream ends of the detector!
      //But! Maybe try putting in a check to allow no backward-going tracklets out of a first padrow
      //and no forward-going tracklets out of a last padrow.
      // if (forwardTrackletIndices.size() == 0 || backwardTrackletIndices.size() == 0) continue;
      //However! This will only work for fully through-going tracks. This means decays and other
      //track topologies will probably be killed...

      //Sandbox for future algorithm development!
      // //Using each forward-going tracklet, search for tracklet with similar slope going backwards.
      // for (modutils::TrackletIndices::iterator forwardIt = forwardTrackletIndices.begin(),
      //forwardEnd = forwardTrackletIndices.end(); forwardIt != forwardEnd; ++forwardIt) {
      // 	modutils::TrackletIndex forwardIndex = *forwardIt;
      // 	//Get forward tracklet slopes.
      // 	double forwardAngleXZ = fTracklets[forwardIndex].GetAngleXZ();
      // 	double forwardAngleYZ = fTracklets[forwardIndex].GetAngleYZ();
      // 	for (modutils::TrackletIndices::iterator backwardIt = backwardTrackletIndices.begin(),
      //backwardEnd = backwardTrackletIndices.end(); backwardIt != backwardEnd; ++backwardIt) {
      // 	  modutils::TrackletIndex backwardIndex = *backwardIt;
      // 	  //Get backward tracklet slopes.
      // 	  double backwardAngleXZ = fTracklets[backwardIndex].GetAngleXZ();
      // 	  double backwardAngleYZ = fTracklets[backwardIndex].GetAngleYZ();
      // 	  Find matching backward track segment. Keep forward and backward tracklets if match is found.
      // 	  bool passedAngleCuts = false;
      // 	  if ( (fabs(forwardAngleYZ - backwardAngleYZ) < fYZAngleTolerance) &&
      //(fabs(forwardAngleXZ - backwardAngleXZ) < fXZAngleTolerance) ) {
      // 	    passedAngleCuts = true;
      // 	  }
      // 	  if (passedAngleCuts) {
      // 	    forwardIndicesToKeep.insert(forwardIndex);
      // 	    backwardIndicesToKeep.insert(backwardIndex);
      // 	  }
      // 	  Infrastructure for subtractive algorithms.
      // 	  else {
      // 	    nodeIt->second.RemoveForwardTrackletIndex(forwardIndex);
      // 	    forwardIndicesToDelete.insert(forwardIndex);
      // 	    nodeIt->second.RemoveBackwardTrackletIndex(backwardIndex);
      // 	    backwardIndicesToDelete.insert(backwardIndex);
      // 	  }
      // 	}
      // }
      // //Purge tracklets with no equal-and-opposite match in YZ-plane or that fail the XZ slope cut.
      // nodeIt->second.KeepTracklets(forwardIndicesToKeep,backwardIndicesToKeep);

      //Infrastructure for subtractive algorithms.
      // RemoveTracklets(forwardIndicesToDelete);
      // RemoveTracklets(backwardIndicesToDelete);
    }
  }
  
  //Function for filtering spurious tracklets. Used in minimal tracklet construction mode.
  void CellularAutomaton::ApplyTrackletPointingRule() {
    //Filtering philosophy: If Node A contains a tracklet pointing to Node B,
    //but Node B does not contain a tracklet pointing to Node A, we delete the tracklet.

    //Loop through all nodes.
    for (NodeIterator nodeIt = fNodes.begin(), nodeEnd = fNodes.end(); nodeIt != nodeEnd; ++nodeIt) {

      //Container for tracklet indices to be deleted.
      modutils::TrackletIndices forwardTrackletIndicesToKeep;
      modutils::TrackletIndices backwardTrackletIndicesToKeep;
      
      //If nodes have more than one forward/backward tracklet, something went wrong in the algorithm
      //or the wrong type of CA tracking was run.
      const modutils::TrackletIndices& forwardTrackletIndices = nodeIt->second.GetForwardTrackletIndices();
      if (forwardTrackletIndices.size() > 1) {
	ERROR("Node has multiple forward tracklets! This filtering rule cannot be applied! "
	      "Use the SlopeFiltering rule instead.");
	return;
      }
      const modutils::TrackletIndices& backwardTrackletIndices = nodeIt->second.GetBackwardTrackletIndices();
      if (backwardTrackletIndices.size() > 1) {
	ERROR("Node has multiple backward tracklets! This filtering rule cannot be applied! "
	      "Use the SlopeFiltering rule instead.");
	return;
      }

      //Get nodes connected by forward-going and backward-going tracklets.
      //Ignore clusters with no tracklets at all.
      if (forwardTrackletIndices.size() == 0 && backwardTrackletIndices.size() == 0)
	continue;
      const modutils::Tracklet& forwardTracklet = fTracklets[*forwardTrackletIndices.begin()];
      const modutils::Tracklet& backwardTracklet = fTracklets[*backwardTrackletIndices.begin()];

      const modutils::Node& forwardNode = fNodes[forwardTracklet.GetLastClusterIndex()];
      const modutils::Node& backwardNode = fNodes[backwardTracklet.GetFirstClusterIndex()];

      //If the forward node doesn't have any backward-going tracklets (as in the case of the
      //downstream edge of a chamber), keep the tracklet.
      const modutils::TrackletIndices& trackletsForwardNode = forwardNode.GetBackwardTrackletIndices();
      if (trackletsForwardNode.size() == 0)
	forwardTrackletIndicesToKeep.insert(forwardTracklet.GetIndex());
      if (SearchTracklet(*trackletsForwardNode.begin())) {
      	//Get index of node pointed to by downstream node's backward-going tracklet.
      	const modutils::ClusterIndex& pointedNodeIndex =
	  fTracklets[(*trackletsForwardNode.begin())].GetFirstClusterIndex();
      	//Compare index of node pointed to with original node's index. If not the same, mark for deletion.
      	if (pointedNodeIndex == nodeIt->second.GetIndex())
	  forwardTrackletIndicesToKeep.insert(forwardTracklet.GetIndex());
      }
      
      //If the backward node doesn't have any forward-going tracklets
      //(as in the case of the upstream edge of a chamber), keep the tracklet.
      const modutils::TrackletIndices& trackletsBackwardNode = backwardNode.GetForwardTrackletIndices();
      if (trackletsBackwardNode.size() == 0)
	backwardTrackletIndicesToKeep.insert(backwardTracklet.GetIndex());
      if (SearchTracklet(*trackletsBackwardNode.begin())) {
      	//Get index of node pointed to by uptream node's forward-going tracklet.
      	const modutils::ClusterIndex& pointedNodeIndex =
	  fTracklets[*trackletsBackwardNode.begin()].GetLastClusterIndex();
      	//Compare index of node pointed to with original node's index. If not the same, mark for deletion.
      	if (pointedNodeIndex == nodeIt->second.GetIndex())
	  backwardTrackletIndicesToKeep.insert(backwardTracklet.GetIndex());
      }
      nodeIt->second.KeepTracklets(forwardTrackletIndicesToKeep,backwardTrackletIndicesToKeep);
    }
  }

  //Function for finding most-upstream clusters on each track segment. Creates list of these
  //for joining segments and fitting. 
  void CellularAutomaton::FindFirstClusterOnTrackSegments() {

    //This function should only run after track candidates have been found. Make sure this has been done.
    if (fTrackCandidates.size() == 0) {
      //WARNING("Track candidates have not been calculated, or there were zero found in this chamber!");
      return;
    }

    //Clear list of first clusters and track segment indices.
    fNTrackSegments = 0;

    for (TrackCandidateIterator candidateIt = fTrackCandidates.begin(),
	   candidateEnd = fTrackCandidates.end(); candidateIt != candidateEnd; candidateIt++) {
      //Get node list for each track candidate.
      const NodeList candidateNodes = (*candidateIt).GetTrackNodes();
      NodeList::iterator firstNode = candidateNodes.begin();
      ClusterIndex firstNodeIndex = (*firstNode).GetIndex();
      //Nodes are organized with ascending index as we move downstream.
      //The lowest-valued index will be the most upstream node on the track segment. AAHH NOT ALWAYS TRUE!!
      //Loop through nodes and mark node index with smallest z-value.
      //      for (NodeList::iterator nodeIt = candidateNodes.begin(), nodeEnd = candidateNodes.end();
      //nodeIt != nodeEnd; nodeIt++) {
      //	if (fNodes[*nodeIt.GetIndex()].GetPosition().GetZ() <
      //            fNodes[firstNodeIndex].GetPosition().GetZ()) firstNodeIndex = *nodeIt;
      //      }
      //Sanity check to make sure there are no unwanted links.
      if (fNodes[firstNodeIndex].GetBackwardTrackletIndices().size() == 0) {
	fNodes[firstNodeIndex].SetFirstNodeOnTrackSegment(true);
	//Add to list of first clusters. Make sure index matches track candidate's index.
	fTrackSegmentFirstClusters.insert(MakeNodeAndTrackSegmentPair(firstNodeIndex,
								      (*candidateIt).GetIndex()));
      }
    }
    //Set total number of track segments.
    fNTrackSegments = fTrackSegmentFirstClusters.size();
  }
  
  //Function to play around with applying more cuts. Change to whatever you'd like to try!
  void CellularAutomaton::ApplyCut(utl::Point targetPosition) {
    
    modutils::TrackletIndices forwardIndicesToKeep;
    modutils::TrackletIndices backwardIndicesToKeep;
    modutils::TrackletIndices forwardIndicesToDelete;
    modutils::TrackletIndices backwardIndicesToDelete;
    
    //Loop through all nodes.
    for (NodeIterator nodeIt = fNodes.begin(), nodeEnd = fNodes.end(); nodeIt != nodeEnd; ++nodeIt) {
      modutils::TrackletIndices forwardTrackletIndices = nodeIt->second.GetForwardTrackletIndices();
      modutils::TrackletIndices backwardTrackletIndices = nodeIt->second.GetBackwardTrackletIndices();
      forwardIndicesToKeep.clear();
      backwardIndicesToKeep.clear();
      forwardIndicesToDelete.clear();
      backwardIndicesToDelete.clear();
      //Using each tracklet, calculate the extrapolated track position at the target.
      for (modutils::TrackletIndices::iterator forwardIt = forwardTrackletIndices.begin(),
	     forwardEnd = forwardTrackletIndices.end(); forwardIt != forwardEnd; ++forwardIt) {
	modutils::TrackletIndex forwardIndex = *forwardIt;
	bool passedForwardCut = false;
	//If the extrapolated z-position is more than .5 meters from the target, reject tracklet.
	utl::Line trackletLine((nodeIt->second).GetPosition(),fTracklets[forwardIndex].GetAngleTangents());
	utl::Point extrapolatedTrackletPosition =
	  trackletLine.Propagate((targetPosition - (nodeIt->second).GetPosition()).GetZ()/
				 trackletLine.GetDirection().GetZ());
	if ( fabs((extrapolatedTrackletPosition-targetPosition).GetY()) < 50*utl::cm)
	  passedForwardCut = true;
	if (passedForwardCut)
	  forwardIndicesToKeep.insert(forwardIndex);
	else {
	  forwardIndicesToDelete.insert(forwardIndex);
	}
      }
      for (modutils::TrackletIndices::iterator backwardIt = backwardTrackletIndices.begin(),
	     backwardEnd = backwardTrackletIndices.end(); backwardIt != backwardEnd; ++backwardIt) {
	modutils::TrackletIndex backwardIndex = *backwardIt;
	bool passedBackwardCut = false;	  
	//If the extrapolated z-position is more than .5 meters from the target, reject tracklet.
	//Don't forget to search for tracklet by index on the previous z-plane!
	utl::Line trackletLine((nodeIt->second).GetPosition(),fTracklets[backwardIndex].GetAngleTangents());
	utl::Point extrapolatedTrackletPosition =
	  trackletLine.Propagate((targetPosition - (nodeIt->second).GetPosition()).GetZ()/
				 trackletLine.GetDirection().GetZ());
	if ( fabs((extrapolatedTrackletPosition-targetPosition).GetY()) < 50*utl::cm)
	  passedBackwardCut = true;
	if (passedBackwardCut)
	  backwardIndicesToKeep.insert(backwardIndex);
	else {
	  backwardIndicesToDelete.insert(backwardIndex);
	}
      }
      //Purge tracklets that don't pass cuts.
      nodeIt->second.KeepTracklets(forwardIndicesToKeep,backwardIndicesToKeep);
      // RemoveTracklets(forwardIndicesToDelete);
      // RemoveTracklets(backwardIndicesToDelete);
    }
  }
  
  void CellularAutomaton::FindTrackCandidates() {
    //Clear out any previous candidates.
    fTrackCandidates.clear();
    fUsedClusters.clear();
    //Resize used cluster list to total number of clusters in event.
    if (fNodes.size() == 0)
      return;
    const ClusterIndex& largestNodeIndex = fNodes.rbegin()->first;
    fUsedClusters.resize(largestNodeIndex + 1,false);

    //Loop through all nodes.
    modutils::TrackSegmentIndex index = 0;
    for (NodeIterator nodeIt = fNodes.begin(), nodeEnd = fNodes.end(); nodeIt != nodeEnd; ++nodeIt) {
      const modutils::Node& node = nodeIt->second;
      //Check to see if cluster has already been considered.
      if (fUsedClusters[node.GetIndex()])
	continue;
      //Check to see if cluster has any associated trackelts.
      if (node.GetForwardTrackletIndices().size() == 0 && node.GetBackwardTrackletIndices().size() == 0)
	continue;
      //Now, begin a new collection of clusters if node has tracklets associated with it.
      modutils::TrackCandidate trackCandidate(node);
      trackCandidate.SetIndex(index);
      index++;
      fTrackCandidates.push_back(trackCandidate);
      if ( node.GetIndex() > fUsedClusters.size() )
	fUsedClusters.resize(node.GetIndex(),false);
      fUsedClusters[node.GetIndex()] = true;
      //Find other nodes connected to this node by its constituent tracklets.
      //Add them to this list of clusters.
      FindConnectedClusters(fTrackCandidates.back(), node);
    }
  }

  //Recursive function to find all clusters connected to a starting cluster by tracklets.
  void CellularAutomaton::FindConnectedClusters(modutils::TrackCandidate& trackCandidate,
						const modutils::Node& node) {

    const modutils::TrackletIndices& forwardTrackletIndices = node.GetForwardTrackletIndices();
    const modutils::TrackletIndices& backwardTrackletIndices = node.GetBackwardTrackletIndices();
    for (modutils::TrackletIndices::iterator forwardIt = forwardTrackletIndices.begin(),
	   forwardEnd = forwardTrackletIndices.end(); forwardIt != forwardEnd; ++forwardIt) {
      //Get tracklets and constituent clusters. 
      //For forward tracklets, the first node belongs to the node we're considering.
      //We only need the last node.
      const modutils::Tracklet& forwardTracklet = fTracklets[*forwardIt];
      const modutils::ClusterIndex& lastIndex = forwardTracklet.GetLastClusterIndex();
      const modutils::Node& lastNode = fNodes[lastIndex];
      //Apply track candidate check. Add node to track if it passes check.
      //For now: don't implement check.
      if (fUsedClusters[lastNode.GetIndex()] == false) {
	trackCandidate.AddNodeToTrack(lastNode,forwardTracklet);
	//clusterList.push_back(lastNode.GetIndex());
	fUsedClusters[lastNode.GetIndex()] = true;
	FindConnectedClusters(trackCandidate,lastNode);
      }
    }
    for (modutils::TrackletIndices::iterator backwardIt = backwardTrackletIndices.begin(),
	   backwardEnd = backwardTrackletIndices.end(); backwardIt != backwardEnd; ++backwardIt) {
      //Get tracklets and constituent clusters. 
      //For backward tracklets, the last node belongs to the node we're considering.
      //We only need the first node.
      const modutils::Tracklet& backwardTracklet = fTracklets[*backwardIt];
      const modutils::ClusterIndex& firstIndex = backwardTracklet.GetFirstClusterIndex();
      const modutils::Node& firstNode = fNodes[firstIndex];
      //Apply track candidate check. Add node to track if it passes check.
      //For now: don't implement check.
      if (fUsedClusters[firstNode.GetIndex()] == false) {
	trackCandidate.AddNodeToTrack(firstNode,backwardTracklet);
	//clusterList.push_back(firstNode.GetIndex());
	fUsedClusters[firstNode.GetIndex()] = true;
	FindConnectedClusters(trackCandidate,firstNode);
      }
    }
  }

  void CellularAutomaton::FilterTrackCandidates() {
    //This function should only run after track candidates have been found. Make sure this has been done.
    if (fTrackCandidates.size() == 0) {
      //WARNING("Track candidates have not been calculated, or there were zero found in this chamber!");
      return;
    }

    //Loop through tracklets in this candidate. Filter according to first-grade outlier test. TODO!!
    for (TrackCandidateIterator candidateIt = fTrackCandidates.begin(), candidateEnd = fTrackCandidates.end(); 
	 candidateIt != candidateEnd; candidateIt++) {
      modutils::TrackCandidate& candidate = *candidateIt;
      
      //Length check for track candidates.
      if (candidate.GetNClusters() < fMinimumClustersOnLocalTrack)
	continue;
      
      //Calculate relevant tracklet parameters.
      candidate.CalculateAngleStatistics();
      // double xzMedian = candidate.GetXZAngleMedian();
      // double xzInterQuartileRange = candidate.GetXZAngleInterQuartileRange();
      double yzMedian = candidate.GetYZAngleMedian();
      double yzInterQuartileRange = candidate.GetYZAngleInterQuartileRange();
      
      //Container for tracklets to delete.
      modutils::TrackletIndices trackletIndicesToDelete;
      
      //Loop through tracklets and get rid of outliers.
      const modutils::TrackletVector& tracklets = candidate.GetTracklets();
      for (modutils::TrackletVector::const_iterator trackletIt = tracklets.begin(),
	     trackletEnd= tracklets.end(); trackletIt != trackletEnd; trackletIt++) {
	// double angleXZ = trackletIt->GetAngleXZ();
	double angleYZ = trackletIt->GetAngleYZ();
	
	//Cuts. Outlier is defined as more than 1.5 * inter-quartile range from the median.
	if (fabs(angleYZ - yzMedian) > 1.5*yzInterQuartileRange)
	  trackletIndicesToDelete.insert(trackletIt->GetIndex());
	
      }
      
      //Remove bad tracklets.
      if (trackletIndicesToDelete.size() > 0)
	RemoveTracklets(trackletIndicesToDelete);
      
    }
  }

  
  void CellularAutomaton::AddForwardTrackletToNode(const modutils::ClusterIndex& index,
						   const modutils::Tracklet& tracklet) {
    //Add a (forward-going) tracklet to an existing node.
    fNodes[index].AddForwardTrackletIndex(tracklet.GetIndex());
  }
  void CellularAutomaton::AddBackwardTrackletToNode(const modutils::ClusterIndex& index,
						    const modutils::Tracklet& tracklet) {
    //Add a (backward-going) tracklet to an existing node.
    fNodes[index].AddBackwardTrackletIndex(tracklet.GetIndex());
  }
  
  void CellularAutomaton::RemoveTracklets(modutils::TrackletIndices trackletIndicesToDelete) {
    //Remove tracklets from event. Search for them first, to be sure we haven't already deleted them.
    for (modutils::TrackletIndices::iterator trackletIndexIt = trackletIndicesToDelete.begin(),
	   trackletIndexEnd = trackletIndicesToDelete.end();
	 trackletIndexIt != trackletIndexEnd; ++trackletIndexIt) {
      TrackletIterator trackletIt = fTracklets.find(*trackletIndexIt);
      if (trackletIt == fTracklets.end()) {
	cout << "(From RemoveTracklets): Didn't get the index correctly!! Tracklet index: "
	     << *trackletIndexIt << endl;
	continue;
      }
      //Erase index from forward and backward node tracklet lists.
      modutils::Tracklet& tracklet = trackletIt->second;
      
      NodeIterator firstNodeIt = fNodes.find(tracklet.GetFirstClusterIndex());
      if (firstNodeIt == fNodes.end()) {
	cout << "(From RemoveTracklets): Tracklet index: "
	     << *trackletIndexIt
	     << " Error removing tracklet from node "
	     << tracklet.GetFirstClusterIndex() << endl;
	continue;
      }
      NodeIterator lastNodeIt = fNodes.find(tracklet.GetLastClusterIndex());
      if (lastNodeIt == fNodes.end()) {
	cout << "(From RemoveTracklets): Tracklet index: "
	     << *trackletIndexIt
	     << " Error removing tracklet from node "
	     << tracklet.GetLastClusterIndex() << endl;
	continue;
      }
      fNodes[tracklet.GetFirstClusterIndex()].RemoveForwardTrackletIndex(*trackletIndexIt);
      fNodes[tracklet.GetLastClusterIndex()].RemoveBackwardTrackletIndex(*trackletIndexIt);
      fTracklets.erase(trackletIt);
    }
  }
  
  void CellularAutomaton::AddTracklet(const modutils::Tracklet& tracklet) {
    //Fill cellular automaton with all tracklets.
    //Make a std::pair consisting of trackletIndex (key) and tracklet and add to the tracklet map.
    fTracklets.insert(std::pair<modutils::TrackletIndex,modutils::Tracklet>(tracklet.GetIndex(),tracklet));
  }
  
  bool CellularAutomaton::SearchTracklet(const modutils::TrackletIndex& index) {
    //Search for a tracklet by index. Return true if it exists.
    TrackletIterator it = fTracklets.find(index);
    return (it != fTracklets.end());
  }

  void CellularAutomaton::AddNode(modutils::Node& node) {
    //Fill cellular automaton with all nodes.
    //Make a std::pair consisting of clusterIndex (key) and node and add to the node map.
    fNodes.insert(std::pair<modutils::ClusterIndex,modutils::Node>(node.GetIndex(),node));
  }
  
  bool CellularAutomaton::SearchNode(const modutils::ClusterIndex& index) {
    //Search for a node on a z-plane by cluster index. Return true if it exists.
    NodeIterator it = fNodes.find(index);
    return (it != fNodes.end());
  }

  modutils::Node& CellularAutomaton::GetNode(const ClusterIndex& index) {
    NodeIterator nodeIt = fNodes.find(index);
    if (nodeIt == fNodes.end()) {
      ostringstream err;
      err << "Couldn't find a node with index " << index
	  << "! Node should exist!! SearchNode(index): " << SearchNode(index)
	  << ". fNodes.size(): " << fNodes.size()
	  << ". Existing nodes: ";
      // for (NodeIterator nodeIt = fNodes.begin(), nodeEnd = fNodes.end(); nodeIt != nodeEnd; ++nodeIt) {
      // 	err << endl << nodeIt->first;
      // }
      ERROR(err);
    }
    return nodeIt->second;
  }
  
  void CellularAutomaton::SetNodeUsed(const ClusterIndex& index, const bool isUsed) {
    NodeIterator nodeIt = fNodes.find(index);
    if (nodeIt == fNodes.end()) {
      ostringstream err;
      err << "Couldn't find a node with index" << index << "! Node should exist!!";
      ERROR(err);
    }
    (nodeIt->second).SetUsed(isUsed);
  }
}
