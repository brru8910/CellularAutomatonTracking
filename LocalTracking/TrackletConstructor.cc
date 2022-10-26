#include "TrackletConstructor.h"

using namespace std;
using namespace evt;
using namespace evt::proc;

namespace modutils {
  
  bool TrackletConstructor::ConstructTracklets(evt::proc::Reconstruction& procReconstruction,
					       const modutils::ClusterPositionsMap& clusterPositionsByIndex,
					       modutils::CellularAutomaton& cellularAutomaton,
					       double areaSearchScaleFactor) {
    //Get detector fDetectorSetup parameters.
    //TPC-specific. Make this applicable to VD and other detectors also!!
    const det::TPCConst::EId& chamberId =
      det::TPCConst::GetId(det::Const::GetName(fDetectorSetup.fGlobalId));
    const unsigned int zPlaneSeparation = 1;
    //We need to reduce the displacement tolerance here because we're only using clusters on adjacent padrows.
    //Tolerance can be adjusted externally so that the algorithm can be ran iteratively.
    const double xDisplacementTolerance = areaSearchScaleFactor*(fDetectorSetup.fXDisplacementTolerance);    
    const double yDisplacementTolerance = areaSearchScaleFactor*(fDetectorSetup.fYDisplacementTolerance);
    const double maxAngleChangeXZ = fDetectorSetup.fMaxAngleChangeXZ;
    const double maxAngleChangeYZ = fDetectorSetup.fMaxAngleChangeYZ;
    
    // Loop over detection planes.
    for (unsigned int zPlaneId = fDetectorSetup.fZPlaneIdRange.GetBegin(),
	   zPlaneEnd = fDetectorSetup.fZPlaneIdRange.GetEnd();
	 zPlaneId < zPlaneEnd; ++zPlaneId) {
      //Get this Z-Plane's detection plane grid.
      DetectionPlaneGrid& grid =
	procReconstruction.GetDetectionPlaneGrid(zPlaneId);
      // const utl::Point& gridPoint = grid.GetZPlane().GetAnchor();
      const double gridZ = grid.GetGridZ();
      
      //Get list of active clusters from grid structure.
      const ClusterVector& activeClusters = grid.GetActiveClusters();
      for (ClusterVector::const_iterator it = activeClusters.begin(), end = activeClusters.end();
	   it != end ; ++it) {
	const modutils::ClusterIndex& clusterIndex = *it;
	const utl::Point& clusterPosition = clusterPositionsByIndex.at(clusterIndex);

	//Container for filling node info.
	modutils::Node node;
	//Check if node corresponding to cluster exists.	      
	const bool nodeExists = cellularAutomaton.SearchNode(clusterIndex);
	if (!nodeExists) {
	  //If it doesn't, create it and prepare for adding tracklets (will come later).
	  node.SetIndex(clusterIndex);
	  node.SetZPlaneId(zPlaneId);
	  node.SetPosition(clusterPosition);
	  cellularAutomaton.AddNode(node);
	}
	else node = cellularAutomaton.GetNode(clusterIndex);

	//If node already has been used (contains best tracklets), skip it.
	if (node.IsUsed()) continue;

	//Only add nodes for first/last z-plane -- no need to look for tracklets.
	if (zPlaneId == fDetectorSetup.fZPlaneIdRange.GetBegin() ||
	    zPlaneId == fDetectorSetup.fZPlaneIdRange.GetEnd() - 1) continue;
	
	//Make a list of potential tracklets .
	modutils::TrackletCandidates trackletCandidates;
	const modutils::ZPlaneId nextZPlaneId = zPlaneId + zPlaneSeparation;
	const modutils::ZPlaneId previousZPlaneId = zPlaneId - zPlaneSeparation;
	
	if (previousZPlaneId < fDetectorSetup.fZPlaneIdRange.GetBegin()) {
	  ERROR("Previous zplane ID is before the chamber range!!");
	  continue;
	}
	if (nextZPlaneId > fDetectorSetup.fZPlaneIdRange.GetEnd()) {
	  ERROR("Next zplane ID is past the chamber range!!");
	  continue;
	}
	//Get next layer's cluster indices, using grid structure to only get clusters near (x,y) points
	//we're interested in.
	const DetectionPlaneGrid& nextLayerGrid =
	  procReconstruction.GetDetectionPlaneGrid(nextZPlaneId);

	//Calculate z-displacement between clusters in order to scale cluster distance cuts.
	const double nextGridZ = nextLayerGrid.GetGridZ();
	const double downstreamDisplacementZ = nextGridZ - gridZ;

	//Get subset of cluster indices that are within acceptable distance from our central cluster.
	const ClusterVector& clusterIndicesNextDetectionLayer =
	  nextLayerGrid.GetClustersFromPositionAndRange(clusterPosition.GetX(),
							clusterPosition.GetY(),
							xDisplacementTolerance*downstreamDisplacementZ,
							yDisplacementTolerance*downstreamDisplacementZ);

	for (ClusterVector::const_iterator nextLayerIt = clusterIndicesNextDetectionLayer.begin(), 
	       nextLayerEnd = clusterIndicesNextDetectionLayer.end();
	     nextLayerIt != nextLayerEnd; ++nextLayerIt) {
	  
	  const modutils::ClusterIndex& clusterIndexNextLayer = *nextLayerIt;

	  const utl::Point& clusterPositionNextLayer = clusterPositionsByIndex.at(clusterIndexNextLayer);
	  //Put in some sanity cuts for forming tracklets (clusters must be inside TPCs),
	  //and also some cuts to keep searches local (computation time down!).
	  const double downstreamDisplacementX = clusterPositionNextLayer.GetX() - clusterPosition.GetX();
	  const double downstreamDisplacementY = clusterPositionNextLayer.GetY() - clusterPosition.GetY();
	  if ( fabs(downstreamDisplacementX) > xDisplacementTolerance*fabs(downstreamDisplacementZ) || 
	       fabs(downstreamDisplacementY) > yDisplacementTolerance*fabs(downstreamDisplacementZ) )
	    continue;

	  //Construct forward tracklet (for obtaning tracklet angles).
	  const modutils::Tracklet forwardTracklet(clusterIndex,clusterPosition,clusterIndexNextLayer,
						   clusterPositionNextLayer,0,chamberId);
	  const double forwardXZAngle = forwardTracklet.GetAngleXZ();
	  const double forwardYZAngle = forwardTracklet.GetAngleYZ();

	  //Now loop through nodes on previous detection plane for backward tracklet angle calculation.
	  const DetectionPlaneGrid& previousLayerGrid =
	    procReconstruction.GetDetectionPlaneGrid(previousZPlaneId);
	  
	  //Calculate z-displacement between clusters in order to scale cluster distance cuts.
	  const double previousGridZ = previousLayerGrid.GetGridZ();
	  const double upstreamDisplacementZ = gridZ - previousGridZ;

	  //Get subset of cluster indices that are within acceptable distance from our central cluster.
	  const ClusterVector& clusterIndicesPreviousDetectionLayer =
	    previousLayerGrid.GetClustersFromPositionAndRange(clusterPosition.GetX() - downstreamDisplacementX,
							      clusterPosition.GetY() - downstreamDisplacementY,
							      xDisplacementTolerance*upstreamDisplacementZ,
							      yDisplacementTolerance*upstreamDisplacementZ);
	  for (ClusterVector::const_iterator previousLayerIt =
		 clusterIndicesPreviousDetectionLayer.begin(), 
		 previousLayerEnd = clusterIndicesPreviousDetectionLayer.end();
	       previousLayerIt != previousLayerEnd; ++previousLayerIt) {
	    
	    const modutils::ClusterIndex& clusterIndexPreviousLayer = *previousLayerIt;
	    const utl::Point& clusterPositionPreviousLayer =
	      clusterPositionsByIndex.at(clusterIndexPreviousLayer);

	    const double upstreamDisplacementX = clusterPosition.GetX() - clusterPositionPreviousLayer.GetX();
	    const double upstreamDisplacementY = clusterPosition.GetY() - clusterPositionPreviousLayer.GetY();
	    if (fabs(upstreamDisplacementX) > xDisplacementTolerance*fabs(upstreamDisplacementZ) || 
		fabs(upstreamDisplacementY) > yDisplacementTolerance*fabs(upstreamDisplacementZ)  )
	      continue;
	    
	    //Construct backward tracklet (for obtaning tracklet angles).
	    const modutils::Tracklet backwardTracklet(clusterIndexPreviousLayer,clusterPositionPreviousLayer,
						      clusterIndex,clusterPosition,0,chamberId);
	    const double backwardXZAngle = backwardTracklet.GetAngleXZ();
	    const double backwardYZAngle = backwardTracklet.GetAngleYZ();

	    //Find the combination of forward/backward tracklets that gives the smallest change in slope.
	    //To do this, we fill a std::set of candidate objects containing the changes in slope, the cluster
	    //indices, and is sorted by smallest change in slope -> largest change in slope.
	    const double xzAngleDifference = forwardXZAngle - backwardXZAngle;
	    const double yzAngleDifference = forwardYZAngle - backwardYZAngle;
	    
	    if (fabs(xzAngleDifference) > maxAngleChangeXZ || 
		fabs(yzAngleDifference) > maxAngleChangeYZ  ) continue;

	    const modutils::TrackletCandidate candidate(xzAngleDifference,yzAngleDifference,
							clusterIndexPreviousLayer,clusterIndex,
							clusterIndexNextLayer,zPlaneId);

	    trackletCandidates.insert(candidate);
	  } //End previous plane loop
	} //End next plane loop
	//If no pair was found, move on.
	if (trackletCandidates.size() == 0)  { 
	  continue;
	}

	//Get the best candidate pair.
	const modutils::TrackletCandidate& bestCandidate = *(trackletCandidates.begin());
	//Now create forward and backward tracklets from best candidate.
	const modutils::ClusterIndex& clusterIndexNextLayer = bestCandidate.GetClusterIndexNextLayer();
	const utl::Point& clusterPositionNextLayer = clusterPositionsByIndex.at(clusterIndexNextLayer);
	const modutils::ClusterIndex& clusterIndexPreviousLayer = bestCandidate.GetClusterIndexPreviousLayer();
	const utl::Point& clusterPositionPreviousLayer = clusterPositionsByIndex.at(clusterIndexPreviousLayer);

	const modutils::Tracklet bestForwardTracklet(clusterIndex,clusterPosition,clusterIndexNextLayer,
						     clusterPositionNextLayer,fTrackletIndex,chamberId);
	cellularAutomaton.AddTracklet(bestForwardTracklet);
	node.AddForwardTrackletIndex(fTrackletIndex);		  
	fTrackletIndex++;
	const modutils::Tracklet bestBackwardTracklet(clusterIndexPreviousLayer,clusterPositionPreviousLayer,
						      clusterIndex,clusterPosition,fTrackletIndex,chamberId);
	cellularAutomaton.AddTracklet(bestBackwardTracklet);
	node.AddBackwardTrackletIndex(fTrackletIndex);		  
	fTrackletIndex++;
	
	//Add node with selected tracklets to CA and set node as used.
	cellularAutomaton.AddForwardTrackletToNode(clusterIndex,bestForwardTracklet);
	cellularAutomaton.AddBackwardTrackletToNode(clusterIndex,bestBackwardTracklet);
	cellularAutomaton.SetNodeUsed(clusterIndex,true);
	//Remove cluster from detection plane (for speed increase!)
	grid.AddClusterToDisabledList(clusterIndex);

      } //End tracklet-making loop
      grid.DisableClusters();
    } //End z-plane loop

    return true;
  }
  
}
