#include "TrackSegmentExtender.h"

#include <utl/Point.h>

#include "TNtuple.h"

using namespace std;
using namespace evt;
using namespace evt::proc;


namespace modutils {
  
  bool TrackSegmentExtender::MergeTrackSegments(evt::proc::Reconstruction& procReconstruction,
						modutils::LocalTrackSegments& localTrackSegments,
						const det::TPCChamber& chamber) {
    //Sort vector of local track segments before working with
    //them. This ensures that they'll be ordered by number of
    //clusters.
    localTrackSegments.sort();
    
    //Container for holding track segments we've merged.
    set<modutils::TrackSegmentIndex> segmentsToErase;
    
    //Start with track with the most clusters.
    for (modutils::LocalTrackSegments::iterator localTrackSegmentIt = localTrackSegments.begin(), 
  	   localTrackSegmentEnd = localTrackSegments.end();
  	 localTrackSegmentIt != localTrackSegmentEnd; ++localTrackSegmentIt) {

      //Get local track segment for modification.
      modutils::LocalTrackSegment& localTrackSegment = *localTrackSegmentIt;

      if (segmentsToErase.find(localTrackSegmentIt->GetIndex()) != segmentsToErase.end())
      	continue;
      
      //Create a list of merge candidates for this track segment.
      modutils::MergeCandidates mergeCandidates;
      
      //Record track's first and last z-plane ID.
      const modutils::ZPlaneId& trackFirstZPlane = localTrackSegment.GetFirstZPlaneId();
      const modutils::ZPlaneId& trackLastZPlane = localTrackSegment.GetLastZPlaneId();
      
      //Some track segments overlap on z-planes due to cluster splitting. To merge these
      //tracks, we begin the merge search just before the end of the track.
      for (unsigned int nextZPlane = fDetectorSetup.fZPlaneIdRange.GetBegin(),
  	     lastZPlane = fDetectorSetup.fZPlaneIdRange.GetEnd();
  	   nextZPlane < lastZPlane; ++nextZPlane ) {
	
  	//Skip z-planes spanned by the track segment.
  	if (nextZPlane > trackFirstZPlane &&
  	    nextZPlane < trackLastZPlane)
  	  continue;

  	proc::DetectionPlaneGrid& nextLayerGrid = 
	  procReconstruction.GetDetectionPlaneGrid(nextZPlane);
	
  	//Extrapolate track to next detection plane and calculate
  	//allowed (x,y) window for clusters to be in for merging.
  	const double nextGridPlaneZ = nextLayerGrid.GetGridZ();
  	modutils::KalmanTrackParameters extrapolatedTrackParameters;

	//Test for MTPCL vs MTPCR.
	if (fDetectorSetup.fFitterType == modutils::eKalmanFilter) {
	  fKalmanFilter.ExtrapolateToGlobalZ(localTrackSegment.GetTrackParameters(),
					     extrapolatedTrackParameters,
					     nextGridPlaneZ);
  	}
  	if (fDetectorSetup.fFitterType == modutils::eStraightTrackFitter) { 
 	  extrapolatedTrackParameters = localTrackSegment.GetTrackParameters();
  	  extrapolatedTrackParameters.ExtrapolateStraightFitToZ(nextGridPlaneZ);
  	}	
  	if (fDetectorSetup.fFitterType == modutils::eParabolicTrackFitter) {
  	  extrapolatedTrackParameters = localTrackSegment.GetTrackParameters();
  	  extrapolatedTrackParameters.ExtrapolateParabolicFitToZ(nextGridPlaneZ);
  	}	

	//Get plane object corresponding to z-plane. If TPC Id is MTPCL and
	//cluster X is negative, the cluster is in MTPCR.
	const utl::Plane& zPlane = (chamber.GetId() == det::TPCConst::eMTPCL &&
				    extrapolatedTrackParameters.GetPosition().GetX() < 0) ?
	  fMTPCR->GetZPlane(nextZPlane) :
	  chamber.GetZPlane(nextZPlane);
  	if (fDetectorSetup.fFitterType == modutils::eKalmanFilter) {
  	  fKalmanFilter.ExtrapolateToPlane(localTrackSegment.GetTrackParameters(),
					   extrapolatedTrackParameters,
					   zPlane);
  	}
  	if (fDetectorSetup.fFitterType == modutils::eStraightTrackFitter) {
  	  extrapolatedTrackParameters = localTrackSegment.GetTrackParameters();
  	  extrapolatedTrackParameters.ExtrapolateStraightFitToPlane(zPlane);
  	}	
  	if (fDetectorSetup.fFitterType == modutils::eParabolicTrackFitter) {
  	  extrapolatedTrackParameters = localTrackSegment.GetTrackParameters();
  	  extrapolatedTrackParameters.ExtrapolateParabolicFitToPlane(zPlane);
  	}	
	
  	const utl::Point& extrapolatedPosition = extrapolatedTrackParameters.GetPosition();

  	//Make sure we haven't left the chamber. MTPCL requires a check of MTPCR, since they
  	//are treated at the same time.
  	if (chamber.GetId() == det::TPCConst::eMTPCL) {
  	  if (!chamber.IsIn(extrapolatedPosition) &&
  	      !fMTPCR->IsIn(extrapolatedPosition) ) 
  	    continue;
  	}
  	else 
  	  if (!chamber.IsIn(extrapolatedPosition)) 
  	    continue;

  	//Check if cluster/track is located in acceptance window (including errors).
  	//Don't let the window get much larger than the maximum allowed cluster-track residual.
	const double maxResidual = fDetectorSetup.fMaxAverageResidual;
  	const double searchRangeX = 10*maxResidual;
  	const double searchRangeY = 10*maxResidual;
	
  	//Record extrapolated parameters.
  	const double trackX = extrapolatedTrackParameters.GetX();
  	const double trackY = extrapolatedTrackParameters.GetY();
  	const double trackA = extrapolatedTrackParameters.GetA();
  	const double trackB = extrapolatedTrackParameters.GetB();
  	const double trackQOverP = extrapolatedTrackParameters.GetQOverP();
	
  	//Get subset of track indices that are within acceptable distance from our central cluster.
  	const proc::TrackSegmentIndexVector& trackIndicesNextDetectionLayer =
  	  nextLayerGrid.GetTrackSegmentsFromPositionAndRange(extrapolatedPosition.GetX(),
  							     extrapolatedPosition.GetY(),
  							     searchRangeX,searchRangeY);
		  
  	//Loop through track segments on next detection plane.
  	for (proc::TrackSegmentIndexVector::const_iterator mergeTrackIt =
  	       trackIndicesNextDetectionLayer.begin(), mergeTrackEnd =
	       trackIndicesNextDetectionLayer.end();
  	     mergeTrackIt != mergeTrackEnd; mergeTrackIt++) {

	  //Skip segments we've merged.
	  if (segmentsToErase.find(*mergeTrackIt) != segmentsToErase.end())
	    continue;

	  //Get the other local track.
  	  const modutils::LocalTrackSegments::iterator pointedTrackIt =
  	    std::find_if(localTrackSegments.begin(), localTrackSegments.end(),
  			 modutils::FindLocalTrackByIndex(*mergeTrackIt));
  	  if (pointedTrackIt == localTrackSegments.end())
  	    continue;
  	  const modutils::LocalTrackSegment& pointedTrack = *pointedTrackIt;

  	  //Don't try to merge this track with itself.
  	  if (localTrackSegment.GetIndex() == pointedTrack.GetIndex())
  	    continue;

  	  modutils::KalmanTrackParameters pointedTrackParameters = pointedTrack.GetTrackParameters();

  	  const utl::Vector& residual = extrapolatedPosition - pointedTrackParameters.GetPosition();
  	  if (fabs(residual.GetX()) > searchRangeX ||
  	      fabs(residual.GetY()) > searchRangeY  )
  	    continue;
	  
  	  const double differenceX = trackX - pointedTrackParameters.GetX();
  	  const double differenceY = trackY - pointedTrackParameters.GetY();
  	  const double differenceA = trackA - pointedTrackParameters.GetA();
  	  const double differenceB = trackB - pointedTrackParameters.GetB();
  	  const double differenceQOverP = trackQOverP - pointedTrackParameters.GetQOverP();

  	  //Create merge candidate using calculated parameter differences and detector-dependent
  	  //weighting parameters.
  	  modutils::MergeCandidate mergeCandidate(pointedTrack,
  	  					  differenceX,differenceY,
  	  					  differenceA,differenceB,
  	  					  differenceQOverP,fDetectorSetup);

  	  // const double metric =
  	  //   sqrt(differenceX*differenceX/extrapolatedTrackParameters.GetCovarianceMatrix()[0][0] +
  	  // 	 differenceY*differenceY/extrapolatedTrackParameters.GetCovarianceMatrix()[1][1] +
  	  // 	 differenceA*differenceA/extrapolatedTrackParameters.GetCovarianceMatrix()[2][2] +
  	  // 	 differenceB*differenceB/extrapolatedTrackParameters.GetCovarianceMatrix()[3][3] );
	  
  	  //Add chi2/ndf of track.
  	  const double pointedChi2PerNDF =
	    pointedTrackParameters.GetChi2()/pointedTrackParameters.GetNdf();
  	  mergeCandidate.SetChi2PerNDF(pointedChi2PerNDF);
  	  //Add to list of potential tracks to merge.
  	  mergeCandidates.insert(mergeCandidate);
  	}
      } //end for(ZPlanes) 
      
      //If we have a non-empty list of merge candidates, see if the
      //best candidate is compatible with the original track.
      if (mergeCandidates.size() == 0)
  	continue;
      
      for (modutils::MergeCandidates::const_iterator candidateIt = mergeCandidates.begin(),
  	     candidateEnd = mergeCandidates.end();
  	   candidateIt != candidateEnd; ++candidateIt) {

	//Get the best compatible track from the list.
  	const modutils::LocalTrackSegment& pointedTrack = candidateIt->GetLocalTrackSegment();
  	const double metric = candidateIt->GetMetric();

  	//Loose-ish metric cut.
  	if (metric > fDetectorSetup.fMaxMergeMetric)
  	  continue;

	//Do a test refit. If the chi2 value is bad, don't perform the merge.
	modutils::KalmanTrackParameters trackParameters;
	std::vector<evt::Index<evt::rec::Cluster> > testClusters = localTrackSegment.GetClusters();
	for (modutils::ClusterVector::const_iterator clusterIt = pointedTrack.GetClusters().begin(),
	       clusterEnd = pointedTrack.GetClusters().end(); clusterIt != clusterEnd; ++clusterIt) {
	  testClusters.push_back(*clusterIt);
	}
	
	//Re-fit here.
	bool goodFit = false;
	const utl::Plane dummyPlane(utl::Point(0,0,0),utl::Vector(0,0,1));
	if (fDetectorSetup.fFitterType == modutils::eKalmanFilter) {
	  goodFit = fKalmanFilter.KalmanFit(testClusters,trackParameters);
	}
	if (fDetectorSetup.fFitterType == modutils::eStraightTrackFitter) {
	  modutils::StraightTrackFitResult fitResult;
	  goodFit = fStraightTrackFitter.FitStraightLine(testClusters, fitResult);
	  trackParameters.SetStraightFitParameters(fitResult,dummyPlane);
	}
	if (fDetectorSetup.fFitterType == modutils::eParabolicTrackFitter) {
	  modutils::ParabolicTrackFitResult fitResult;
	  goodFit = fParabolicTrackFitter.FitTrack(testClusters, fitResult);
	  trackParameters.SetParabolicFitParameters(fitResult,dummyPlane);
	}
	if (!goodFit)
	  continue;
	
	const double averageTrackResidualRefit = sqrt(trackParameters.GetAverageXResidual()*
						      trackParameters.GetAverageXResidual() +
						      trackParameters.GetAverageYResidual()*
						      trackParameters.GetAverageYResidual());
	
	//Make sure new average track residual is not too large.
	const double maxResidual = fDetectorSetup.fMaxAverageResidual;
	if (averageTrackResidualRefit > maxResidual)
	  continue;

	//Remove the tracks from both grids.
	proc::DetectionPlaneGrid& trackGrid =
	  procReconstruction.GetDetectionPlaneGrid(localTrackSegment.GetFirstZPlaneId());
        trackGrid.AddTrackSegmentToDisabledList(localTrackSegment.GetIndex());
        trackGrid.DisableTrackSegments();
	proc::DetectionPlaneGrid& pointedGrid =
	  procReconstruction.GetDetectionPlaneGrid(pointedTrack.GetFirstZPlaneId());
        pointedGrid.AddTrackSegmentToDisabledList(pointedTrack.GetIndex());
        pointedGrid.DisableTrackSegments();
	  
  	//Copy clusters into the parent track.
  	CopyClusters(localTrackSegment,pointedTrack);

  	//Sort clusters by z-position.
  	localTrackSegment.SortClusters();

  	//Get upstream grid Z.
  	const proc::DetectionPlaneGrid& firstGrid =
	  procReconstruction.GetDetectionPlaneGrid(localTrackSegment.GetFirstZPlaneId());
	
	const utl::Plane& firstZPlane =
	  (chamber.GetId() == det::TPCConst::eMTPCL &&
	   fClusterPositionsByIndex[localTrackSegment.GetFirstClusterIndex()].GetX() < 0) ?
	  fMTPCR->GetZPlane(firstGrid.GetZPlaneId()) :
	  chamber.GetZPlane(firstGrid.GetZPlaneId());
	
  	//We will merge these tracks.
  	segmentsToErase.insert(pointedTrack.GetIndex());
	
	//Final track refit.
	//Make sure track parameters are stored on z-plane.
  	if (fDetectorSetup.fFitterType == modutils::eKalmanFilter) {
	  fKalmanFilter.KalmanFit(localTrackSegment.GetClusters(), trackParameters);
	  fKalmanFilter.ExtrapolateToPlane(trackParameters,
					   trackParameters,
					   firstZPlane);
	  localTrackSegment.SetTrackParameters(trackParameters);
  	}
  	if (fDetectorSetup.fFitterType == modutils::eStraightTrackFitter) {
  	  modutils::StraightTrackFitResult fitResult;
  	  fStraightTrackFitter.FitStraightLine(localTrackSegment.GetClusters(), fitResult);
  	  trackParameters.SetStraightFitParameters(fitResult,firstZPlane);
  	  localTrackSegment.SetTrackParameters(trackParameters);
  	}
  	if (fDetectorSetup.fFitterType == modutils::eParabolicTrackFitter) {
  	  modutils::ParabolicTrackFitResult fitResult;
  	  fParabolicTrackFitter.FitTrack(localTrackSegment.GetClusters(), fitResult);
  	  trackParameters.SetParabolicFitParameters(fitResult,firstZPlane);
  	  localTrackSegment.SetTrackParameters(trackParameters);
  	}

	//Add to upstream grid.	
	proc::DetectionPlaneGrid& newGrid =
	  procReconstruction.GetDetectionPlaneGrid(localTrackSegment.GetFirstZPlaneId());
	newGrid.AddTrackSegmentToGrid(localTrackSegment.GetIndex(),trackParameters.GetPosition());
	
      } //end for(mergeCandidates)
    } //end for(trackSegments)

    //Remove merged segments.
    for (set<modutils::TrackSegmentIndex>::iterator indexIt = segmentsToErase.begin(),
    	   indexEnd = segmentsToErase.end(); indexIt != indexEnd; ++indexIt) {
      const modutils::LocalTrackSegments::iterator pointedTrackIt =
    	std::find_if(localTrackSegments.begin(), localTrackSegments.end(),
    		     modutils::FindLocalTrackByIndex(*indexIt));
      localTrackSegments.erase(pointedTrackIt);
    }
    return true;
  }

  bool TrackSegmentExtender::PerformClusterPickup(evt::proc::Reconstruction& procReconstruction,
  						  modutils::LocalTrackSegments& localTrackSegments,
  						  modutils::ClusterList& clustersOnTracks,
  						  const det::TPCChamber& chamber) {
    //Sort vector of local track segments before working with them.
    //This ensures that they'll be ordered by starting z-position.
    localTrackSegments.sort();
    
    //Start from track whose upstream end is furthest from the target.
    for (modutils::LocalTrackSegments::iterator localTrackSegmentIt = localTrackSegments.begin(),
  	   localTrackSegmentEnd = localTrackSegments.end(); localTrackSegmentIt != localTrackSegmentEnd;
  	 localTrackSegmentIt++) {
      
      //Get local track segment for modification.
      modutils::LocalTrackSegment& localTrackSegment = *localTrackSegmentIt;

      //Residual cuts. Tighter here.
      const double maxResidualX = 3*localTrackSegment.GetTrackParameters().GetAverageXResidual();
      const double maxResidualY = 3*localTrackSegment.GetTrackParameters().GetAverageYResidual();
      
      const modutils::ZPlaneId trackFirstZPlane = localTrackSegment.GetFirstZPlaneId();
	const utl::Plane& firstZPlane =
	  (chamber.GetId() == det::TPCConst::eMTPCL &&
	   fClusterPositionsByIndex[localTrackSegment.GetFirstClusterIndex()].GetX() < 0) ?
	  fMTPCR->GetZPlane(trackFirstZPlane) :
	  chamber.GetZPlane(trackFirstZPlane);
	

      //Flag for re-sorting clusters if we pick any up. Sorting and re-fitting performed after
      //picking up all compatible clusters.
      bool pickedUpClusters = false;

      //Loop through all planes in detector.
      for (unsigned int nextZPlane = fDetectorSetup.fZPlaneIdRange.GetBegin(),
  	     lastZPlane = fDetectorSetup.fZPlaneIdRange.GetEnd();
  	   nextZPlane < lastZPlane; ++nextZPlane ) {
	
	const proc::DetectionPlaneGrid& nextLayerGrid =
	  procReconstruction.GetDetectionPlaneGrid(nextZPlane);

  	//Extrapolate track to next detection plane and calculate allowed (x,y) window for clusters to be
  	//in for merging.
  	const double nextGridPlaneZ = nextLayerGrid.GetGridZ();
  	modutils::KalmanTrackParameters extrapolatedTrackParameters;

	//Test for MTPCL vs MTPCR.
	if (fDetectorSetup.fFitterType == modutils::eKalmanFilter) {
	  fKalmanFilter.ExtrapolateToGlobalZ(localTrackSegment.GetTrackParameters(),
					     extrapolatedTrackParameters,
					     nextGridPlaneZ);
  	}
  	if (fDetectorSetup.fFitterType == modutils::eStraightTrackFitter) { 
 	  extrapolatedTrackParameters = localTrackSegment.GetTrackParameters();
  	  extrapolatedTrackParameters.ExtrapolateStraightFitToZ(nextGridPlaneZ);
  	}	
  	if (fDetectorSetup.fFitterType == modutils::eParabolicTrackFitter) {
  	  extrapolatedTrackParameters = localTrackSegment.GetTrackParameters();
  	  extrapolatedTrackParameters.ExtrapolateParabolicFitToZ(nextGridPlaneZ);
  	}	

	//Get plane object corresponding to z-plane. If TPC Id is MTPCL and
	//cluster X is negative, the cluster is in MTPCR.
	const utl::Plane& zPlane = (chamber.GetId() == det::TPCConst::eMTPCL &&
				    extrapolatedTrackParameters.GetPosition().GetX() < 0) ?
	  fMTPCR->GetZPlane(nextZPlane) :
	  chamber.GetZPlane(nextZPlane);
  	if (fDetectorSetup.fFitterType == modutils::eKalmanFilter) {
  	  fKalmanFilter.ExtrapolateToPlane(localTrackSegment.GetTrackParameters(),
					   extrapolatedTrackParameters,
					   zPlane);
  	}
  	if (fDetectorSetup.fFitterType == modutils::eStraightTrackFitter) {
  	  extrapolatedTrackParameters = localTrackSegment.GetTrackParameters();
  	  extrapolatedTrackParameters.ExtrapolateStraightFitToPlane(zPlane);
  	}	
  	if (fDetectorSetup.fFitterType == modutils::eParabolicTrackFitter) {
  	  extrapolatedTrackParameters = localTrackSegment.GetTrackParameters();
  	  extrapolatedTrackParameters.ExtrapolateParabolicFitToPlane(zPlane);
  	}	
	
  	const utl::Point& extrapolatedPosition = extrapolatedTrackParameters.GetPosition();
	
	//Make sure we haven't left the chamber. MTPCL requires a check of MTPCR, since they
	//are treated at the same time.
	if (chamber.GetId() == det::TPCConst::eMTPCL) {
	  if (!chamber.IsIn(extrapolatedPosition) &&
	      !fMTPCR->IsIn(extrapolatedPosition) ) 
	    continue;
	}
	else {
	  if (!chamber.IsIn(extrapolatedPosition)) 
	    continue;
	}
	
	//Perform lone cluster pickup. Only do this for a few padrows -- extrapolation from
	//badly measured tracks should not span the entire chamber.
	//Get subset of cluster indices that are within acceptable distance from our central cluster.
	const proc::ClusterVector& clusterIndicesNextDetectionLayer =
	  nextLayerGrid.GetClustersFromPositionAndRange(extrapolatedPosition.GetX(),
							extrapolatedPosition.GetY(),
							maxResidualX,maxResidualY);
	
	//Loop through clusters on next detection plane.
	for ( proc::tpc::Clusters::ClusterIndices::const_iterator nextLayerIt =
		clusterIndicesNextDetectionLayer.begin(), nextLayerEnd =
		clusterIndicesNextDetectionLayer.end(); nextLayerIt != nextLayerEnd; ++nextLayerIt ) {

	  const utl::Point& clusterPosition = fClusterPositionsByIndex[*nextLayerIt];	  
	  const modutils::ZPlaneId& zPlane = fClusterZPlanesByIndex[*nextLayerIt];	  

	  //Don't add cluster if it's already on a track segment.
	  modutils::ClusterList::iterator clusterIndexIt = clustersOnTracks.find(*nextLayerIt);
	  const bool isUsed = (clusterIndexIt != clustersOnTracks.end()) ? 
	    true : false;
	  if (isUsed)
	    continue;
	  
	  const double zOffset = abs(extrapolatedPosition.GetZ() - clusterPosition.GetZ());

	  if (zOffset > 1*utl::cm) {
	    
	    if (fDetectorSetup.fFitterType == modutils::eKalmanFilter) {
	      fKalmanFilter.ExtrapolateToPlane(localTrackSegment.GetTrackParameters(),
					       extrapolatedTrackParameters,
					       firstZPlane);
	    }
	    if (fDetectorSetup.fFitterType == modutils::eStraightTrackFitter) { 
	      extrapolatedTrackParameters = localTrackSegment.GetTrackParameters();
	      extrapolatedTrackParameters.ExtrapolateStraightFitToPlane(firstZPlane);
	    }	
	    if (fDetectorSetup.fFitterType == modutils::eParabolicTrackFitter) {
	      extrapolatedTrackParameters = localTrackSegment.GetTrackParameters();
	      extrapolatedTrackParameters.ExtrapolateParabolicFitToPlane(firstZPlane);
	    }	
	  }
	  
	  //Now compare extrapolated track position at cluster Z to cluster position.
	  const double residualX = fabs(extrapolatedPosition.GetX() - clusterPosition.GetX());
	  const double residualY = fabs(extrapolatedPosition.GetY() - clusterPosition.GetY());

	  //If cluster position is no good, move on.
	  if (residualX > maxResidualX ||
	      residualY > maxResidualY )
	    continue;
	  
	  //Do a test refit. If the chi2 value is bad, don't perform the merge.
	  modutils::KalmanTrackParameters trackParameters;
	  std::vector<evt::Index<evt::rec::Cluster> > testClusters = localTrackSegment.GetClusters();
	  testClusters.push_back(*nextLayerIt);
	  
	  //Re-fit here.
	  bool goodFit = false;
	  if (fDetectorSetup.fFitterType == modutils::eKalmanFilter) {
	    goodFit = fKalmanFilter.KalmanFit(testClusters,trackParameters);
	  }
	  if (fDetectorSetup.fFitterType == modutils::eStraightTrackFitter) {
	    modutils::StraightTrackFitResult fitResult;
	    goodFit = fStraightTrackFitter.FitStraightLine(testClusters, fitResult);
	    trackParameters.SetStraightFitParameters(fitResult,firstZPlane);
	  }
	  if (fDetectorSetup.fFitterType == modutils::eParabolicTrackFitter) {
	    modutils::ParabolicTrackFitResult fitResult;
	    goodFit = fParabolicTrackFitter.FitTrack(testClusters, fitResult);
	    trackParameters.SetParabolicFitParameters(fitResult,firstZPlane);
	  }
	  if (!goodFit)
	    continue;
	  const double averageTrackResidualRefit = sqrt(trackParameters.GetAverageXResidual()*
							trackParameters.GetAverageXResidual() +
							trackParameters.GetAverageYResidual()*
							trackParameters.GetAverageYResidual());
	  
	  //Make sure re-fit did not return a bad track.
	  const double maxResidual = fDetectorSetup.fMaxAverageResidual;
	  if (averageTrackResidualRefit > maxResidual)
	    continue;
	  
	  clustersOnTracks.insert(*nextLayerIt);
	  localTrackSegment.AddCluster(*nextLayerIt,clusterPosition,zPlane);
	  
	  pickedUpClusters = true;
	} //End clusters on next layer loop.
      } //End z-plane loop.
      
      //If we've picked up any clusters, re-calculate track parameters and place on the correct
      //detection plane grid.
      if (pickedUpClusters) {

	//Sort clusters and get information about new first/last clusters.
	localTrackSegment.SortClusters();
	const proc::ClusterVector& clusters = localTrackSegment.GetClusters();
	// const utl::Point& firstClusterPosition =
	//   fClusterPositionsByIndex[localTrackSegment.GetFirstClusterIndex()];
	const modutils::ZPlaneId firstZPlaneId = localTrackSegment.GetFirstZPlaneId();

	const utl::Plane& firstZPlane =
	  (chamber.GetId() == det::TPCConst::eMTPCL &&
	   fClusterPositionsByIndex[localTrackSegment.GetFirstClusterIndex()].GetX() < 0) ?
	  fMTPCR->GetZPlane(firstZPlaneId) :
	  chamber.GetZPlane(firstZPlaneId);

	//Remove the track from its grid.
	proc::DetectionPlaneGrid& trackGrid =
	  procReconstruction.GetDetectionPlaneGrid(trackFirstZPlane);
        trackGrid.AddTrackSegmentToDisabledList(localTrackSegment.GetIndex());
        trackGrid.DisableTrackSegments();
	
	//Final track refit.
	//Make sure track parameters are on z-plane.
	modutils::KalmanTrackParameters trackParameters;
	if (fDetectorSetup.fFitterType == modutils::eKalmanFilter) {
	  fKalmanFilter.KalmanFit(clusters, trackParameters);
	  fKalmanFilter.ExtrapolateToPlane(trackParameters,
					   trackParameters,
					   firstZPlane);
	  localTrackSegment.SetTrackParameters(trackParameters);
	}
	if (fDetectorSetup.fFitterType == modutils::eStraightTrackFitter) {
	  modutils::StraightTrackFitResult fitResult;
	  fStraightTrackFitter.FitStraightLine(clusters, fitResult);
	  trackParameters.SetStraightFitParameters(fitResult,firstZPlane);
	  localTrackSegment.SetTrackParameters(trackParameters);
	}
	if (fDetectorSetup.fFitterType == modutils::eParabolicTrackFitter) {
	  modutils::ParabolicTrackFitResult fitResult;
	  fParabolicTrackFitter.FitTrack(clusters, fitResult);
	  trackParameters.SetParabolicFitParameters(fitResult,firstZPlane);
	  localTrackSegment.SetTrackParameters(trackParameters);
	}

	//Now add track to new grid with re-fit position.	
	proc::DetectionPlaneGrid& newGrid =
	  procReconstruction.GetDetectionPlaneGrid(firstZPlaneId);
	newGrid.AddTrackSegmentToGrid(localTrackSegment.GetIndex(),trackParameters.GetPosition());
      } //End if (pickedUpClusters).
    } //End track segments loop.
    // INFO("That's all, folks!!");

    return true;
  }
  
  //Copy track clusters/cluster info. For track merging. 
  void TrackSegmentExtender::CopyClusters(modutils::LocalTrackSegment& trackToKeep,
					  const modutils::LocalTrackSegment& trackToRemove) {
    
    //Add clusters to keeper track. Sort when finished.
    modutils::ClusterVector trackNewClusters = trackToRemove.GetClusters();
    for (modutils::ClusterVector::const_iterator clusterIt = trackNewClusters.begin(), 
	   clusterEnd = trackNewClusters.end(); clusterIt != clusterEnd; ++clusterIt) {
      trackToKeep.AddCluster(*clusterIt,fClusterPositionsByIndex[*clusterIt],
			     fClusterZPlanesByIndex[*clusterIt]);
    }
    trackToKeep.SortClusters();
  }
  
}
