#include "KalmanFilterWB.h"
#include <det/Detector.h>
#include <det/TPC.h>
#include <det/MagneticField.h>
#include <fwk/CentralConfig.h>
#include <utl/Branch.h>
#include <utl/ErrorLogger.h>
#include <utl/GeometryUtilities.h>
#include <utl/ShineUnits.h>
#include <utl/PhysicalConst.h>
#include <utl/TrackParametersConversion.h>

#include <sstream>
#include <iomanip>

using namespace fwk;
using namespace utl;
using namespace std;

using namespace evt;
using namespace evt::rec;
using namespace det;

namespace modutils {

  void KalmanFilter::SetTrackParameters(TrackFit& track,
					modutils::KalmanTrackParameters& trackParameters) {
    
    utl::CovarianceMatrix cov(5);
    
    //Covariance matrix translation.
    using namespace modutils;
    cov[KalmanTrackParameters::eX][KalmanTrackParameters::eX     ] = track.C.C00;
    cov[KalmanTrackParameters::eX][KalmanTrackParameters::eY     ] = track.C.C10;
    cov[KalmanTrackParameters::eX][KalmanTrackParameters::eA     ] = track.C.C20;
    cov[KalmanTrackParameters::eX][KalmanTrackParameters::eB     ] = track.C.C30;
    cov[KalmanTrackParameters::eX][KalmanTrackParameters::eQOverP] = track.C.C40;
    
    cov[KalmanTrackParameters::eY][KalmanTrackParameters::eY     ] = track.C.C11;
    cov[KalmanTrackParameters::eY][KalmanTrackParameters::eA     ] = track.C.C21;
    cov[KalmanTrackParameters::eY][KalmanTrackParameters::eB     ] = track.C.C31;
    cov[KalmanTrackParameters::eY][KalmanTrackParameters::eQOverP] = track.C.C41;

    cov[KalmanTrackParameters::eA][KalmanTrackParameters::eA     ] = track.C.C22;
    cov[KalmanTrackParameters::eA][KalmanTrackParameters::eB     ] = track.C.C32;
    cov[KalmanTrackParameters::eA][KalmanTrackParameters::eQOverP] = track.C.C42;
    
    cov[KalmanTrackParameters::eB][KalmanTrackParameters::eB     ] = track.C.C33;
    cov[KalmanTrackParameters::eB][KalmanTrackParameters::eQOverP] = track.C.C43;
    
    cov[KalmanTrackParameters::eQOverP][KalmanTrackParameters::eQOverP] = track.C.C44;
    
    trackParameters.SetCovarianceMatrix(cov);
    
    //Track parameters.
    trackParameters.SetX(track.T[0]);
    trackParameters.SetY(track.T[1]);
    trackParameters.SetA(track.T[2]);
    trackParameters.SetB(track.T[3]);
    trackParameters.SetQOverP(track.T[4]);

    //Auxilliary parameters and chi2/ndf.
    trackParameters.SetZ(track.T[5]);
    trackParameters.SetTrackLength(track.T[6]);
    trackParameters.SetChi2(track.Chi2);
    trackParameters.SetNdf(track.NDF);
    trackParameters.SetAverageXResidual(track.fXResidualAverage);
    trackParameters.SetAverageYResidual(track.fYResidualAverage);

    //Trimmed clusters from chi2 calculation.
    trackParameters.SetTrimmedClusters(track.fTrimmedClusters);
  }

  bool KalmanFilter::KalmanFit(const ClusterVector & clusterIndices,
			       modutils::KalmanTrackParameters& trackParameters)
  {

    //Make sure RecEvent pointer is set.
    if (fRecEvent == 0) {
      FATAL("Use Init() function first to initialize the event pointer!\n");
      return false;
    }
    const evt::RecEvent& recEvent = *fRecEvent;
    
    if(clusterIndices.size() < 2) {
      WARNING("Number of clusters < 2. Do nothing!");
      return false;
    }

    //Sort clusters by z-position here. Later we can move this elsewhere if it
    //proves too time-consuming (like in local tracking), and make this function
    //expect a sorted list of clusters. This would prevent unnecessary re-sorting
    //when re-fitting tracks.
    ClusterVector sortedClusterIndices = SortClustersByZ(clusterIndices);

    //Initialize vector of measurements.
    vector<Measurement> measurements;

    //Get magnetic field interface.
    const det::MagneticField& magneticField = det::Detector::GetInstance().GetMagneticField(); 
    
    //Add Kalman Fitler measurement objects.
    for (ClusterVector::const_iterator it = sortedClusterIndices.begin(),
	   itEnd = sortedClusterIndices.end(); it != itEnd; ++it) {

      const Cluster& cluster = recEvent.Get(*it);
      Measurement measurement;
      measurement.fX = cluster.GetPosition().GetX();
      measurement.fY = cluster.GetPosition().GetY();
      measurement.fZ = cluster.GetPosition().GetZ();
      measurement.fClusterIndex = cluster.GetIndex();
      measurement.fWeight = 1;

      //Get y-part of magnetic field value. This is for the Kalman Filter
      //guess step -- it needs the integral of By over z.
      const utl::Vector& fieldValue = magneticField.GetField(cluster.GetPosition());
      measurement.fField.Y = fieldValue.GetY()/kilogauss;

      using namespace ClusterConst;
      const double clusterUncertaintyX = cluster.GetPositionUncertainty(eX);
      const double clusterUncertaintyY = cluster.GetPositionUncertainty(eY);

      measurement.fSigmaX = clusterUncertaintyX;
      measurement.fSigma2X = clusterUncertaintyX*clusterUncertaintyX;
      measurement.fSigmaY = clusterUncertaintyY;
      measurement.fSigma2Y = clusterUncertaintyY*clusterUncertaintyY; 
      measurement.fSigma2XY = fabs(cluster.GetShapeCovariance(eX,eY));

      measurements.push_back(measurement);
    } // end of loop over clusters

    //Field integral for initial momentum guess.
    double previousZ = measurements.back().fZ;
    double Sy = 0;
    double sy = 0;
    for (vector<Measurement>::reverse_iterator it = measurements.rbegin(),
	   itEnd = measurements.rend(); it != itEnd; ++it) {
      Measurement& measurement = *it;
      double dz = measurement.fZ-previousZ;
      double Hy = measurement.fField.Y;
      Sy += dz*sy + dz*dz*Hy/2.f;
      sy += dz*Hy;
      measurement.fSy = Sy;
      previousZ = measurement.fZ;
    }
    
    //Track fit object.
    TrackFit track;

    //Perform fit.    
    FitTrack(track,
	     measurements,
	     eForward,
	     false);
    
    //Store SHINE track parameters.
    SetTrackParameters(track, trackParameters);

    //Find out if track parameter 0 (x-position) is NaN. Return false if it is.
    return (trackParameters.GetX() == trackParameters.GetX() &&
	    trackParameters.GetMomentum().GetX() == trackParameters.GetMomentum().GetX());
  }
  
  bool KalmanFilter::KalmanFitWithVertex(const ClusterVector & clusterIndices,
  					 const evt::rec::Vertex& vertex,
  					 modutils::KalmanTrackParameters& seedParameters)
  {
    
    //Make sure RecEvent pointer is set.
    if (fRecEvent == 0) {
      FATAL("Use Init() function first to initialize the event pointer!\n");
      return false;
    }
    const evt::RecEvent& recEvent = *fRecEvent;
    
    if(clusterIndices.size() < 2) {
      WARNING("Number of clusters < 2. Do nothing!");
      return false;
    }
    
    //Initialize vector of measurements.
    vector<Measurement> measurements;
    
    //Get magnetic field interface.
    const det::MagneticField& magneticField = det::Detector::GetInstance().GetMagneticField(); 

    //First, add vertex as a measurement.
    Measurement vertexMeasurement;
    vertexMeasurement.fX = vertex.GetPosition().GetX();
    vertexMeasurement.fY = vertex.GetPosition().GetY();
    vertexMeasurement.fZ = vertex.GetPosition().GetZ();
    vertexMeasurement.fWeight = 1;
    //Get y-part of magnetic field value. This is for the Kalman Filter
    //guess step -- it needs the integral of By over z.
    const utl::Vector& fieldValue = magneticField.GetField(vertex.GetPosition());
    vertexMeasurement.fField.Y = fieldValue.GetY()/kilogauss;
    
    //Get uncertanties.
    const double vertexUncertaintyX = sqrt(vertex.GetCovarianceMatrix()[0][0]);
    const double vertexUncertaintyY = sqrt(vertex.GetCovarianceMatrix()[1][1]);
    if (vertexUncertaintyX == 0 || vertexUncertaintyY == 0) {
      WARNING("Track refit with vertex failed. X or Y position uncertainty was set to 0.");
      return false;
    }
    vertexMeasurement.fSigmaX = vertexUncertaintyX;
    vertexMeasurement.fSigma2X = vertexUncertaintyX*vertexUncertaintyX;
    vertexMeasurement.fSigmaY = vertexUncertaintyY;
    vertexMeasurement.fSigma2Y = vertexUncertaintyY*vertexUncertaintyY;
    vertexMeasurement.fSigma2XY = vertexMeasurement.fSigma2X + vertexMeasurement.fSigma2Y;

    //Add to vector.
    measurements.push_back(vertexMeasurement);
    
    //Sort clusters by z-position here. Later we can move this elsewhere if it
    //proves too time-consuming (like in local tracking), and make this function
    //expect a sorted list of clusters. This would prevent unnecessary re-sorting
    //when re-fitting tracks.
    ClusterVector sortedClusterIndices = SortClustersByZ(clusterIndices);

    //Add Kalman Fitler measurement objects using clusters.
    for (ClusterVector::const_iterator it = sortedClusterIndices.begin(),
	   itEnd = sortedClusterIndices.end(); it != itEnd; ++it) {

      const Cluster& cluster = recEvent.Get(*it);
      Measurement measurement;
      measurement.fX = cluster.GetPosition().GetX();
      measurement.fY = cluster.GetPosition().GetY();
      measurement.fZ = cluster.GetPosition().GetZ();
      measurement.fClusterIndex = cluster.GetIndex();
      measurement.fWeight = 1;
      //Get y-part of magnetic field value. This is for the Kalman Filter
      //guess step -- it needs the integral of By over z.
      const utl::Vector& fieldValue = magneticField.GetField(cluster.GetPosition());
      measurement.fField.Y = fieldValue.GetY()/kilogauss;

      using namespace ClusterConst;
      const double clusterUncertaintyX = cluster.GetPositionUncertainty(eX);
      const double clusterUncertaintyY = cluster.GetPositionUncertainty(eY);

      measurement.fSigmaX = clusterUncertaintyX;
      measurement.fSigma2X = clusterUncertaintyX*clusterUncertaintyX;
      measurement.fSigmaY = clusterUncertaintyY;
      measurement.fSigma2Y = clusterUncertaintyY*clusterUncertaintyY;
      measurement.fSigma2XY = fabs(cluster.GetShapeCovariance(eX,eY));

      //Add to vector.
      measurements.push_back(measurement);
    } // end of loop over clusters

    //Field integral for initial momentum guess.
    double previousZ = measurements.back().fZ;
    double Sy = 0;
    double sy = 0;
    for (vector<Measurement>::reverse_iterator it = measurements.rbegin(),
	   itEnd = measurements.rend(); it != itEnd; ++it) {
      Measurement& measurement = *it;
      double dz = measurement.fZ-previousZ;
      double Hy = measurement.fField.Y;
      Sy += dz*sy + dz*dz*Hy/2.f;
      sy += dz*Hy;
      measurement.fSy = Sy;
      previousZ = measurement.fZ;
    }
    
    //Track fit object. Pass previously-calculated parameters as a seed for this fit.
    TrackFit track;
    track.T[0] = seedParameters.GetX();
    track.T[1] = seedParameters.GetY();
    track.T[2] = seedParameters.GetA();
    track.T[3] = seedParameters.GetB();
    track.T[4] = seedParameters.GetQOverP();
    track.T[5] = seedParameters.GetZ();
    //Extrapolate to end of track.
    ExtrapolateToGlobalZ(track.T,track.C,measurements.back().fZ);
	
    //Perform fit.    
    FitTrack(track,
	     measurements,
	     eForward,
	     true);

    //Store SHINE track parameters.
    SetTrackParameters(track, seedParameters);

    //Find out if track parameter 0 (x-position) is NaN. Return false if it is.
    return (seedParameters.GetX() == seedParameters.GetX() &&
	    seedParameters.GetMomentum().GetX() == seedParameters.GetMomentum().GetX());
  }
  
  void
  KalmanFilter::ExtrapolateToGlobalZ(double T[],
				     Cov& C,
				     const double targetZ)
  {
    //Get magnetic field interface.
    const det::MagneticField& magneticField = det::Detector::GetInstance().GetMagneticField();
    //Create containers for transfering magnetic field information to Extrapolation function.
    FieldRegion f;
    FieldVector H0, H1, H2;

    //Begin with maximum step size. Adjust later as necessary.
    double stepSize = fMaximumStepSize;

    //If we're close enough to the target Z, just use the difference in positions.
    if (fabs(T[5] - targetZ) < stepSize) stepSize = fabs(T[5] - targetZ);

    //Step in the other direction if we need to.
    if(targetZ < T[5]) stepSize *= -1;

    //Return if we're already close enough.
    if (fabs(targetZ - T[5]) < fExtrapolationTolerance) return;

    //Define positions to sample the magnetic field. nextZ is one full stepSize
    //away, and intermediateZ is half of this distance. Calculate approximate
    //x and y at these points using track parameters at the first point.
    double currentX, currentY, currentZ,
      intermediateX, intermediateY, intermediateZ,
      nextX, nextY, nextZ, linearizedQOverP;
    currentX      = T[0];
    currentY      = T[1];
    currentZ      = T[5];
    nextX         = currentX + T[2]*stepSize;
    nextY         = currentY + T[3]*stepSize;
    intermediateX = currentX + (nextX - currentX)/2.;
    intermediateY = currentY + (nextY - currentY)/2.;

    while(fabs(currentZ - targetZ) > fExtrapolationTolerance) {

      if(fabs(currentZ - targetZ) > fabs(stepSize)) nextZ = currentZ + stepSize;
      else{
	stepSize =  targetZ - currentZ;
	nextZ = targetZ;
      }

      linearizedQOverP = T[4];

      nextX         = currentX + T[2]*stepSize;
      nextY         = currentY + T[3]*stepSize;
      intermediateX = currentX + (nextX - currentX)/2.;
      intermediateY = currentY + (nextY - currentY)/2.;

      intermediateZ = (currentZ + nextZ)/2;

      const utl::Vector& currentField =
	magneticField.GetField(utl::Point(currentX,currentY,currentZ));
      H0.X = currentField.GetX()/utl::kilogauss;
      H0.Y = currentField.GetY()/utl::kilogauss;
      H0.Z = currentField.GetZ()/utl::kilogauss;
      
      const utl::Vector& nextField =
	magneticField.GetField(utl::Point(nextX,nextY,nextZ));
      H2.X = nextField.GetX()/utl::kilogauss;
      H2.Y = nextField.GetY()/utl::kilogauss;
      H2.Z = nextField.GetZ()/utl::kilogauss;

      const utl::Vector& intermediateField =
	magneticField.GetField(utl::Point(intermediateX,intermediateY,intermediateZ));
      H1.X = intermediateField.GetX()/utl::kilogauss;
      H1.Y = intermediateField.GetY()/utl::kilogauss;
      H1.Z = intermediateField.GetZ()/utl::kilogauss;

      //FIXME: THINK ABOUT THIS!!! VERY IMPORTANT!!
      const double maxGradient = 0.1;

      const double gradient = (H2.Y - H0.Y)/(nextZ - currentZ);

      if (abs(gradient) > maxGradient && stepSize >= fMinimumStepSize) {
	stepSize /= 2;
	continue;
      }
      else {
	stepSize *= 2;
	if (abs(stepSize) > fMaximumStepSize) {
	  stepSize = Sign(stepSize)*fMaximumStepSize;
	}
	if (abs(stepSize) < fMinimumStepSize) {
	  stepSize = Sign(stepSize)*fMinimumStepSize;
	}
      }

      
      // f.Set( H0, currentZ, H1, intermediateZ, H2, nextZ);
      // f.Set( H2, nextZ, H1, intermediateZ, H0, currentZ);

      // //Best one so far!
      // f.Set( H2, nextZ, H0, currentZ, H1, intermediateZ );

      //Test.
      f.Set( H2, nextZ, H1, intermediateZ, H0, currentZ);
      
      ExtrapolateALight(T, C, nextZ, linearizedQOverP, f);
      
      currentX      = nextX;
      currentY      = nextY;
      currentZ      = nextZ;
    }
  }

  //Extrapolate track parameters using KalmanTrackParameters.
  bool KalmanFilter::ExtrapolateToGlobalZ(const modutils::KalmanTrackParameters& trackParametersIn,
					  modutils::KalmanTrackParameters& trackParametersOut,
					  const double z)
  {
    //Copy track parameters and covariance matrix into a structure that can be used by Extrapolate function.
    double trackParameterBuffer[6];
    trackParameterBuffer[0] = trackParametersIn.GetX();
    trackParameterBuffer[1] = trackParametersIn.GetY();
    trackParameterBuffer[5] = trackParametersIn.GetZ();
    trackParameterBuffer[2] = trackParametersIn.GetA();
    trackParameterBuffer[3] = trackParametersIn.GetB();
    trackParameterBuffer[4] = trackParametersIn.GetQOverP();

    Cov covarianceMatrixBuffer;
    const utl::CovarianceMatrix& cov = trackParametersIn.GetCovarianceMatrix();

    covarianceMatrixBuffer.C00 = cov[KalmanTrackParameters::eX][KalmanTrackParameters::eX     ];
    covarianceMatrixBuffer.C10 = cov[KalmanTrackParameters::eX][KalmanTrackParameters::eY     ];
    covarianceMatrixBuffer.C20 = cov[KalmanTrackParameters::eX][KalmanTrackParameters::eA     ];
    covarianceMatrixBuffer.C30 = cov[KalmanTrackParameters::eX][KalmanTrackParameters::eB     ];
    covarianceMatrixBuffer.C40 = cov[KalmanTrackParameters::eX][KalmanTrackParameters::eQOverP];

    covarianceMatrixBuffer.C11 = cov[KalmanTrackParameters::eY][KalmanTrackParameters::eY     ];
    covarianceMatrixBuffer.C21 = cov[KalmanTrackParameters::eY][KalmanTrackParameters::eA     ];
    covarianceMatrixBuffer.C31 = cov[KalmanTrackParameters::eY][KalmanTrackParameters::eB     ];
    covarianceMatrixBuffer.C41 = cov[KalmanTrackParameters::eY][KalmanTrackParameters::eQOverP];

    covarianceMatrixBuffer.C22 = cov[KalmanTrackParameters::eA][KalmanTrackParameters::eA     ];
    covarianceMatrixBuffer.C32 = cov[KalmanTrackParameters::eA][KalmanTrackParameters::eB     ];
    covarianceMatrixBuffer.C42 = cov[KalmanTrackParameters::eA][KalmanTrackParameters::eQOverP];

    covarianceMatrixBuffer.C33 = cov[KalmanTrackParameters::eB][KalmanTrackParameters::eB     ];
    covarianceMatrixBuffer.C43 = cov[KalmanTrackParameters::eB][KalmanTrackParameters::eQOverP];

    covarianceMatrixBuffer.C44 = cov[KalmanTrackParameters::eQOverP][KalmanTrackParameters::eQOverP];
    
    //Perform the extrapolation!
    ExtrapolateToGlobalZ(trackParameterBuffer, covarianceMatrixBuffer, z);

    //Translate results to SHINE structures.
    utl::CovarianceMatrix covOut(5);
    
    trackParametersOut.SetX(trackParameterBuffer[0]);
    trackParametersOut.SetY(trackParameterBuffer[1]);
    trackParametersOut.SetZ(z);
    trackParametersOut.SetA(trackParameterBuffer[2]);
    trackParametersOut.SetB(trackParameterBuffer[3]);
    trackParametersOut.SetQOverP(trackParameterBuffer[4]);

    //Covariance matrix filling.
    using namespace modutils;
    covOut[KalmanTrackParameters::eX][KalmanTrackParameters::eX     ] = covarianceMatrixBuffer.C00;
    covOut[KalmanTrackParameters::eX][KalmanTrackParameters::eY     ] = covarianceMatrixBuffer.C10;
    covOut[KalmanTrackParameters::eX][KalmanTrackParameters::eA     ] = covarianceMatrixBuffer.C20;
    covOut[KalmanTrackParameters::eX][KalmanTrackParameters::eB     ] = covarianceMatrixBuffer.C30;
    covOut[KalmanTrackParameters::eX][KalmanTrackParameters::eQOverP] = covarianceMatrixBuffer.C40;
    
    covOut[KalmanTrackParameters::eY][KalmanTrackParameters::eY     ] = covarianceMatrixBuffer.C11;
    covOut[KalmanTrackParameters::eY][KalmanTrackParameters::eA     ] = covarianceMatrixBuffer.C21;
    covOut[KalmanTrackParameters::eY][KalmanTrackParameters::eB     ] = covarianceMatrixBuffer.C31;
    covOut[KalmanTrackParameters::eY][KalmanTrackParameters::eQOverP] = covarianceMatrixBuffer.C41;

    covOut[KalmanTrackParameters::eA][KalmanTrackParameters::eA     ] = covarianceMatrixBuffer.C22;
    covOut[KalmanTrackParameters::eA][KalmanTrackParameters::eB     ] = covarianceMatrixBuffer.C32;
    covOut[KalmanTrackParameters::eA][KalmanTrackParameters::eQOverP] = covarianceMatrixBuffer.C42;
    
    covOut[KalmanTrackParameters::eB][KalmanTrackParameters::eB     ] = covarianceMatrixBuffer.C33;
    covOut[KalmanTrackParameters::eB][KalmanTrackParameters::eQOverP] = covarianceMatrixBuffer.C43;
    
    covOut[KalmanTrackParameters::eQOverP][KalmanTrackParameters::eQOverP] = covarianceMatrixBuffer.C44;

    trackParametersOut.SetCovarianceMatrix(covOut);

    //For now, just copy track's fit chi2 value. Think about this some more!!
    trackParametersOut.SetChi2(trackParametersIn.GetChi2());
    trackParametersOut.SetNdf(trackParametersIn.GetNdf());

    //Find out if track parameter 0 (x-position) is NaN. Return false if it is.
    return (trackParameterBuffer[0] == trackParameterBuffer[0] &&
	    trackParameterBuffer[2] == trackParameterBuffer[2]);
  }


  //Extrapolate track parameters to a plane.
  bool KalmanFilter::ExtrapolateToPlane(const modutils::KalmanTrackParameters& trackParametersIn,
					modutils::KalmanTrackParameters& trackParametersOut,
					const utl::Plane& plane,
					const double tolerance)
  {
    unsigned int nAttempts = 0;
    const unsigned int maxAttempts = 10;
    //Copy track parameters.
    trackParametersOut = trackParametersIn;
    //Iterate until we are close enough to the desired plane.
    while (fabs(GeometryUtilities::Distance(trackParametersOut.GetPosition(),plane)) > tolerance &&
	   nAttempts < maxAttempts) {
      const utl::Point& position = trackParametersOut.GetPosition();
      const utl::Line trajectory(position,trackParametersOut.GetMomentum());
      //Extrapolate towards plane using straight line.
      const utl::Point& intersection = GeometryUtilities::Intersection(plane,trajectory);
      ExtrapolateToGlobalZ(trackParametersOut,trackParametersOut,intersection.GetZ());
      ++nAttempts;
    }
    if (nAttempts == maxAttempts)
      return false;
    return true;
  }
  
  //Function for sorting clusters by z-position (ascending).
  ClusterVector KalmanFilter::SortClustersByZ(const ClusterVector& clusterIndices) {
    
    //Map for storage and sorting. Key value is z-position.
    std::map<double, evt::Index<evt::rec::Cluster> > orderedClusterIndices;
    //Fill map.
    for(ClusterVector::const_iterator clusterIt = clusterIndices.begin();
	clusterIt!= clusterIndices.end(); ++clusterIt) {
      const double clusterZ = fRecEvent->Get(*clusterIt).GetPosition().GetZ();
      orderedClusterIndices[clusterZ] = *clusterIt;
    }
    //Fill new vector with sorted indices.
    ClusterVector sortedClusterIndices;
    for (std::map<double, ClusterIndex>::const_iterator it = orderedClusterIndices.begin(),
	   itEnd = orderedClusterIndices.end(); it != itEnd; ++it) {
      sortedClusterIndices.push_back(it->second);
    }
    return sortedClusterIndices;
  }

} // end of namespace modutils
