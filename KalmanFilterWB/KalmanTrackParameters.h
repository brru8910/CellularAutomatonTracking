#ifndef _utl_KalmanTrackParameters_h_
#define _utl_KalmanTrackParameters_h_

#include <utl/CovarianceMatrix.h>
#include <evt/rec/RecEventConst.h>
#include <evt/rec/Track.h>
#include <evt/rec/Vertex.h>
#include <evt/rec/VertexTrack.h>
#include <utl/ErrorLogger.h>
#include <utl/Fitted2DLine.h>
#include <utl/GeometryUtilities.h>
#include <utl/Plane.h>
#include <utl/Point.h>
#include <utl/Vector.h>
#include <utl/ShineUnits.h>
#include <utl/TrackParameters.h>

#include <modutils/StraightTrackingTypes.h>
#include <modutils/ParabolicTrackingTypes.h>

namespace evt{
  namespace rec{
    class Track;
    class Vertex;
    class VertexTrack;
  }
  class RecEvent;
}

namespace modutils {

  typedef evt::Index<evt::rec::Cluster> ClusterIndex;

  /// This class accompanies the KalmanFilterWB class for track parameter storage and interaction. Track parameters are stored
  /// in the {x, y, a, b, q/|p|} representation (a = x-z slope, b = y-z slope), to prevent unnecessary conversion to/from the 
  /// NA61 representation during SHINE reconstruction. The track parameters will be stored in RecEvent for access between local 
  /// tracking, global tracking and vertexing modules (TODO!).
  class KalmanTrackParameters {
    
  public:
    enum EKalmanTrackParameters {
				 eX = 0,
				 eY,
				 eA,
				 eB,
				 eQOverP,
				 eNTrackParameters,
				 eZ,
				 eNFreeParameters = eNTrackParameters
    };
    
    ///Default constructor.
    KalmanTrackParameters() : 
      fX(0),
      fY(0),
      fZ(0),
      fA(0),
      fB(0),
      fQOverP(0),
      fChi2(0),
      fNdf(0)
    { 
      fCovarianceMatrix.SetExtent(5);
      fCovarianceMatrixRepresentation = utl::TrackParameters::eKalmanTrackParameters;
    }
    
    ///Constructor using explicit parameters.
    KalmanTrackParameters(const double x,
			  const double y,
			  const double z,
			  const double a,
			  const double b,
			  const double qOverP,
			  const utl::CovarianceMatrix& covarianceMatrix,
			  const utl::TrackParameters::ECoordinateSystem coordinateSystem) : 
      fX(x),
      fY(y),
      fZ(z),
      fA(a),
      fB(b),
      fQOverP(qOverP),
      fChi2(0),
      fNdf(0),
      fCovarianceMatrix(covarianceMatrix),
      fCovarianceMatrixRepresentation(coordinateSystem)
    { }

    ///Constructor using SHINE track object.
    KalmanTrackParameters(const evt::rec::Track& track)
    { 
      //Set momentum at the start of the track.
      const utl::Point& position = track.GetMomentumPoint();
      const utl::Vector& momentum = track.GetMomentum();
      fX = position.GetX();
      fY = position.GetY();
      fZ = position.GetZ();
      fQOverP = track.GetCharge()/momentum.GetMag();
      fA = momentum.GetX()/momentum.GetZ();
      fB = momentum.GetY()/momentum.GetZ();
      fCovarianceMatrix = track.GetCovarianceMatrix();
      fChi2 = track.GetChi2();
      fNdf = track.GetNdf();
      fCovarianceMatrixRepresentation = utl::TrackParameters::eNA61;
    }

    ///Constructor using SHINE VertexTrack and Vertex object.
    KalmanTrackParameters(const evt::rec::VertexTrack& vertexTrack, const evt::rec::Vertex& vertex)
    { 
      //Set momentum at the start of the vertexTrack.
      const utl::Point& position = vertex.GetPosition();
      const utl::Vector& momentum = vertexTrack.GetMomentum();
      fX = position.GetX();
      fY = position.GetY();
      fZ = position.GetZ();
      fQOverP = vertexTrack.GetCharge()/momentum.GetMag();
      fA = momentum.GetX()/momentum.GetZ();
      fB = momentum.GetY()/momentum.GetZ();
      fCovarianceMatrix = vertexTrack.GetCovarianceMatrix();
      fChi2 = vertexTrack.GetChi2();
      fNdf = vertexTrack.GetNdf();
      fCovarianceMatrixRepresentation = utl::TrackParameters::eNA61;
    }

    ///Setters.
    void SetX(const double x) { fX = x; }
    void SetY(const double y) { fY = y; }
    void SetZ(const double z) { fZ = z; }
    void SetA(const double a) { fA = a; }
    void SetB(const double b) { fB = b; }
    void SetQOverP(const double qOverP) { fQOverP = qOverP; }
    void SetChi2(const double chi2) { fChi2 = chi2; }
    void SetNdf(const double ndf) { fNdf = ndf; }
    void SetCovarianceMatrix(const utl::CovarianceMatrix& cov) { fCovarianceMatrix = cov; }
    void SetTrackLength(const double length) { fTrackLength = length; }
    void SetAverageXResidual(const double average) { fAverageXResidual = average; }
    void SetAverageYResidual(const double average) { fAverageYResidual = average; }
    void SetTrimmedClusters(const std::vector<modutils::ClusterIndex>& clusters)
    { fTrimmedClusters = clusters; }
    
    //Translate straight track fit result to KalmanTrackParameters.
    void SetStraightFitParameters(const modutils::StraightTrackFitResult& fitResult,
				  const utl::Plane& referencePlane) {
      //Set position.
      fX = fitResult.fTrackParameters.fMX;
      fY = fitResult.fTrackParameters.fMY;
      fZ = 0;
      //Record slopes as calculated in straight track fit.
      fA = fitResult.fTrackParameters.fNX;
      fB = fitResult.fTrackParameters.fNY;
      //We don't know anything about charge or momentum.
      fQOverP = 0;

      //Copy covariance matrix elements.
      fCovarianceMatrix[eX][eX] = fabs(fitResult.fTrackParametersCovar.fMXMX);
      fCovarianceMatrix[eX][eA] = fabs(fitResult.fTrackParametersCovar.fMXNX);
      fCovarianceMatrix[eA][eA] = fabs(fitResult.fTrackParametersCovar.fNXNX);

      fCovarianceMatrix[eY][eY] = fabs(fitResult.fTrackParametersCovar.fMYMY);
      fCovarianceMatrix[eY][eB] = fabs(fitResult.fTrackParametersCovar.fMYNY);
      fCovarianceMatrix[eB][eB] = fabs(fitResult.fTrackParametersCovar.fNYNY);

      //Copy chi2 and NDF.
      fChi2 = fitResult.fChi2;
      fNdf = 2*fitResult.fNClusters - 4;

      //Copy average residuals.
      SetAverageXResidual(fitResult.fAverageXResidual);
      SetAverageYResidual(fitResult.fAverageYResidual);

      //Parameters exist on reference plane.
      ExtrapolateStraightFitToPlane(referencePlane);
    }
    
    //Translate parabolic track fit result to KalmanTrackParameters.
    void SetParabolicFitParameters(const modutils::ParabolicTrackFitResult& fitResult,
				   const utl::Plane& referencePlane) {
      //Set position.
      fX = fitResult.fTrackParameters.fCX;
      fY = fitResult.fTrackParameters.fBX;
      fZ = 0;
      //Record slopes as calculated in parabolic track fit.
      fA = 2*fitResult.fTrackParameters.fAX*fZ + fitResult.fTrackParameters.fBX;
      fB = fitResult.fTrackParameters.fAY;
      //We don't know anything about charge or momentum.
      fQOverP = 0;

      //Copy covariance matrix elements. FIXME: Re-calculate X-Z covariances!
      fCovarianceMatrix[eX][eX] = fabs(fitResult.fTrackParametersCovar.fAXAX);
      fCovarianceMatrix[eX][eA] = fabs(fitResult.fTrackParametersCovar.fAXBX);
      fCovarianceMatrix[eA][eA] = fabs(fitResult.fTrackParametersCovar.fBXBX);

      fCovarianceMatrix[eY][eY] = fabs(fitResult.fTrackParametersCovar.fAYAY);
      fCovarianceMatrix[eY][eB] = fabs(fitResult.fTrackParametersCovar.fAYBY);
      fCovarianceMatrix[eB][eB] = fabs(fitResult.fTrackParametersCovar.fBYBY);

      //Copy chi2 and NDF.
      fChi2 = fitResult.fChi2;
      fNdf = 2*fitResult.fNClusters - 4;

      //Store parabolic fit result.
      fParabolicFitResult = fitResult;
      
      //Copy average residuals.
      SetAverageXResidual(fitResult.fAverageXResidual);
      SetAverageYResidual(fitResult.fAverageYResidual);

      //Parameters exist on reference plane.
      ExtrapolateParabolicFitToPlane(referencePlane);
    }

    ///Getters.
    double GetX() const { return fX; }
    double GetY() const { return fY; }
    double GetZ() const { return fZ; }
    double GetA() const { return fA; }
    double GetB() const { return fB; }
    double GetQOverP() const { return fQOverP; }
    double GetTrackLength() const { return fTrackLength; }
    double GetXUncertainty() { return sqrt(fabs(fCovarianceMatrix[eX][eX])); }
    double GetYUncertainty() { return sqrt(fabs(fCovarianceMatrix[eY][eY])); }
    double GetAUncertainty() { return sqrt(fabs(fCovarianceMatrix[eA][eA])); }
    double GetBUncertainty() { return sqrt(fabs(fCovarianceMatrix[eB][eB])); }
    double GetQOverPUncertainty() { return sqrt(fabs(fCovarianceMatrix[eQOverP][eQOverP])); }
    double GetMomentumUncertainty() { return GetPMag()*GetPMag()*sqrt(fabs(fCovarianceMatrix[eQOverP][eQOverP])); }
    const utl::Point GetPosition() const { return utl::Point(fX,fY,fZ); }
    const utl::Vector GetMomentum() const { 
      const double largeMomentum = 100000*utl::GeV;
      const double pMag = (fQOverP == 0) ? largeMomentum : fabs(1/fQOverP);
      const double pZ = pMag*(1/sqrt(1 + fA*fA + fB*fB));
      const double pX = fA*pZ;
      const double pY = fB*pZ;
      return utl::Vector(pX,pY,pZ); 
    }
    int GetCharge() const {
      const double largeMomentum = 100000*utl::GeV;
      const double pMag = (fQOverP == 0) ? largeMomentum : fabs(1/fQOverP);
      return round(fQOverP*pMag);
    }
    double GetPMag() {
      const double largeMomentum = 100000*utl::GeV;
      const double pMag = (fQOverP == 0) ? largeMomentum : fabs(1/fQOverP);
      return pMag;
    }
    double GetChi2() const { return fChi2; }
    unsigned int GetNdf() const { return fNdf; }
    double GetAverageXResidual() const { return fAverageXResidual; }
    double GetAverageYResidual() const { return fAverageYResidual; }
    const std::vector<modutils::ClusterIndex>& GetTrimmedClusters() const { return fTrimmedClusters; }
    const utl::CovarianceMatrix
    GetCovarianceMatrix(const utl::TrackParameters::ECoordinateSystem coordinateSystem =
			utl::TrackParameters::eKalmanTrackParameters) const;
    
    //Get straight track fit result from KalmanTrackParameters. One for yz, one for xz.
    //XZ.
    const utl::Fitted2DLine GetStraightFitResultXZ() const {
      utl::Fitted2DLine xzLine;
      //Administrative things.
      xzLine.SetChi2(fChi2);
      xzLine.SetNdf(fNdf);
      //Covariance matrix translation.
      utl::CovarianceMatrix covariance(2);
      using namespace utl;
      covariance[EFitted2DLinePars::eIntercept][EFitted2DLinePars::eIntercept] = fCovarianceMatrix[eX][eX];
      covariance[EFitted2DLinePars::eSlope][EFitted2DLinePars::eIntercept]     = fCovarianceMatrix[eA][eX];
      covariance[EFitted2DLinePars::eSlope][EFitted2DLinePars::eSlope]         = fCovarianceMatrix[eA][eA];
      xzLine.SetCovarianceMatrix(covariance);
      //Set slope.
      xzLine.SetSlope(fA);
      //Calculate and set intercept.
      double intercept = fX - fA*fZ;
      xzLine.SetIntercept(intercept);
      return xzLine;
    }
    //YZ.
    const utl::Fitted2DLine GetStraightFitResultYZ() const {
      utl::Fitted2DLine yzLine;
      //Administrative things.
      yzLine.SetChi2(fChi2);
      yzLine.SetNdf(fNdf);
      //Covariance matrix translation.
      utl::CovarianceMatrix covariance(2);
      using namespace utl;
      covariance[EFitted2DLinePars::eIntercept][EFitted2DLinePars::eIntercept] = fCovarianceMatrix[eY][eY];
      covariance[EFitted2DLinePars::eSlope][EFitted2DLinePars::eIntercept]     = fCovarianceMatrix[eB][eY];
      covariance[EFitted2DLinePars::eSlope][EFitted2DLinePars::eSlope]         = fCovarianceMatrix[eB][eB];
      yzLine.SetCovarianceMatrix(covariance);
      //Set slope.
      yzLine.SetSlope(fB);
      //Calculate and set intercept.
      double intercept = fY - fB*fZ;
      yzLine.SetIntercept(intercept);
      return yzLine;
    }

    //Extrapolate using straight fit results.
    void ExtrapolateStraightFitToZ(const double z) {
      //Get the fitted lines.
      const utl::Fitted2DLine xzLine = GetStraightFitResultXZ();
      const utl::Fitted2DLine yzLine = GetStraightFitResultYZ();
      //Extrapolate to this z.
      const double extrapolatedX = xzLine.GetValueAtZ(z);
      const double extrapolatedY = yzLine.GetValueAtZ(z);
      //Extrapolate uncertainties.
      const double varianceX = xzLine.GetVarianceAtZ(z-fZ);
      const double varianceY = yzLine.GetVarianceAtZ(z-fZ);
      //Update track position.
      fX = extrapolatedX;
      fY = extrapolatedY;
      fZ = z;
      //Update covariance matrix.
      fCovarianceMatrix[eX][eX] = varianceX;
      fCovarianceMatrix[eY][eY] = varianceY;
    }

    //Extrapolate using straight fit results.
    void ExtrapolateStraightFitToPlane(const utl::Plane& plane) {
      //Find intersection point of line and plane.
      const utl::Line trajectory(GetPosition(),GetMomentum());
      const utl::Point& intersection = utl::GeometryUtilities::Intersection(plane, trajectory);
      const double z = intersection.GetZ();
      //Get the fitted lines.
      const utl::Fitted2DLine xzLine = GetStraightFitResultXZ();
      const utl::Fitted2DLine yzLine = GetStraightFitResultYZ();
      //Extrapolate to this z.
      const double extrapolatedX = xzLine.GetValueAtZ(z);
      const double extrapolatedY = yzLine.GetValueAtZ(z);
      //Extrapolate uncertainties.
      const double varianceX = xzLine.GetVarianceAtZ(z-fZ);
      const double varianceY = yzLine.GetVarianceAtZ(z-fZ);
      //Update track position.
      fX = extrapolatedX;
      fY = extrapolatedY;
      fZ = z;
      //Update covariance matrix.
      fCovarianceMatrix[eX][eX] = varianceX;
      fCovarianceMatrix[eY][eY] = varianceY;
    }

    //Extrapolate using parabolic fit results.
    void ExtrapolateParabolicFitToZ(const double z) {
      //Extrapolate uncertainties.
      //FIXME!! Do this right!
      //Update covariance matrix.
      fCovarianceMatrix[eX][eX] += fCovarianceMatrix[eX][eX]*(z-fZ)*(z-fZ);
      fCovarianceMatrix[eY][eY] += fCovarianceMatrix[eY][eY]*(z-fZ)*(z-fZ);
      //Extrapolate position using the ParabolicFitResult method.
      fX = fParabolicFitResult.ExtrapolateXToZ(z);
      fY = fParabolicFitResult.ExtrapolateYToZ(z);
      fZ = z;
    }
    
    //Extrapolate using parabolic fit result to plane.
    void ExtrapolateParabolicFitToPlane(const utl::Plane& plane,
					const double tolerance = 1e-3*utl::cm) {
      while (fabs(utl::GeometryUtilities::Distance(GetPosition(),plane)) > tolerance) {
	//Find intersection point of line and plane.
	const utl::Line trajectory(GetPosition(),GetMomentum());
	const utl::Point& intersection = utl::GeometryUtilities::Intersection(plane, trajectory);
	const double z = intersection.GetZ();
	//Extrapolate uncertainties.
	//FIXME!! Do this right!
	//Update covariance matrix.
	fCovarianceMatrix[eX][eX] += fCovarianceMatrix[eX][eX]*(z-fZ)*(z-fZ);
	fCovarianceMatrix[eY][eY] += fCovarianceMatrix[eY][eY]*(z-fZ)*(z-fZ);
	//Extrapolate position using the ParabolicFitResult method.
	fX = fParabolicFitResult.ExtrapolateXToZ(z);
	fY = fParabolicFitResult.ExtrapolateYToZ(z);
	fZ = z;
      }
    }

    //Printing function. 
    std::string const Print() const{
      std::stringstream msg;
      msg << std::setw(6);
      msg << "X = " << fX << " "
	  << "Y = " << fY << " "
	  << "Z = " << fZ << " "
	  << "A = " << fA << " "
	  << "B = " << fB << " "
	  << "QOverP = " << fQOverP << " ";
      msg << std::endl;
      for(unsigned int i=0; i!= fCovarianceMatrix.GetExtent(); i++)
	{
	  for(unsigned int j=0; j!= fCovarianceMatrix.GetExtent(); j++)
	    {
	      msg << fCovarianceMatrix[i][j] << " ";
	    }
	  msg << std::endl;
	}
      return msg.str();
    }
    
    
  private:
    //Covariance matrix inversion.
    void InvertCovarianceMatrix(const utl::CovarianceMatrix& covarianceMatrix);

    ///Momentum point global x-position.
    double fX;
    ///Momentum point global y-position.
    double fY;
    ///Momentum point global z-position.
    double fZ;
    ///Track x-z slope at momentum point.
    double fA;
    ///Track y-z slope at momentum point.
    double fB;
    ///Track charge over momentum magnitude.
    double fQOverP;
    ///Parabolic track fit storage.
    modutils::ParabolicTrackFitResult fParabolicFitResult;
    ///Fit chi2.
    double fChi2;
    ///Fit ndf.
    unsigned int fNdf;
    ///Track length.
    double fTrackLength;
    ///Average track residuals.
    double fAverageXResidual;
    double fAverageYResidual;
    ///Trimmed cluster indices.
    std::vector<ClusterIndex> fTrimmedClusters;
    ///Fit covariance matrix.
    mutable utl::CovarianceMatrix fCovarianceMatrix;
    ///Covariance matrix for inverting fit covariance matrix.
    mutable utl::CovarianceMatrix fInverseCovarianceMatrix;
    ///Current covariance matrix representation.
    mutable utl::TrackParameters::ECoordinateSystem fCovarianceMatrixRepresentation;
    ///Covariance matrix is singular flag.
    bool fCovarianceMatrixIsSingular;
    
  };
}

#endif
