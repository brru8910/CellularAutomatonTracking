#include "KalmanTrackParameters.h"

#include <TMatrixD.h>

using namespace std;

namespace modutils {
  
  void KalmanTrackParameters::InvertCovarianceMatrix(const utl::CovarianceMatrix& covarianceMatrix)
  {  
    double determinant = 0;
    TMatrixD cov(eNFreeParameters, eNFreeParameters);
    for (unsigned int i = 0; i < eNFreeParameters; ++i)
      for (unsigned int j = 0; j < eNFreeParameters; ++j)
        cov(i, j) = covarianceMatrix[i][j];
    
    cov.InvertFast(&determinant);
    if (fabs(determinant) < 1e-60) {
      fCovarianceMatrixIsSingular = true;
      cerr << " cov is singular!! " << determinant << endl;
      return;
    }
    
    fInverseCovarianceMatrix.SetExtent(eNFreeParameters);
    for (unsigned int i = 0; i < eNFreeParameters; ++i)
      for (unsigned int j = i; j < eNFreeParameters; ++j)
        fInverseCovarianceMatrix[i][j] = cov(i, j);
  }
  
  using namespace utl;
  const CovarianceMatrix
  KalmanTrackParameters::GetCovarianceMatrix(const TrackParameters::ECoordinateSystem newSys) const
  {
    
    if (newSys == fCovarianceMatrixRepresentation)
      return fCovarianceMatrix;
    //Convert from SHINE basis to KalmanTrackParamters basis.
    else if (newSys == TrackParameters::eKalmanTrackParameters &&
	     fCovarianceMatrixRepresentation == TrackParameters::eNA61) {

      //Helper quantities.
      const double c0 = 1 + fA*fA;

      //Convert to NA61 coordinates.
      const double q = fQOverP/fabs(fQOverP);
      const double qOverPXZ = q*fQOverP*sqrt((c0+fB*fB)/(c0));
      const double tanLambda = fB/sqrt(c0);
      const double psi = atan(1/fA);

      //Helper quantities.
      const double c1 = 1 + tanLambda*tanLambda;
      const double c2 = sin(psi);
      
      //Compute Jacobian. Full matrix:
      TMatrixD jacobian(eNFreeParameters, eNFreeParameters);
      
      //Diagonal terms.
      jacobian(TrackParameters::eX,KalmanTrackParameters::eX) = 1;
      jacobian(TrackParameters::eY,KalmanTrackParameters::eY) = 1;
      
      //Momentum terms: a.
      jacobian(TrackParameters::ePsi,KalmanTrackParameters::eA) =
	-1./(c2*c2);
      
      //Momentum terms: b.
      jacobian(TrackParameters::eTanLambda,KalmanTrackParameters::eB) =
	1./c2;
      jacobian(TrackParameters::ePsi,KalmanTrackParameters::eB) =
	-1.*tanLambda*cos(psi)/(c2*c2);
      
      //Momentum terms: q/|p|.      
      jacobian(TrackParameters::eQOverPXZ,KalmanTrackParameters::eQOverP) =
	q/sqrt(c1);
      jacobian(TrackParameters::eTanLambda,KalmanTrackParameters::eQOverP) =
        -1.*q*fabs(qOverPXZ)/sqrt(c1*c1*c1);
      

      //Rotate.
      TMatrixD covariance(eNFreeParameters, eNFreeParameters);
      for (unsigned int i = 0; i < eNFreeParameters; ++i)
	for (unsigned int j = 0; j < eNFreeParameters; ++j)
	  covariance(i, j) = fCovarianceMatrix[i][j];
      
      TMatrixD jacobianTransposed(TMatrixD::kTransposed,jacobian);
      TMatrixD rotatedCovariance = jacobianTransposed*covariance*jacobian;

      //Copy and return result.
      for (unsigned int i = 0; i < eNFreeParameters; ++i)
      	for (unsigned int j = i; j < eNFreeParameters; ++j)
      	  fCovarianceMatrix[i][j] = rotatedCovariance(i, j);

      //Set representation flag.
      fCovarianceMatrixRepresentation = TrackParameters::eKalmanTrackParameters;
      
      return fCovarianceMatrix;
    }
    //Convert from KalmanTrackParamters basis to SHINE basis.
    else if (newSys == TrackParameters::eNA61 && 
	     fCovarianceMatrixRepresentation == TrackParameters::eKalmanTrackParameters) {
      
      //Helper quantities.
      const double c0 = fA*fB;
      const double c1 = 1 + fA*fA;
      const double c2 = 1 + fA*fA + fB*fB;
      
      //Compute Jacobian. Full matrix:
      TMatrixD jacobian(eNFreeParameters, eNFreeParameters);
      
      //Diagonal terms.
      jacobian(KalmanTrackParameters::eX,TrackParameters::eX) = 1;
      jacobian(KalmanTrackParameters::eY,TrackParameters::eY) = 1;
      //Momentum terms: q/p_{xz}.
      jacobian(KalmanTrackParameters::eA,TrackParameters::eQOverPXZ) =
	(-1.*c0*fB*fQOverP/(c1*c1))*sqrt(c1/c2);
      jacobian(KalmanTrackParameters::eB,TrackParameters::eQOverPXZ) =
	(fB*fQOverP/c1)*sqrt(c1/c2);
      jacobian(KalmanTrackParameters::eQOverP,TrackParameters::eQOverPXZ) =
	sqrt(c2/c1);
      
      //Momentum terms: tan(lambda).
      jacobian(KalmanTrackParameters::eA,TrackParameters::eTanLambda) =
	-1.*c0/(c1*sqrt(c1));
      jacobian(KalmanTrackParameters::eB,TrackParameters::eTanLambda) =
	1./sqrt(c1);
      
      //Momentum terms: psi.
      jacobian(KalmanTrackParameters::eA,TrackParameters::ePsi) =
	-1./c1;

      //Rotate.
      TMatrixD covariance(eNFreeParameters, eNFreeParameters);
      for (unsigned int i = 0; i < eNFreeParameters; ++i)
	for (unsigned int j = 0; j < eNFreeParameters; ++j)
	  covariance(i, j) = fCovarianceMatrix[i][j];
      
      TMatrixD jacobianTransposed(TMatrixD::kTransposed,jacobian);
      TMatrixD rotatedCovariance = jacobianTransposed*covariance*jacobian;

      //Copy and return result.
      for (unsigned int i = 0; i < eNFreeParameters; ++i)
      	for (unsigned int j = i; j < eNFreeParameters; ++j)
      	  fCovarianceMatrix[i][j] = rotatedCovariance(i, j);

      //Set representation flag.
      fCovarianceMatrixRepresentation = TrackParameters::eKalmanTrackParameters;
      
      return fCovarianceMatrix;
    }
    else
      throw DoesNotComputeException("Improper track parameter representation provided!");
  }
  
};
