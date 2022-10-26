#ifndef _modutils_KalmanFilter_h_
#define _modutils_KalmanFilter_h_


/****************************************
 * Primary Author: Maksym Zyzak         *
 * Developed by: Wojciech Brylinski     *
 * Date: 9.03.2018                      *
 ****************************************/


#include "FitClasses.h"
#include "KalmanTrackParameters.h"
#include <evt/RecEvent.h>
#include <det/Detector.h>
#include <det/MagneticField.h>
#include <modutils/LocalTrackingTypes.h>
#include <utl/Plane.h>

#include <fwk/VModule.h>
#include <boost/utility.hpp>

#include <iomanip>
#include <math.h>

namespace modutils {
  
  enum EFitDirection {
    eForward = 0,
    eBackward = 1
  };
  
  //Calculate physical constants.
  //Get atmospheric pressure and calculate MIP dE/dx.
  const double mipAr = 1.521*utl::MeV*utl::cm2/utl::g;
  const double mipCO2 = 1.822*utl::MeV*utl::cm2/utl::g;
  //TPC gas density at NTP.
  const double density = 1.205E-3*(utl::g/utl::cm3);
  const double kMIPDEDXArCO2 = (0.9*mipAr + 0.1*mipCO2)*density;
  //estimated radiation length for Ar/CO2 95/5, assuming normal conditions (NTP)
  const double radLengthArCO2 = 20.49*utl::g/utl::cm2;
  const double kRadLengthArCO2 = radLengthArCO2/density;
      
  class KalmanFilter {
    
  public:
    
  KalmanFilter() : 
    fRecEvent(0) {}
		
    ~KalmanFilter() {}

    /// Set initial values
    void Init(const evt::RecEvent& recEvent) { 
      fRecEvent = &recEvent; 
      fMaximumStepSize = 10*utl::cm;
      fMinimumStepSize = 0.1*utl::cm;
      fExtrapolationTolerance = 0.1*utl::cm;
    }
    void SetTrackParameters(TrackFit& track,
			    modutils::KalmanTrackParameters& trackParameters);
    void ExtrapolateToGlobalZ(double T[],
			      Cov& C,
			      const double z_out);
    bool ExtrapolateToGlobalZ(const modutils::KalmanTrackParameters& trackParametersIn,
			      modutils::KalmanTrackParameters& trackParametersOut,
			      const double z);		
    bool ExtrapolateToPlane(const modutils::KalmanTrackParameters& trackParametersIn,
			    modutils::KalmanTrackParameters& trackParametersOut,
			    const utl::Plane& plane,
			    const double tolerance = 1e-3*utl::cm);
    bool KalmanFit(const ClusterVector& clusterIndices,
		   modutils::KalmanTrackParameters& trackParameters);
    bool KalmanFitWithVertex(const ClusterVector& clusterIndices,
    			     const evt::rec::Vertex& vertex,
    			     modutils::KalmanTrackParameters& seedParameters);
    ClusterVector SortClustersByZ(const ClusterVector& clusterIndices);
    
    //Fit function declaration.
    double rcp(double x);	
    void ExtrapolateALight(double T [],
			   Cov& C,
			   const double z_out,
			   double qp0,
			   FieldRegion& F);
    void Filter(TrackFit& track,
		const Measurement& measurement);
    void Smooth(TrackFit& t,
		std::vector<Measurement>& measurements,
		const EFitDirection fitDirection);
    void SmoothMeasurement(TrackFit& track,
			   const Measurement& measurement);
    void FilterFirst(TrackFit& track,
		     const Measurement& measurement,
		     const Measurement& measurementNext);
    void GuessVec(TrackFit& t,
		  std::vector<Measurement>& measurements,
		  const EFitDirection direction);
    void FitTrack(TrackFit& t,
		  std::vector<Measurement>& measurements,
		  const EFitDirection direction,
		  const bool seedFit = false);
    void Fit(TrackFit& t,
	     std::vector<Measurement>& measurements,
	     const EFitDirection direction,
	     const bool seedFit = false);
    void CalculateChi2(TrackFit& track,
		       const std::vector<Measurement>& measurements);
    void SetMaximumStepSize(const double& maxStep)
    { fMaximumStepSize = maxStep; }
    void SetMinimumStepSize(const double& minStep)
    { fMinimumStepSize = minStep; }
    void SetExtrapolationTolerance(const double& tolerance)
    { fExtrapolationTolerance = tolerance; }

  private:
    /// pointer to the reconstructed event
    const evt::RecEvent* fRecEvent;
    ///Maximum step size for Kalman Filter extrapolation.
    //Default value is set in Init() function.
    double fMaximumStepSize;
    ///Minimum step size for Kalman Filter extrapolation.
    //Default value is set in Init() function.
    double fMinimumStepSize;
    ///Tolerance for reaching target z during Kalman Filter extrapolation.
    //Default value is set in Init() function.
    double fExtrapolationTolerance;

  };

} // end of namespace modutils

#endif // _modutils_KalmanFitter_h_
