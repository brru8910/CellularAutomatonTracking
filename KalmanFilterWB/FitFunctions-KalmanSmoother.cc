#include <det/Detector.h>
#include <utl/ErrorLogger.h>
#include <utl/ShineUnits.h>
#include <det/MagneticField.h>
#include <sstream>
#include "FitClasses.h"
#include "KalmanFilterWB.h"

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TString.h"
#include "TMultiGraph.h"
#include "TMarker.h"
#include "TLine.h"
#include "TGraphErrors.h"

#include <Math/SMatrix.h>

/****************************************
 * Primary Author: Maksym Zyzak         *
 * Developed by: Wojciech Brylinski     *
 * Date: 9.03.2018                      *
 ****************************************/

namespace modutils {

  using namespace std;

  double KalmanFilter::rcp( double x ){ return 1./x; }	

  void KalmanFilter::ExtrapolateALight
  ( 
   double T [], // input track parameters (x,y,tx,ty,Q/p)
   Cov &C,     // input covariance matrix
   const double z_out  , // extrapolate to this z position
   double       qp0    , // use Q/p linearisation at this value
   FieldRegion &F
    )
  {
    //
    //  Part of the analytic extrapolation formula with error (0.000299792458*B*dz)^4/4!
    //  
    double  
      c1 = 1., c2 = 2., c3 = 3., c4 = 4., c6 = 6., c9 = 9., c15 = 15., c18 = 18., c45 = 45.,
      c2i = 1./2., c3i = 1./3., c6i = 1./6., c12i = 1./12.;

    double qp = T[4];

    // if (qp != qp) ERROR("Initial q/p is nan!!");
    // if (qp0 != qp0) ERROR("Initial q/p is nan!!");
	  
    double dz = (z_out - T[5]);
    double dz2 = dz*dz;
    double dz3 = dz2*dz;
	  
    double T0 = T[0];
    double T1 = T[1];
    double T2 = T[2];
    double T3 = T[3];
    // double T4 = T[4];
    double T5 = T[5];

    // construct coefficients 

    double x   = T[2]; // tx !!
    double y   = T[3]; // ty !!
    double xx  = x*x;
    double xy = x*y;
    double yy = y*y;
    double y2 = y*c2;
    double x2 = x*c2;
    double x4 = x*c4;
    double xx31 = xx*c3+c1;
    double xx159 = xx*c15+c9;

    double Ay = -xx-c1;
    double Ayy = x*(xx*c3+c3);
    double Ayz = -c2*xy; 
    double Ayyy = -(c15*xx*xx+c18*xx+c3);

    double Ayy_x = c3*xx31;
    double Ayyy_x = -x4*xx159;

    double Bx = yy+c1; 
    double Byy = y*xx31; 
    double Byz = c2*xx+c1;
    double Byyy = -xy*xx159;

    double Byy_x = c6*xy;
    double Byyy_x = -y*(c45*xx+c9);
    double Byyy_y = -x*xx159;

    // end of coefficients calculation

    double t2   = c1 + xx + yy;
    double t    = sqrt( t2 );
    double h    = qp0*0.000299792458;
    double ht   = h*t;

    // get field integrals
    double ddz = T5-F.z;
    double Fx0 = F.x0 + F.x1*ddz + F.x2*ddz*ddz;
    double Fx1 = (F.x1 + c2*F.x2*ddz)*dz;
    double Fx2 = F.x2*dz2;
    double Fy0 = F.y0 + F.y1*ddz + F.y2*ddz*ddz;
    double Fy1 = (F.y1 + c2*F.y2*ddz)*dz;
    double Fy2 = F.y2*dz2;
    double Fz0 = F.z0 + F.z1*ddz + F.z2*ddz*ddz;
    double Fz1 = (F.z1 + c2*F.z2*ddz)*dz;
    double Fz2 = F.z2*dz2;

    double sx = ( Fx0 + Fx1*c2i + Fx2*c3i );
    double sy = ( Fy0 + Fy1*c2i + Fy2*c3i );
    double sz = ( Fz0 + Fz1*c2i + Fz2*c3i );

    double Sx = ( Fx0*c2i + Fx1*c6i + Fx2*c12i );
    double Sy = ( Fy0*c2i + Fy1*c6i + Fy2*c12i );
    double Sz = ( Fz0*c2i + Fz1*c6i + Fz2*c12i );

    double syz;
    { 
      double 
	d = 1./360., 
	c00 = 30.*6.*d, c01 = 30.*2.*d,   c02 = 30.*d,
	c10 = 3.*40.*d, c11 = 3.*15.*d,   c12 = 3.*8.*d, 
	c20 = 2.*45.*d, c21 = 2.*2.*9.*d, c22 = 2.*2.*5.*d;
      syz = Fy0*( c00*Fz0 + c01*Fz1 + c02*Fz2) 
	+   Fy1*( c10*Fz0 + c11*Fz1 + c12*Fz2) 
	+   Fy2*( c20*Fz0 + c21*Fz1 + c22*Fz2) ;
    }

    double Syz;
    {
      double 
	d = 1./2520., 
	c00 = 21.*20.*d, c01 = 21.*5.*d, c02 = 21.*2.*d,
	c10 =  7.*30.*d, c11 =  7.*9.*d, c12 =  7.*4.*d, 
	c20 =  2.*63.*d, c21 = 2.*21.*d, c22 = 2.*10.*d;
      Syz = Fy0*( c00*Fz0 + c01*Fz1 + c02*Fz2 ) 
	+   Fy1*( c10*Fz0 + c11*Fz1 + c12*Fz2 ) 
	+   Fy2*( c20*Fz0 + c21*Fz1 + c22*Fz2 ) ;
    }

    double syy  = sy*sy*c2i;
    double syyy = syy*sy*c3i;

    double Syy ;   
    {
      double  
	d= 1./2520., c00= 420.*d, c01= 21.*15.*d, c02= 21.*8.*d,
	c03= 63.*d, c04= 70.*d, c05= 20.*d;
      Syy =  Fy0*(c00*Fy0+c01*Fy1+c02*Fy2) + Fy1*(c03*Fy1+c04*Fy2) + c05*Fy2*Fy2 ;
    }
	  
    double Syyy;
    {
      double 
	d = 1./181440., 
	c000 =   7560*d, c001 = 9*1008*d, c002 = 5*1008*d, 
	c011 = 21*180*d, c012 = 24*180*d, c022 =  7*180*d, 
	c111 =    540*d, c112 =    945*d, c122 =    560*d, c222 = 112*d;
      double Fy22 = Fy2*Fy2;
      Syyy = Fy0*( Fy0*(c000*Fy0+c001*Fy1+c002*Fy2)+ Fy1*(c011*Fy1+c012*Fy2)+c022*Fy22 )
	+    Fy1*( Fy1*(c111*Fy1+c112*Fy2)+c122*Fy22) + c222*Fy22*Fy2                  ;
    }
	  
	  
    double sA1   = sx*xy   + sy*Ay   + sz*y ;
    double sA1_x = sx*y - sy*x2 ;
    double sA1_y = sx*x + sz ;

    double sB1   = sx*Bx   - sy*xy   - sz*x ;
    double sB1_x = -sy*y - sz ;
    double sB1_y = sx*y2 - sy*x ;

    double SA1   = Sx*xy   + Sy*Ay   + Sz*y ;
    double SA1_x = Sx*y - Sy*x2 ;
    double SA1_y = Sx*x + Sz;

    double SB1   = Sx*Bx   - Sy*xy   - Sz*x ;
    double SB1_x = -Sy*y - Sz;
    double SB1_y = Sx*y2 - Sy*x;


    double sA2   = syy*Ayy   + syz*Ayz ;
    double sA2_x = syy*Ayy_x - syz*y2 ;
    double sA2_y = -syz*x2 ;
    double sB2   = syy*Byy   + syz*Byz  ;
    double sB2_x = syy*Byy_x + syz*x4 ;
    double sB2_y = syy*xx31 ;
	  
    double SA2   = Syy*Ayy   + Syz*Ayz ;
    double SA2_x = Syy*Ayy_x - Syz*y2 ;
    double SA2_y = -Syz*x2 ;
    double SB2   = Syy*Byy   + Syz*Byz ;
    double SB2_x = Syy*Byy_x + Syz*x4 ;
    double SB2_y = Syy*xx31 ;

    double sA3   = syyy*Ayyy  ;
    double sA3_x = syyy*Ayyy_x;
    double sB3   = syyy*Byyy  ;
    double sB3_x = syyy*Byyy_x;
    double sB3_y = syyy*Byyy_y;

	 
    double SA3   = Syyy*Ayyy  ;
    double SA3_x = Syyy*Ayyy_x;
    double SB3   = Syyy*Byyy  ;
    double SB3_x = Syyy*Byyy_x;
    double SB3_y = Syyy*Byyy_y;

    double ht1 = ht*dz;
    double ht2 = ht*ht*dz2;
    double ht3 = ht*ht*ht*dz3;
    double ht1sA1 = ht1*sA1;
    double ht1sB1 = ht1*sB1;
    double ht1SA1 = ht1*SA1;
    double ht1SB1 = ht1*SB1;
    double ht2sA2 = ht2*sA2;
    double ht2SA2 = ht2*SA2;
    double ht2sB2 = ht2*sB2;
    double ht2SB2 = ht2*SB2;
    double ht3sA3 = ht3*sA3;
    double ht3sB3 = ht3*sB3; 
    double ht3SA3 = ht3*SA3;
    double ht3SB3 = ht3*SB3;

    T[0] = T0 + (x + ht1SA1 + ht2SA2 + ht3SA3)*dz ;
    T[1] = T1 + (y + ht1SB1 + ht2SB2 + ht3SB3)*dz ;
    T[2] = T2 + ht1sA1 + ht2sA2 + ht3sA3;
    T[3] = T3 + ht1sB1 + ht2sB2 + ht3sB3;
    T[5]+= dz;
	  
    double ctdz  = 0.000299792458*t*dz;
    double ctdz2 = 0.000299792458*t*dz2;

    double dqp = qp - qp0;
    double t2i = c1*rcp(t2);// /t2;
    double xt2i = x*t2i;
    double yt2i = y*t2i;
    double tmp0 = ht1SA1 + c2*ht2SA2 + c3*ht3SA3;
    double tmp1 = ht1SB1 + c2*ht2SB2 + c3*ht3SB3;
    double tmp2 = ht1sA1 + c2*ht2sA2 + c3*ht3sA3;
    double tmp3 = ht1sB1 + c2*ht2sB2 + c3*ht3sB3;

    double j02 = dz*(c1 + xt2i*tmp0 + ht1*SA1_x + ht2*SA2_x + ht3*SA3_x);
    double j12 = dz*(     xt2i*tmp1 + ht1*SB1_x + ht2*SB2_x + ht3*SB3_x);
    double j22 =     c1 + xt2i*tmp2 + ht1*sA1_x + ht2*sA2_x + ht3*sA3_x ;
    double j32 =          xt2i*tmp3 + ht1*sB1_x + ht2*sB2_x + ht3*sB3_x ;
		
    double j03 = dz*(     yt2i*tmp0 + ht1*SA1_y + ht2*SA2_y );
    double j13 = dz*(c1 + yt2i*tmp1 + ht1*SB1_y + ht2*SB2_y + ht3*SB3_y );
    double j23 =          yt2i*tmp2 + ht1*sA1_y + ht2*sA2_y  ;
    double j33 =     c1 + yt2i*tmp3 + ht1*sB1_y + ht2*sB2_y + ht3*sB3_y ;
		
    double j04 = ctdz2*( SA1 + c2*ht1*SA2 + c3*ht2*SA3 );
    double j14 = ctdz2*( SB1 + c2*ht1*SB2 + c3*ht2*SB3 );
    double j24 = ctdz *( sA1 + c2*ht1*sA2 + c3*ht2*sA3 );
    double j34 = ctdz *( sB1 + c2*ht1*sB2 + c3*ht2*sB3 );
	  
    // extrapolate inverse momentum
	  
    T[0]+=j04*dqp;
    T[1]+=j14*dqp;
    T[2]+=j24*dqp;
    T[3]+=j34*dqp;

    //          covariance matrix transport 
    // Commenting out unused matrix entries. Fixes warnings and may speed up algorithm...sightly.
	 
    double c42 = C.C42, c43 = C.C43;

    double cj00 = C.C00 + C.C20*j02 + C.C30*j03 + C.C40*j04;
    // double cj10 = C.C10 + C.C21*j02 + C.C31*j03 + C.C41*j04;
    double cj20 = C.C20 + C.C22*j02 + C.C32*j03 + c42*j04;
    double cj30 = C.C30 + C.C32*j02 + C.C33*j03 + c43*j04;
	 
    double cj01 = C.C10 + C.C20*j12 + C.C30*j13 + C.C40*j14;
    double cj11 = C.C11 + C.C21*j12 + C.C31*j13 + C.C41*j14;
    double cj21 = C.C21 + C.C22*j12 + C.C32*j13 + c42*j14;
    double cj31 = C.C31 + C.C32*j12 + C.C33*j13 + c43*j14;
	 
    // double cj02 = C.C20*j22 + C.C30*j23 + C.C40*j24;
    // double cj12 = C.C21*j22 + C.C31*j23 + C.C41*j24;
    double cj22 = C.C22*j22 + C.C32*j23 + c42*j24;
    double cj32 = C.C32*j22 + C.C33*j23 + c43*j24;
	 
    // double cj03 = C.C20*j32 + C.C30*j33 + C.C40*j34;
    // double cj13 = C.C21*j32 + C.C31*j33 + C.C41*j34;
    double cj23 = C.C22*j32 + C.C32*j33 + c42*j34;
    double cj33 = C.C32*j32 + C.C33*j33 + c43*j34;

    C.C40+= c42*j02 + c43*j03 + C.C44*j04; // cj40
    C.C41+= c42*j12 + c43*j13 + C.C44*j14; // cj41
    C.C42 = c42*j22 + c43*j23 + C.C44*j24; // cj42
    C.C43 = c42*j32 + c43*j33 + C.C44*j34; // cj43

    C.C00 = cj00 + j02*cj20 + j03*cj30 + j04*C.C40;
    C.C10 = cj01 + j02*cj21 + j03*cj31 + j04*C.C41;
    C.C11 = cj11 + j12*cj21 + j13*cj31 + j14*C.C41;

    C.C20 = j22*cj20 + j23*cj30 + j24*C.C40 ;
    C.C30 = j32*cj20 + j33*cj30 + j34*C.C40 ;
    C.C21 = j22*cj21 + j23*cj31 + j24*C.C41 ;
    C.C31 = j32*cj21 + j33*cj31 + j34*C.C41 ;
    C.C22 = j22*cj22 + j23*cj32 + j24*C.C42 ;
    C.C32 = j32*cj22 + j33*cj32 + j34*C.C42 ;
    C.C33 = j32*cj23 + j33*cj33 + j34*C.C43 ;

  }




  void KalmanFilter::Filter(TrackFit &track,
			    HitInfo &info,
			    double u,
			    double w ,
			    const bool calculateChi2)
  {

    // convert input
    double *T = track.T;
    Cov &C = track.C;

    register double wi, zeta = 0, zetawi = 0, HCH;
	  
    register double F0, F1, F2, F3, F4;
    register double  K1, K2, K3, K4;

    //   double wi, zeta, zetawi, HCH;
    //   
    //   double F0, F1, F2, F3, F4;
    //   double  K1, K2, K3, K4;
	  
    zeta = info.fCosPhi*T[0] + info.fSinPhi*T[1] - u;
    // F = CH'

    F0 = info.fCosPhi*C.C00 + info.fSinPhi*C.C10;
    F1 = info.fCosPhi*C.C10 + info.fSinPhi*C.C11;
    
    HCH = ( F0*info.fCosPhi + F1*info.fSinPhi );

    F2 = info.fCosPhi*C.C20 + info.fSinPhi*C.C21;
    F3 = info.fCosPhi*C.C30 + info.fSinPhi*C.C31;
    F4 = info.fCosPhi*C.C40 + info.fSinPhi*C.C41;
    
    //This factor favors the measurement error over the covariance matrix error.
    //It's needed to keep the parameter error from getting too small...
    double errorFactor = sqrt(2*log(2));
    wi = w * rcp(errorFactor*info.fSigma2 + HCH);
    
    zetawi = zeta * wi;

    if (calculateChi2) {
      track.Chi2 += zeta * zetawi;
      track.NDF  += w;
    }

    K1 = F1*wi;
    K2 = F2*wi;
    K3 = F3*wi;
    K4 = F4*wi;

    T[0]-= F0*zetawi;
    T[1]-= F1*zetawi;
    T[2]-= F2*zetawi;
    T[3]-= F3*zetawi;
    T[4]-= F4*zetawi;

    C.C00-= F0*F0*wi;
    C.C10-= K1*F0;
    C.C11-= K1*F1;
    C.C20-= K2*F0;
    C.C21-= K2*F1;
    C.C22-= K2*F2;
    C.C30-= K3*F0;
    C.C31-= K3*F1;
    C.C32-= K3*F2;
    C.C33-= K3*F3;
    C.C40-= K4*F0;
    C.C41-= K4*F1;
    C.C42-= K4*F2;
    C.C43-= K4*F3;
    C.C44-= K4*F4;
  }

  void KalmanFilter::FilterFirst(TrackFit& track,
				 const Measurement& measurement)
  {

    Cov& C = track.C;
    const double w1 = 1.-measurement.fWeight;
    const double sigma2X = measurement.fWeight*measurement.fSigma2X;
    const double sigma2Y = measurement.fWeight*measurement.fSigma2Y;
    const double sigma2XY = measurement.fWeight*measurement.fSigma2XY;

    // initialize covariance matrix. Use reasonably large guesses.
    C.C00= sigma2X*100; 
    C.C10= sigma2XY*100;   C.C11= sigma2Y*100; 
    C.C20= 0.0;      C.C21= 0.0;      C.C22= 0.01;
    C.C30= 0.0;      C.C31= 0.0;      C.C32= 0.0; C.C33= 0.01;
    C.C40= 0.0;      C.C41= 0.0;      C.C42= 0.0; C.C43= 0.0; C.C44= 0.01;

    track.T[0] = measurement.fWeight*measurement.fX + w1*track.T[0];
    track.T[1] = measurement.fWeight*measurement.fY + w1*track.T[1];
    track.NDF = -3.0;
    track.Chi2 = 0.0;
  }

  void KalmanFilter::AddMaterial(TrackFit& track,
				 const Measurement& measurement,
				 const EFitDirection fitDirection,
				 double& qp0)
  {
    //Track parameters. Multiple scattering adds noise on xz & yz angles only.
    double tx = track.T[2];
    double ty = track.T[3];
    
    //Calculate expected variance of scattering angle using Highland formula (see PDG).
    //Assume track to be a pion, since we can't do PID before tracking.
    const double pionMass = 139.6*utl::MeV;
    const double pionMass2 = pionMass*pionMass;
    const double momentum = fabs(1/qp0);
    const double energy = sqrt(momentum*momentum + pionMass2);
    const double gamma = energy/pionMass;
    const double beta = momentum/(gamma*pionMass);
    
    const double theta =
      (13.6*utl::MeV/(beta*momentum))*
      sqrt(measurement.fRadThickness)*(1 + 0.038*measurement.fLogRadThickness);

    //Calculate covariance increase for angle terms.
    const double pOverPZ2 = 1 + tx*tx + ty*ty;
    const double varC22 = theta*theta*(1 + tx*tx)*pOverPZ2;
    const double varC32 = theta*theta*tx*ty*pOverPZ2;
    const double varC33 = theta*theta*(1 + ty*ty)*pOverPZ2;

    //Store.
    Cov &C = track.C;
    C.C22 += varC22;
    C.C32 += varC32; C.C33 += varC33; 
    
    //Modify track energy estimate based on dE/dx in this material. Assume
    //MIP, since this needs to be done very quickly.
    //If forward fitting is re-implemented, this must change to subtracting
    //off the energy.
    if (fitDirection == eBackward) {
      track.T[4] = (qp0 > 0) ?
	1./(momentum + measurement.fMIPDEDX*measurement.fThickness) :
	-1./(momentum + measurement.fMIPDEDX*measurement.fThickness);
    }
    else if (fitDirection == eForward) {
      track.T[4] = (qp0 > 0) ?
	1./(momentum - measurement.fMIPDEDX*measurement.fThickness) :
	-1./(momentum - measurement.fMIPDEDX*measurement.fThickness);
    }
    
  }

  void KalmanFilter::GuessVec(TrackFit& track,
			      const vector<Measurement>& measurements)
  {
    double *T = track.T;

    double A0 = 0, A1=0.0, A2=0.0, A3=0.0, A4=0.0, A5=0.0, a0=0, a1=0.0, a2=0.0,
      b0 =0, b1=0.0, b2=0.0;
    double z0=0, x, y, z, S, w, wz, wS;

    const Measurement& lastMeasurement = measurements.back(); 
    z0 = lastMeasurement.fZ;
    w = lastMeasurement.fWeight;
    A0 = w;
    a0 = w*lastMeasurement.fX;
    b0 = w*lastMeasurement.fY;
      
    ostringstream s1;
    for(unsigned int i = 0; i < measurements.size(); ++i){
      Measurement measurement = measurements[i];
      //Only include measurements within magnetic field region.
      //Otherwise, ROC estimate is underestimated (momentum is
      //overestimated) due to MTPC cluster positions.
      //FIXME: Improve. MTPC clusters still influence guess.
      if (fabs(measurement.fSy) == 0)
      	continue;
      x = measurement.fX;
      y = measurement.fY;
      w = measurement.fWeight;
      z = measurement.fZ - z0;
      S = measurement.fSy;
      wz = w*z;
      wS = w*S;
      A0+=w; 
      A1+=wz;  A2+=wz*z;
      A3+=wS;  A4+=wS*z; A5+=wS*S;
      a0+=w*x; a1+=wz*x; a2+=wS*x;
      b0+=w*y; b1+=wz*y; b2+=wS*y;
      s1 << "Measurement Z: " << measurement.fZ << ". Sy: " << S << endl;
    }
    INFO(s1);
    double A3A3 = A3*A3;
    double A3A4 = A3*A4;
    double A1A5 = A1*A5;
    double A2A5 = A2*A5;
    double A4A4 = A4*A4;
	  
    double det;
    double detInverse = -A2*A3A3 + A1*( A3A4+A3A4 - A1A5) + A0*(A2A5-A4A4);
    double largeValue = 100000;
    if (detInverse != 0) det = rcp(detInverse);
    else { det = largeValue; }

    double Ai0 = ( -A4A4 + A2A5 );
    double Ai1 = (  A3A4 - A1A5 );
    double Ai2 = ( -A3A3 + A0*A5 );
    double Ai3 = ( -A2*A3 + A1*A4 );
    double Ai4 = (  A1*A3 - A0*A4 );
    double Ai5 = ( -A1*A1 + A0*A2 );

    double L, L1;
    T[0] = (Ai0*a0 + Ai1*a1 + Ai3*a2)*det;
    T[2] = (Ai1*a0 + Ai2*a1 + Ai4*a2)*det;
    double txtx1 = 1.f+T[2]*T[2];
    L    = (Ai3*a0 + Ai4*a1 + Ai5*a2)*det *rcp(txtx1);

    L1 = L*T[2];
    A1 = A1 + A3*L1;
    A2 = A2 + ( A4 + A4 + A5*L1 )*L1;
    b1+= b2 * L1;
    det = rcp(-A1*A1+A0*A2);


    T[1] = (  A2*b0 - A1*b1 )*det;
    T[3] = ( -A1*b0 + A0*b1 )*det;
    T[4] = -L*(1./0.000299792458)*sqrt(txtx1 +T[3]*T[3]);
    T[5] = z0;

    ostringstream s;
    s << "GuessVec results: p = " << 1/fabs(T[4]) << ". [" << T[0]
      << ", " << T[1]
      << ", " << T[2]
      << ", " << T[3]
      << ", " << T[4] << "]";
    			 INFO(s);
  }   

  void
  KalmanFilter::Fit(TFile* file,
		    TrackFit& t,
		    vector<Measurement>& measurements,
		    const bool seedFit,
		    int eventNumber,
		    int trackIndex)
  {
    vector<double> xs;
    vector<double> ys;
    vector<double> zs;
    vector<double> xUncertaintiesLow;
    vector<double> xUncertaintiesHigh;
    vector<double> yUncertaintiesLow;
    vector<double> yUncertaintiesHigh;

    vector<double> xsFiltered;
    vector<double> ysFiltered;
    vector<double> zsFiltered;
    vector<double> xUncertaintiesFilteredLow;
    vector<double> xUncertaintiesFilteredHigh;
    vector<double> yUncertaintiesFilteredLow;
    vector<double> yUncertaintiesFilteredHigh;
    
    vector<double> clusterXs;
    vector<double> clusterYs;
    vector<double> clusterXUncertainties;
    vector<double> clusterYUncertainties;

    //Track length.
    double length = 0;

    double c16 = 16.;
    // upstream

    //Criteria for using previous fit as a seed:
    //   --User selected seedFit = true
    //   --Track q/|p| was not 0 or NaN
    //   --Track momentum was not 1 TeV (case for MTPC tracks)
    if (!seedFit              ||
	t.T[4] == 0           ||
	t.T[4] != t.T[4]      ||
	t.T[4] == 1/1*utl::TeV ) {
      GuessVec(t, measurements);
    }
    
    if (t.T[4] != t.T[4]) ERROR("Initial q/p is nan!!");
	  
    //Begin backward fitting.
    for (unsigned int i = measurements.size() - 1; i > 0; --i ){
      HitInfo xInfo, yInfo;
      xInfo.fCosPhi = 1.0;
      xInfo.fSinPhi = 0.0;
      xInfo.fSigma2  = measurements[i-1].fSigma2X;
      xInfo.fSigma216 = xInfo.fSigma2*c16;
      yInfo.fCosPhi = 0.0;
      yInfo.fSinPhi = 1.0;
      yInfo.fSigma2  = measurements[i-1].fSigma2Y;
      yInfo.fSigma216 = yInfo.fSigma2*c16;

      //Don't use hits with 0 uncertainty.
      if (yInfo.fSigma2 == 0 || xInfo.fSigma2 == 0) continue;
	      
      Measurement& measurement = measurements[i];
      Measurement& measurementNext = measurements[i-1];
	
      if( i == measurements.size() ) {
    	FilterFirst( t, measurement );

    	measurement.fTrackEstimate = t;
	  
    	xs.push_back(t.T[0]);
    	ys.push_back(t.T[1]);
    	xUncertaintiesLow.push_back(t.T[0]-sqrt(t.C.C00));
    	xUncertaintiesHigh.push_back(t.T[0]+sqrt(t.C.C00));
    	yUncertaintiesLow.push_back(t.T[1]-sqrt(t.C.C11));
    	yUncertaintiesHigh.push_back(t.T[1]+sqrt(t.C.C11));
	  
    	clusterXs.push_back(measurement.fX);
    	clusterYs.push_back(measurement.fY);
    	zs.push_back(measurement.fZ);
	  
    	clusterXUncertainties.push_back(measurement.fSigmaX);
    	clusterYUncertainties.push_back(measurement.fSigmaY);
      }
	
      double x = t.T[0];
      double y = t.T[1];
      double z = t.T[5];
	  
      ExtrapolateToGlobalZ(t.T, t.C, measurementNext.fZ);
      AddMaterial( t, measurement, eBackward, t.T[4] );

      double nextX = t.T[0];
      double nextY = t.T[1];
      double nextZ = t.T[5];
	
      xs.push_back(nextX);
      ys.push_back(nextY);
      zs.push_back(nextZ);

      xUncertaintiesLow.push_back(nextX-sqrt(t.C.C00));
      xUncertaintiesHigh.push_back(nextX+sqrt(t.C.C00));
      yUncertaintiesLow.push_back(nextY-sqrt(t.C.C11));
      yUncertaintiesHigh.push_back(nextY+sqrt(t.C.C11));

      clusterXs.push_back(measurementNext.fX);
      clusterYs.push_back(measurementNext.fY);
	
      clusterXUncertainties.push_back(measurement.fSigmaX);
      clusterYUncertainties.push_back(measurement.fSigmaY);

    	length += sqrt((x - nextX)*(x - nextX) + 
    		       (y - nextY)*(y - nextY) + 
    		       (z - nextZ)*(z - nextZ));

      bool calculateChi2 = true;
      Filter( t, xInfo, measurementNext.fX, measurementNext.fWeight, calculateChi2);
      Filter( t, yInfo, measurementNext.fY, measurementNext.fWeight, calculateChi2);
	
      measurementNext.fTrackEstimate = t;

      xsFiltered.push_back(t.T[0]);
      ysFiltered.push_back(t.T[1]);
      zsFiltered.push_back(t.T[5]);
      xUncertaintiesFilteredLow.push_back(t.T[0]-sqrt(t.C.C00));
      xUncertaintiesFilteredHigh.push_back(t.T[0]+sqrt(t.C.C00));
      yUncertaintiesFilteredLow.push_back(t.T[1]-sqrt(t.C.C11));
      yUncertaintiesFilteredHigh.push_back(t.T[1]+sqrt(t.C.C11));
	
    } //End backward fitting.
    
    ostringstream s;
    for(unsigned int i = 0; i < measurements.size(); ++i){
      Measurement measurement = measurements[i];
      Cov& C = measurement.fTrackEstimate.C;
      s << endl << "Measurement " << i << endl
    	<<  "Track parameters: (q*p = " << 1./measurement.fTrackEstimate.T[4]
    	<< ") [" << measurement.fTrackEstimate.T[0]
    	<< " " << measurement.fTrackEstimate.T[1]
    	<< " " << measurement.fTrackEstimate.T[2]
    	<< " " << measurement.fTrackEstimate.T[3]
    	<< " " << measurement.fTrackEstimate.T[4]
    	<< "] " << endl;
      s << "Covariance:"
    	<< endl
    	<< C.C00 << "\n"
    	<< C.C10 << "   " << C.C11 << "\n"
    	<< C.C20 << "   " << C.C21 << "   " << C.C22 << "\n"
    	<< C.C30 << "   " << C.C31 << "   " << C.C32 << "   " << C.C33 << "\n"
    	<< C.C40 << "   " << C.C41 << "   " << C.C42 << "   " << C.C43 << "   " << C.C44 << "\n";
      
    }
    INFO(s);


    // //Begin forward smoothing.
    // for (unsigned int i = 0; i < measurements.size() - 1; ++i ){
    //   HitInfo xInfo, yInfo;
    //   xInfo.fCosPhi = 1.0;
    //   xInfo.fSinPhi = 0.0;
    //   xInfo.fSigma2  = measurements[i+1].fSigma2X;
    //   xInfo.fSigma216 = xInfo.fSigma2*c16;
    //   yInfo.fCosPhi = 0.0;
    //   yInfo.fSinPhi = 1.0;
    //   yInfo.fSigma2  = measurements[i+1].fSigma2Y;
    //   yInfo.fSigma216 = yInfo.fSigma2*c16;

    //   //Don't use hits with 0 uncertainty.
    //   if (yInfo.fSigma2 == 0 || xInfo.fSigma2 == 0) continue;
	      
    //   Measurement& measurement = measurements[i];
    //   Measurement& measurementNext = measurements[i+1];
	
    //   if(i == measurements.size() - 1) {
    // 	measurement.fTrackEstimate = t;
    //   }
	  
    //   ExtrapolateToGlobalZ(t.T, t.C, measurementNext.fZ);
    //   AddMaterial( t, measurementNext, eForward, t.T[4] );

    //   ostringstream s1;
    //   s1 << endl << "Extrapolated forward to measurment " << i << "." << endl
    // 	 <<  "Track parameters: (q*p = " << 1./measurement.fTrackEstimate.T[4]
    // 	 << ") [" << measurement.fTrackEstimate.T[0]
    // 	 << " " << measurement.fTrackEstimate.T[1]
    // 	 << " " << measurement.fTrackEstimate.T[2]
    // 	 << " " << measurement.fTrackEstimate.T[3]
    // 	 << " " << measurement.fTrackEstimate.T[4]
    // 	 << "] " << endl;
    //   s1 << "Covariance:"
    // 	 << endl
    // 	 << t.C.C00 << "\n"
    // 	 << t.C.C10 << "   " << t.C.C11 << "\n"
    // 	 << t.C.C20 << "   " << t.C.C21 << "   " << t.C.C22 << "\n"
    // 	 << t.C.C30 << "   " << t.C.C31 << "   " << t.C.C32 << "   " << t.C.C33 << "\n"
    // 	 << t.C.C40 << "   " << t.C.C41 << "   " << t.C.C42 << "   " << t.C.C43 << "   " << t.C.C44 << "\n";
    //   INFO(s1);

      
    //   Smooth( t, measurementNext);

    //   measurementNext.fTrackEstimate = t;
      
    // } //End forward smoothing.
    



    // // Begin KF debug fit drawing.
    
    // vector<double> trackXs;
    // vector<double> trackYs;
    // vector<double> trackZs;
    // const double startZ = zs.front();
    // const double stopZ = zs.back();
    // const double nTrackSteps = 5000;
    // const double stepSize = (stopZ - startZ)/nTrackSteps;
    // for (unsigned int step = 0; step < nTrackSteps; ++step) {
    //   const double currentZ = startZ + step*stepSize;
    //   ExtrapolateToGlobalZ(t.T, t.C, currentZ);
    //   trackXs.push_back(t.T[0]);
    //   trackYs.push_back(t.T[1]);
    //   trackZs.push_back(currentZ);
    // }
    
    // TLegend* legend = new TLegend(0.67, 0.8, 0.89, 0.89);
    // TMarker* clusterMarker = new TMarker(0,0,7);
    // clusterMarker->SetMarkerColor(kGreen);
    // TMarker* fitMarker = new TMarker(0,0,7);
    // fitMarker->SetMarkerColor(kBlue);
    // TMarker* filteredMarker = new TMarker(0,0,7);
    // filteredMarker->SetMarkerColor(kMagenta);
    // TLine* trackLine = new TLine;
    // trackLine->SetLineStyle(3);
    // trackLine->SetLineColor(kRed);
    // TLine* extrapolationLine = new TLine;
    // extrapolationLine->SetLineStyle(3);
    // extrapolationLine->SetLineColor(kBlue);
    // TLine* filteredLine = new TLine;
    // filteredLine->SetLineStyle(3);
    // filteredLine->SetLineColor(kMagenta);
    // legend->AddEntry(clusterMarker, "Clusters", "p");
    // legend->AddEntry(fitMarker, "KF Extrapolation", "p");
    // legend->AddEntry(filteredMarker, "KF Filtered Estimate", "p");
    // legend->AddEntry(extrapolationLine, "Extrapolation Uncertainty", "l");
    // legend->AddEntry(filteredLine, "Filtered Uncertainty", "l");
    // legend->AddEntry(trackLine, "Final Fit", "l");
    
    // TGraph xzFit(zs.size(),
    // 		 zs.data(),xs.data());
    // xzFit.SetMarkerColor(kBlue);
    // xzFit.SetMarkerStyle(8);
    // xzFit.Draw("ALP");
    // xzFit.SetMinimum(-400);
    // xzFit.SetMaximum(400);
    // xzFit.GetXaxis()->SetLimits(-600.,1000.);
    // TGraph xzFiltered(zsFiltered.size(),
    // 		      zsFiltered.data(),xsFiltered.data());
    // xzFiltered.SetMarkerColor(kMagenta);
    // xzFiltered.SetMarkerStyle(8);
    // xzFiltered.Draw("ALP");
    // xzFiltered.SetMinimum(-400);
    // xzFiltered.SetMaximum(400);
    // xzFiltered.GetXaxis()->SetLimits(-600.,1000.);
    // TGraph xzErrorHigh(zs.size(),
    // 		       zs.data(),xUncertaintiesHigh.data());
    // xzErrorHigh.SetLineColor(kBlue);
    // xzErrorHigh.SetLineStyle(2);
    // xzErrorHigh.SetLineWidth(2);
    // xzErrorHigh.Draw("L");
    // xzErrorHigh.SetMinimum(-400);
    // xzErrorHigh.SetMaximum(400);
    // xzErrorHigh.GetXaxis()->SetLimits(-600.,1000.);
    // TGraph xzErrorLow(zs.size(),
    // 		      zs.data(),xUncertaintiesLow.data());
    // xzErrorLow.SetLineColor(kBlue);
    // xzErrorLow.SetLineStyle(2);
    // xzErrorLow.SetLineWidth(2);
    // xzErrorLow.Draw("L");
    // xzErrorLow.SetMinimum(-400);
    // xzErrorLow.SetMaximum(400);
    // xzErrorLow.GetXaxis()->SetLimits(-600.,1000.);
    
    // TGraph xzErrorFilteredHigh(zsFiltered.size(),
    // 			       zsFiltered.data(),xUncertaintiesFilteredHigh.data());
    // xzErrorFilteredHigh.SetLineColor(kMagenta);
    // xzErrorFilteredHigh.SetLineStyle(2);
    // xzErrorFilteredHigh.SetLineWidth(2);
    // xzErrorFilteredHigh.Draw("L");
    // xzErrorFilteredHigh.SetMinimum(-400);
    // xzErrorFilteredHigh.SetMaximum(400);
    // xzErrorFilteredHigh.GetXaxis()->SetLimits(-600.,1000.);
    // TGraph xzErrorFilteredLow(zsFiltered.size(),
    // 			      zsFiltered.data(),xUncertaintiesFilteredLow.data());
    // xzErrorFilteredLow.SetLineColor(kMagenta);
    // xzErrorFilteredLow.SetLineStyle(2);
    // xzErrorFilteredLow.SetLineWidth(2);
    // xzErrorFilteredLow.Draw("L");
    // xzErrorFilteredLow.SetMinimum(-400);
    // xzErrorFilteredLow.SetMaximum(400);
    // xzErrorLow.GetXaxis()->SetLimits(-600.,1000.);
    
    // TGraphErrors xzTrack(trackZs.size(),
    // 			 trackZs.data(),trackXs.data());
    // xzTrack.SetMarkerColor(kRed);
    // xzTrack.SetMarkerStyle(1);
    // xzTrack.Draw("P");
    // xzTrack.SetMinimum(-400);
    // xzTrack.SetMaximum(400);
    // xzTrack.GetXaxis()->SetLimits(-600.,1000.);
    // TGraphErrors xzClusters(zs.size(),
    // 			    zs.data(),clusterXs.data(),
    // 			    vector<double>(zs.size(),0).data(),clusterXUncertainties.data());
    // xzClusters.SetMarkerColor(kGreen);
    // xzClusters.SetMarkerStyle(8);
    // xzClusters.Draw("P");
    // xzClusters.SetMinimum(-400);
    // xzClusters.SetMaximum(400);
    // xzClusters.GetXaxis()->SetLimits(-600.,1000.);
    
    // TMultiGraph xGraph;
    // TString nameStringXZ = Form("event%iTrack%iBackwardFitXZ",eventNumber,trackIndex);
    // TString titleStringXZ("XZ Track Fit;Z (cm);X (cm)");
    // xGraph.SetName(nameStringXZ);
    // xGraph.SetTitle(titleStringXZ);
    // xGraph.Add(&xzFit,"P");
    // xGraph.Add(&xzFiltered,"P");
    // xGraph.Add(&xzClusters,"P");
    // xGraph.Add(&xzTrack,"P");
    // xGraph.Add(&xzErrorHigh,"l");
    // xGraph.Add(&xzErrorLow,"l");
    // xGraph.Add(&xzErrorFilteredHigh,"l");
    // xGraph.Add(&xzErrorFilteredLow,"l");
    
    // TGraph yzFit(zs.size(),
    // 		 zs.data(),ys.data());
    // yzFit.SetMarkerColor(kBlue);
    // yzFit.SetMarkerStyle(8);
    // yzFit.Draw("ALP");
    // yzFit.SetMinimum(-400);
    // yzFit.SetMaximum(400);
    // yzFit.GetXaxis()->SetLimits(-600.,1000.);
    // TGraph yzFiltered(zsFiltered.size(),
    // 		      zsFiltered.data(),ysFiltered.data());
    // yzFiltered.SetMarkerColor(kMagenta);
    // yzFiltered.SetMarkerStyle(8);
    // yzFiltered.Draw("ALP");
    // yzFiltered.SetMinimum(-400);
    // yzFiltered.SetMaximum(400);
    // yzFiltered.GetXaxis()->SetLimits(-600.,1000.);
    // TGraph yzErrorHigh(zs.size(),
    // 		       zs.data(),yUncertaintiesHigh.data());
    // yzErrorHigh.SetLineColor(kBlue);
    // yzErrorHigh.SetLineStyle(2);
    // yzErrorHigh.SetLineWidth(2);
    // yzErrorHigh.Draw("L");
    // yzErrorHigh.SetMinimum(-400);
    // yzErrorHigh.SetMaximum(400);
    // yzErrorHigh.GetXaxis()->SetLimits(-600.,1000.);
    // TGraph yzErrorLow(zs.size(),
    // 		      zs.data(),yUncertaintiesLow.data());
    // yzErrorLow.SetLineColor(kBlue);
    // yzErrorLow.SetLineStyle(2);
    // yzErrorLow.SetLineWidth(2);
    // yzErrorLow.Draw("L");
    // yzErrorLow.SetMinimum(-400);
    // yzErrorLow.SetMaximum(400);
    // yzErrorLow.GetXaxis()->SetLimits(-600.,1000.);
    
    // TGraph yzErrorFilteredHigh(zsFiltered.size(),
    // 			       zsFiltered.data(),yUncertaintiesFilteredHigh.data());
    // yzErrorFilteredHigh.SetLineColor(kMagenta);
    // yzErrorFilteredHigh.SetLineStyle(2);
    // yzErrorFilteredHigh.SetLineWidth(2);
    // yzErrorFilteredHigh.Draw("L");
    // yzErrorFilteredHigh.SetMinimum(-400);
    // yzErrorFilteredHigh.SetMaximum(400);
    // yzErrorFilteredHigh.GetXaxis()->SetLimits(-600.,1000.);
    // TGraph yzErrorFilteredLow(zsFiltered.size(),
    // 			      zsFiltered.data(),yUncertaintiesFilteredLow.data());
    // yzErrorFilteredLow.SetLineColor(kMagenta);
    // yzErrorFilteredLow.SetLineStyle(2);
    // yzErrorFilteredLow.SetLineWidth(2);
    // yzErrorFilteredLow.Draw("L");
    // yzErrorFilteredLow.SetMinimum(-400);
    // yzErrorFilteredLow.SetMaximum(400);
    // yzErrorLow.GetXaxis()->SetLimits(-600.,1000.);
    
    // TGraphErrors yzTrack(trackZs.size(),
    // 			 trackZs.data(),trackYs.data());
    // yzTrack.SetMarkerColor(kRed);
    // yzTrack.SetMarkerStyle(1);
    // yzTrack.Draw("P");
    // yzTrack.SetMinimum(-400);
    // yzTrack.SetMaximum(400);
    // yzTrack.GetYaxis()->SetLimits(-600.,1000.);
    // TGraphErrors yzClusters(zs.size(),
    // 			    zs.data(),clusterYs.data(),
    // 			    vector<double>(zs.size(),0).data(),clusterYUncertainties.data());
    // yzClusters.SetMarkerColor(kGreen);
    // yzClusters.SetMarkerStyle(8);
    // yzClusters.Draw("P");
    // yzClusters.SetMinimum(-400);
    // yzClusters.SetMaximum(400);
    // yzClusters.GetYaxis()->SetLimits(-600.,1000.);
    // TMultiGraph yGraph;
    // TString nameStringYZ = Form("event%iTrack%iBackwardFitYZ",eventNumber,trackIndex);
    // TString titleStringYZ("YZ Track Fit;Z (cm);Y (cm)");
    // yGraph.SetName(nameStringYZ);
    // yGraph.SetTitle(titleStringYZ);
    // yGraph.Add(&yzFit,"P");
    // yGraph.Add(&yzFiltered,"P");
    // yGraph.Add(&yzClusters,"P");
    // yGraph.Add(&yzTrack,"P");
    // yGraph.Add(&yzErrorHigh,"l");
    // yGraph.Add(&yzErrorLow,"l");
    // yGraph.Add(&yzErrorFilteredHigh,"l");
    // yGraph.Add(&yzErrorFilteredLow,"l");
    
    // // if (debug) {
    // file->cd();
    // TCanvas xz;
    // xz.SetName("XZFit");
    // xz.SetTitle("XZFit");
    // xGraph.Draw("ALP");
    // legend->Draw("SAME");
    // TCanvas yz;
    // yz.SetName("YZFit");
    // yz.SetTitle("YZFit");
    // yGraph.Draw("ALP");
    // legend->Draw("SAME");

    // xz.Write();
    // yz.Write();
    // xGraph.Write();
    // yGraph.Write();
      
    // // }

    //   //End KF debug fit drawing.


    // //Begin forward fitting.
    // for (unsigned int i = 0; i < measurements.size() - 1; ++i ){
    //   HitInfo xInfo, yInfo;
    //   xInfo.fCosPhi = 1.0;
    //   xInfo.fSinPhi = 0.0;
    //   xInfo.fSigma2  = measurements[i+1].fSigma2X;
    //   xInfo.fSigma216 = xInfo.fSigma2*c16;
    //   yInfo.fCosPhi = 0.0;
    //   yInfo.fSinPhi = 1.0;
    //   yInfo.fSigma2  = measurements[i+1].fSigma2Y;
    //   yInfo.fSigma216 = yInfo.fSigma2*c16;

    //   //Don't use hits with 0 uncertainty.
    //   if (yInfo.fSigma2 == 0 || xInfo.fSigma2 == 0) continue;
	      
    //   Measurement& measurement = measurements[i];
    //   Measurement& measurementNext = measurements[i+1];
	
    //   if(i == 0) {
    // 	FilterFirst( t, measurement );
    //   }
	  
    //   ExtrapolateToGlobalZ(t.T, t.C, measurementNext.fZ);
    //   AddMaterial( t, measurementNext, eForward, t.T[4] );

    //   ostringstream s1;
    //   s1 << endl << "Extrapolated forward to measurment " << i << "." << endl
    // 	 <<  "Track parameters: (q*p = " << 1./measurement.fTrackEstimate.T[4]
    // 	 << ") [" << measurement.fTrackEstimate.T[0]
    // 	 << " " << measurement.fTrackEstimate.T[1]
    // 	 << " " << measurement.fTrackEstimate.T[2]
    // 	 << " " << measurement.fTrackEstimate.T[3]
    // 	 << " " << measurement.fTrackEstimate.T[4]
    // 	 << "] " << endl;
    //   s1 << "Covariance:"
    // 	 << endl
    // 	 << t.C.C00 << "\n"
    // 	 << t.C.C10 << "   " << t.C.C11 << "\n"
    // 	 << t.C.C20 << "   " << t.C.C21 << "   " << t.C.C22 << "\n"
    // 	 << t.C.C30 << "   " << t.C.C31 << "   " << t.C.C32 << "   " << t.C.C33 << "\n"
    // 	 << t.C.C40 << "   " << t.C.C41 << "   " << t.C.C42 << "   " << t.C.C43 << "   " << t.C.C44 << "\n";
    //   INFO(s1);


    //   bool calculateChi2 = true;
    //   Filter( t, xInfo, measurementNext.fX, measurementNext.fWeight, calculateChi2);
    //   Filter( t, yInfo, measurementNext.fY, measurementNext.fWeight, calculateChi2);
	
    //   // Smooth( t, measurementNext);

    //   measurementNext.fTrackEstimate = t;
      
    // } //End forward fitting.
    

    // //Begin backward smoothing.
    // for (unsigned int i = measurements.size() - 1; i > 0; --i ){
    //   HitInfo xInfo, yInfo;
    //   xInfo.fCosPhi = 1.0;
    //   xInfo.fSinPhi = 0.0;
    //   xInfo.fSigma2  = measurements[i+1].fSigma2X;
    //   xInfo.fSigma216 = xInfo.fSigma2*c16;
    //   yInfo.fCosPhi = 0.0;
    //   yInfo.fSinPhi = 1.0;
    //   yInfo.fSigma2  = measurements[i+1].fSigma2Y;
    //   yInfo.fSigma216 = yInfo.fSigma2*c16;

    //   //Don't use hits with 0 uncertainty.
    //   if (yInfo.fSigma2 == 0 || xInfo.fSigma2 == 0) continue;
	      
    //   Measurement& measurement = measurements[i];
    //   Measurement& measurementNext = measurements[i+1];
	
    //   if(i == measurements.size()) {
    // 	measurement.fTrackEstimate = t;
    //   }
	  
    //   ExtrapolateToGlobalZ(t.T, t.C, measurementNext.fZ);
    //   AddMaterial( t, measurementNext, eForward, t.T[4] );

    //   ostringstream s1;
    //   s1 << endl << "Extrapolated forward to measurment " << i << "." << endl
    // 	 <<  "Track parameters: (q*p = " << 1./measurement.fTrackEstimate.T[4]
    // 	 << ") [" << measurement.fTrackEstimate.T[0]
    // 	 << " " << measurement.fTrackEstimate.T[1]
    // 	 << " " << measurement.fTrackEstimate.T[2]
    // 	 << " " << measurement.fTrackEstimate.T[3]
    // 	 << " " << measurement.fTrackEstimate.T[4]
    // 	 << "] " << endl;
    //   s1 << "Covariance:"
    // 	 << endl
    // 	 << t.C.C00 << "\n"
    // 	 << t.C.C10 << "   " << t.C.C11 << "\n"
    // 	 << t.C.C20 << "   " << t.C.C21 << "   " << t.C.C22 << "\n"
    // 	 << t.C.C30 << "   " << t.C.C31 << "   " << t.C.C32 << "   " << t.C.C33 << "\n"
    // 	 << t.C.C40 << "   " << t.C.C41 << "   " << t.C.C42 << "   " << t.C.C43 << "   " << t.C.C44 << "\n";
    //   INFO(s1);
	
    //   Smooth( t, measurementNext);

    //   measurementNext.fTrackEstimate = t;
      
    // } //End backward smoothing.



    
    //Record length.
    t.T[6] = length;

  }

  void KalmanFilter::Smooth(TrackFit& track,
			    const Measurement& measurement) {
    
    const Cov& backwardC = measurement.fTrackEstimate.C;

    ROOT::Math::SVector<double,5> forwardState;
    ROOT::Math::SVector<double,5> backwardState;
    ROOT::Math::SMatrix<double,5,5,ROOT::Math::MatRepSym<double,5> > forwardCov;
    ROOT::Math::SMatrix<double,5,5,ROOT::Math::MatRepSym<double,5> > backwardCov;
    
    forwardState(0) = track.T[0];
    forwardState(1) = track.T[1];
    forwardState(2) = track.T[2];
    forwardState(3) = track.T[3];
    forwardState(4) = track.T[4];
    
    backwardState(0) = measurement.fTrackEstimate.T[0];
    backwardState(1) = measurement.fTrackEstimate.T[1];
    backwardState(2) = measurement.fTrackEstimate.T[2];
    backwardState(3) = measurement.fTrackEstimate.T[3];
    backwardState(4) = measurement.fTrackEstimate.T[4]; 

    forwardCov(0,0) = track.C.C00;
    forwardCov(1,0) = track.C.C10;
    forwardCov(2,0) = track.C.C20;
    forwardCov(3,0) = track.C.C30;
    forwardCov(4,0) = track.C.C40;
    forwardCov(1,1) = track.C.C11;
    forwardCov(2,1) = track.C.C21;
    forwardCov(3,1) = track.C.C31;
    forwardCov(4,1) = track.C.C41;
    forwardCov(2,2) = track.C.C22;
    forwardCov(3,2) = track.C.C32;
    forwardCov(4,2) = track.C.C42;
    forwardCov(3,3) = track.C.C33;
    forwardCov(4,3) = track.C.C43;
    forwardCov(4,4) = track.C.C44;
    
    backwardCov(0,0) = backwardC.C00;
    backwardCov(1,0) = backwardC.C10;
    backwardCov(2,0) = backwardC.C20;
    backwardCov(3,0) = backwardC.C30;
    backwardCov(4,0) = backwardC.C40;
    backwardCov(1,1) = backwardC.C11;
    backwardCov(2,1) = backwardC.C21;
    backwardCov(3,1) = backwardC.C31;
    backwardCov(4,1) = backwardC.C41;
    backwardCov(2,2) = backwardC.C22;
    backwardCov(3,2) = backwardC.C32;
    backwardCov(4,2) = backwardC.C42;
    backwardCov(3,3) = backwardC.C33;
    backwardCov(4,3) = backwardC.C43;
    backwardCov(4,4) = backwardC.C44;

    ROOT::Math::SMatrix<double,5,5,ROOT::Math::MatRepSym<double,5> > covSum = forwardCov + backwardCov;    
    
    int inversionSuccess = 1;
    ROOT::Math::SMatrix<double,5,5,ROOT::Math::MatRepSym<double,5> > covSumInverse =
			covSum.Inverse(inversionSuccess);
    
    ROOT::Math::SMatrix<double,5,5> kalmanGain = forwardCov * covSumInverse;

    ROOT::Math::SVector<double,5> smoothedState =
      forwardState + kalmanGain * (backwardState - forwardState);    

    ROOT::Math::SMatrix<double,5,5> smoothedCov =
      forwardCov - kalmanGain * forwardCov;
    
    track.T[0] = smoothedState(0);
    track.T[1] = smoothedState(1);
    track.T[2] = smoothedState(2);
    track.T[3] = smoothedState(3);
    track.T[4] = smoothedState(4);   
    
    track.C.C00 = smoothedCov(0,0);
    track.C.C10 = smoothedCov(1,0);
    track.C.C20 = smoothedCov(2,0);
    track.C.C30 = smoothedCov(3,0);
    track.C.C40 = smoothedCov(4,0);
    track.C.C11 = smoothedCov(1,1);
    track.C.C21 = smoothedCov(2,1);
    track.C.C31 = smoothedCov(3,1);
    track.C.C41 = smoothedCov(4,1);
    track.C.C22 = smoothedCov(2,2);
    track.C.C32 = smoothedCov(3,2);
    track.C.C42 = smoothedCov(4,2);
    track.C.C33 = smoothedCov(3,3);
    track.C.C43 = smoothedCov(4,3);
    track.C.C44 = smoothedCov(4,4);
    


    // ostringstream s;
    // s << endl << "Backward state: "
    //   << endl << backwardState
    //   << endl << ". Forward state: "
    //   << endl << forwardState
    //   << endl << ". Smoothed state: "
    //   << endl << smoothedState
    //   << endl << ". Backward cov: "
    //   << endl << backwardCov
    //   << endl << ". Forward cov: "
    //   << endl << forwardCov
    //   << endl << ". Smoothed cov: "
    //   << endl << smoothedCov
    //   << endl << ". Gain matrix: "
    //   << endl << kalmanGain;
    // INFO(s);
    
  }

  
}


