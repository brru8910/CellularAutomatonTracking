#ifndef _modutils_FitClasses_h_
#define _modutils_FitClasses_h_

#include <math.h>
#include <vector>
#include <evt/RecEvent.h>
#include <utl/ShineUnits.h>

/****************************************
 * Primary Author: Maksym Zyzak         *
 * Developed by: Wojciech Brylinski     *
 * Date: 9.03.2018                      *
 ****************************************/



namespace modutils {

	struct HitInfo{
	  float fCosPhi, fSinPhi, fSigma2, fSigma216, fSigma2XY;
	};


	struct FieldVector{
	  float X, Y, Z; //field at current point
	  void Combine( FieldVector &H, const float &w ){
		X+= w*(H.X-X);
		Y+= w*(H.Y-Y);
		Z+= w*(H.Z-Z);
	  }
	};


	struct FieldSlice{
	  float X[10], Y[10], Z[10]; 
	  //float z;
	  FieldSlice(){ for( int i=0; i<10; i++ ) X[i]=Y[i]=Z[i]=0; }

	  void GetField( const float &x, const float &y, float &Hx, float &Hy, float &Hz ){

		float x2 = x*x;
		float y2 = y*y;
		float xy = x*y;
		float x3 = x2*x;
		float y3 = y2*y;
		float xy2 = x*y2;
		float x2y = x2*y;
		//float x4 = x3*x;
		//float xy3 = x*y3;
		//float x2y2 = x*xy2;
		//float x3y = x*x2y;
		//float y4 = y*y3;

		Hx = X[0] +X[1]*x +X[2]*y +X[3]*x2 +X[4]*xy +X[5]*y2 +X[6]*x3 +X[7]*x2y +X[8]*xy2 +X[9]*y3;
		//+ X[10]*x4 + X[11]*x3y +X[12]*x2y2 +X[13]*xy3 + X[14]*y4;
		Hy = Y[0] +Y[1]*x +Y[2]*y +Y[3]*x2 +Y[4]*xy +Y[5]*y2 +Y[6]*x3 +Y[7]*x2y +Y[8]*xy2 +Y[9]*y3;
		//+ Y[10]*x4 + Y[11]*x3y +Y[12]*x2y2 +Y[13]*xy3 + Y[14]*y4;
		Hz = Z[0] +Z[1]*x +Z[2]*y +Z[3]*x2 +Z[4]*xy +Z[5]*y2 +Z[6]*x3 +Z[7]*x2y +Z[8]*xy2 +Z[9]*y3;
		//+ Z[10]*x4 + Z[11]*x3y +Z[12]*x2y2 +Z[13]*xy3 + Z[14]*y4;
	  }

	  void GetField( const float &x, const float &y, FieldVector &H ){
		GetField( x, y, H.X, H.Y, H.Z );
	  }
	};


	struct FieldRegion{
	  float x0, x1, x2 ; // Hx(Z) = x0 + x1*(Z-z) + x2*(Z-z)^2
	  float y0, y1, y2 ; // Hy(Z) = y0 + y1*(Z-z) + y2*(Z-z)^2
	  float z0, z1, z2 ; // Hz(Z) = z0 + z1*(Z-z) + z2*(Z-z)^2
	  float z;

	  FieldRegion(){
		x0 = x1 = x2 = y0 = y1 = y2 = z0 = z1 = z2 = z = 0.;
	  }

	  void Set( const FieldVector &H0, const float &H0z, 
			const FieldVector &H1, const float &H1z, 
			const FieldVector &H2, const float &H2z){
		z = H0z;
		float dz1 = H1z-H0z, dz2 = H2z-H0z;
		float det = (float)1./(dz1*dz2*(dz2-dz1));
		float w21 = -dz2*det;
		float w22 = dz1*det;
		float w11 = -dz2*w21;
		float w12 = -dz1*w22;
		
		float dH1 = H1.X - H0.X;
		float dH2 = H2.X - H0.X;
		x0 = H0.X;
		x1 = dH1*w11 + dH2*w12 ;
		x2 = dH1*w21 + dH2*w22 ;
		  
		dH1 = H1.Y - H0.Y;
		dH2 = H2.Y - H0.Y;
		y0 = H0.Y;
		y1 = dH1*w11 + dH2*w12 ;
		y2 = dH1*w21 + dH2*w22  ;

		dH1 = H1.Z - H0.Z;
		dH2 = H2.Z - H0.Z;
		z0 = H0.Z;
		z1 = dH1*w11 + dH2*w12 ;
		z2 = dH1*w21 + dH2*w22 ;         
	  }
	  
	  void Shift( const float &Z0){
		float dz = Z0-z;
		float x2dz = x2*dz;
		float y2dz = y2*dz;
		float z2dz = z2*dz;
		z = Z0;
		x0+= (x1 + x2dz)*dz;
		x1+= x2dz+x2dz; 
		y0+= (y1 + y2dz)*dz;
		y1+= y2dz+y2dz; 
		z0+= (z1 + z2dz)*dz;
		z1+= z2dz+z2dz; 
	  }

	};

	struct Cov{
	  double C00, 
	    C10, C11, 
	    C20, C21, C22, 
	    C30, C31, C32, C33, 
	    C40, C41, C42, C43, C44;
	};

        struct TrackFit {
	  double T[7]; // x, y, tx, ty, qp, z, length
	  Cov   C;    // cov matr.
	  double Chi2;
	  double NDF;
	  double fXResidualAverage;
	  double fXResidualSigma;
	  double fYResidualAverage;
	  double fYResidualSigma;
	  std::vector<evt::Index<evt::rec::Cluster> > fTrimmedClusters;
	  std::vector<int> fTrimmedClustersInt;
	};

	struct Measurement{
	  double fX, fY, fZ;
	  double fWeight;
	  double fSigmaX, fSigma2X;
	  double fSigmaY, fSigma2Y;
	  double fSigma2XY;
	  double fSy;
	  TrackFit fTrackEstimate;
	  unsigned int fClusterIndex;
	  
	  FieldSlice fMap;
	  FieldVector fField;
	  Measurement() { 
	    fSigmaX = 200*utl::um; fSigma2X = fSigmaX*fSigmaX;
	    fSigmaY = 200*utl::um; fSigma2Y = fSigmaY*fSigmaY;
	  }
	};
	
}

#endif
