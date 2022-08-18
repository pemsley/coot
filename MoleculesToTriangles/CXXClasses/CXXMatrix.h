/*
 *  CXXMatrix.h
 *  iPhoneRibbons
 *
 *  Created by Martin Noble on 01/09/2008.
 *  Copyright 2008 LMB, Oxford University. All rights reserved.
 *
 */

#ifndef CXXMatrix_H
#define CXXMatrix_H

#include "MoleculesToTriangles/CXXSurface/CXXCoord.h"
//#include "Quaternion.h"

class CXXMatrix  {
private:
	FCXXCoord rows[4];
public:
	CXXMatrix(){
		for (int i=0; i<4; i++) rows[i][i] = 1.;
	};
    /*
	CXXMatrix(double d00, double d10, double d20, double d30, 
	          double d01, double d11, double d21, double d31,
	          double d02, double d12, double d22, double d32,
	          double d03, double d13, double d23, double d33){
		rows[0][0] = d00;
		rows[0][1] = d01;
		rows[0][2] = d02;
		rows[0][3] = d03;
		rows[1][0] = d10;
		rows[1][1] = d11;
		rows[1][2] = d12;
		rows[1][3] = d13;
		rows[2][0] = d20;
		rows[2][1] = d21;
		rows[2][2] = d22;
		rows[2][3] = d23;
		rows[3][0] = d30;
		rows[3][1] = d31;
		rows[3][2] = d32;
		rows[3][3] = d33;
	};
	*/
	CXXMatrix(double mat44[4][4]){
		for (int i=0; i<4; i++){
			for (int j=0; j<4; j++){
				rows[i][j] = mat44[i][j];
			}
		}
	};
    
	CXXMatrix(float *values){
		int i, k;
		for (i=0, k=0; i<4; i++){
			for (int j=0; j<4; j++){
				rows[i][j] = values[k++];
			}
		}
	};
    
	static CXXMatrix matrixWithQuaternion(const FCXXCoord &quaternion){
		double x=quaternion[0];
		double y=quaternion[1];
		double z=quaternion[2];
		double w=quaternion[3];
		return matrixWithQuaternion(x, y, z, w);
	};
    
	static CXXMatrix matrixWithQuaternion(double *quaternion){
		double &x=quaternion[0];
		double &y=quaternion[1];
		double &z=quaternion[2];
		double &w=quaternion[3];
		return matrixWithQuaternion(x, y, z, w);
	};
    
	static CXXMatrix matrixWithQuaternion(float *quaternion){
		double x=quaternion[0];
		double y=quaternion[1];
		double z=quaternion[2];
		double w=quaternion[3];
		return matrixWithQuaternion(x, y, z, w);
	};
    
	static CXXMatrix matrixWithQuaternion(double &x, double &y, double&z, double&w){
		
		CXXMatrix result;
		
		result[0][0] = w*w + x*x - y*y - z*z;
		result[0][1] = 2*x*y + 2*w*z;
		result[0][2] = 2*x*z - 2*w*y;
		
		result[1][0] = 2*x*y - 2*w*z;
		result[1][1] = w*w - x*x + y*y - z*z;
		result[1][2] = 2*y*z + 2*w*x;
		
		result[2][0] = 2*x*z + 2*w*y;
		result[2][1] = 2*y*z - 2*w*x;
		result[2][2] = w*w - x*x - y*y + z*z;
		return result;
	};	

	static CXXMatrix matrixWithGLRotation(double *gLRotation ){
		double &w=gLRotation[0];
		double &x=gLRotation[1];
		double &y=gLRotation[2];
		double &z=gLRotation[3];
		return matrixWithGLRotation(w, x, y, z);
	};
    
	static CXXMatrix matrixWithGLRotation(float *gLRotation){
		double w=gLRotation[0];
		double x=gLRotation[1];
		double y=gLRotation[2];
		double z=gLRotation[3];
		return matrixWithGLRotation(w, x, y, z);
	};
    
	static CXXMatrix matrixWithGLRotation(const FCXXCoord &glRotation){
        return matrixWithGLRotation(glRotation[0], glRotation[1], glRotation[2], glRotation[3]);
    };
    
	static CXXMatrix matrixWithGLRotation(const double &w, const double &x, const double &y, const double&z){
		CXXMatrix result;
		int i, j, k;
		CXXMatrix identityMatrix;
		FCXXCoord u;
		CXXMatrix uut;
		CXXMatrix s;
		float costheta, sintheta;
		
		u[0] = x; u[1] = y; u[2] = z; u[3] = 0.;
		u.normalise();
		
		for (i=0; i<4; i++){
			for (j=0; j<4; j++){
				uut[i][j] =  u[i]*u[j];
			}
		}
                // problems compiling CXX_DEGTORAD - don't understand
                // dyld: Symbol not found: __ZN8CXXCoordIfE12CXX_DEGTORADE
                // Referenced from: /Users/pemsley/autobuild/Darwin-phemius.local-gtk4/lib/libMoleculesToTrianglesCXXClasses.0.dylib
                // Expected in: flat namespace
                // in /Users/pemsley/autobuild/Darwin-phemius.local-gtk4/lib/libMoleculesToTrianglesCXXClasses.0.dylib
                // let's remove it.
		// costheta = cos(FCXXCoord ::CXX_DEGTORAD*w);
		// costheta = cos(FCXXCoord ::CXX_DEGTORAD*w);
		costheta = cos(w * 3.1415926/180.0);
		sintheta = sin(w * 3.1415926/180.0);

		s[0][0] =    0.; s[0][1] = -u[2]; s[0][2] =  u[1];
		s[1][0] =  u[2]; s[1][1] =    0.; s[1][2] = -u[0];
		s[2][0] = -u[1]; s[2][1] =  u[0]; s[2][2] =    0.;
		
		for (i=0, k=0; i<4; i++){
			for (j=0; j<4; j++, k++){
				result[i][j] = 0.;
			}
		}
		result[3][3] = 1.;
		
		for (i=0; i<3; i++){
			for (j=0; j<3; j++){
				result[i][j] = uut[i][j] + 
				costheta*(identityMatrix[i][j]-uut[i][j]) + 
				sintheta*s[i][j];
			}
		}
		return result;
	};	
	
	static CXXMatrix matrixWithScale(const double &x) {
		return matrixWithScale(x, x, x);
	};
	
	static CXXMatrix matrixWithScale(const double &x, const double &y, const double &z){
		CXXMatrix result;
		result[0][0] *= x; 
		result[1][1] *= y;
		result[2][2] *= z;
		return result;
	};
	
	static CXXMatrix matrixWithTranslation(const FCXXCoord &value){
        return matrixWithTranslation(value[0], value[1], value[2]);
    };
    
	static CXXMatrix matrixWithTranslation(const double &x, const double &y, const double &z){
		CXXMatrix result;
		result[3][0] = x; 
		result[3][1] = y;
		result[3][2] = z;
		return result;
	};
    
    static CXXMatrix lookAtMatrix(const FCXXCoord &rotation, const FCXXCoord &translation) {
        FCXXCoord minusTranslation = translation * -1.;
        CXXMatrix R(matrixWithGLRotation(rotation));
        CXXMatrix T(matrixWithTranslation(minusTranslation));
        CXXMatrix result(R*T);
        return result;
    };
    
	CXXMatrix transpose() const{
		CXXMatrix result;
		for (int i=0; i<4; i++){
			for (int j=0; j<4; j++){
				result[i][j] = rows[j][i];
			}
		}
		return result;
	};
    
	float scaleComponent(){
		double det = det3x3();
		// a negative determinant 
		double signbearer = copysign(1., det);
		det = (det<0.?-1.*det:det);
		double positiveScaleFactor = pow(det, 0.333333333);
		double result = copysign(positiveScaleFactor, signbearer);
		return result;
	};
    
	CXXMatrix inverseGLOperator(){
		
		CXXMatrix R, S, T, identityMatrix;
		decomposeTSR(T, S, R);
		
		CXXMatrix Tinv(identityMatrix);
		for (int i=0; i<3; i++) Tinv[3][i] = T[3][i]*-1.;
        
		CXXMatrix Rinv(R.transpose());
        
		CXXMatrix Sinv(identityMatrix * (1./S.scaleComponent()));
		Sinv[3][3] = 1.; 
        
		CXXMatrix inverseOperator(Rinv*(Sinv*Tinv));
        
		return inverseOperator;
	};
    
	FCXXCoord &operator [] (unsigned element)  {
		return rows[element];
	};
    
	const FCXXCoord operator [] (unsigned element) const {
		return rows[element];
	};
    
	FCXXCoord operator * (const FCXXCoord &target) const {
        //Actually returns target * *this 
        CXXMatrix myTranspose(transpose());
		return FCXXCoord (myTranspose[0]*target, 
                        myTranspose[1]*target, 
                        myTranspose[2]*target, 
                        myTranspose[3]*target);
	};
	
    CXXMatrix operator / (double invFactor){
		double factor = 1./invFactor;
		return (*this * factor);
    };
	
	CXXMatrix operator * (const CXXMatrix &target) const {
		CXXMatrix result;
		CXXMatrix transposed = target.transpose();
		for (int i=0; i<4; i++){
			for (int j=0; j<4; j++){
				result[i][j] = rows[i]*transposed[j];
			}
		}
		return result;
	};
	
	CXXMatrix operator * (const double &factor) const {
		CXXMatrix result;
		for (int i=0; i<4; i++){
			for (int j=0; j<4; j++){
				result[i][j] = rows[i][j]*factor;
			}
		}
		return result;
	};
	
	CXXMatrix operator - (const CXXMatrix &target) const {
		CXXMatrix result;
		for (int i=0; i<4; i++){
			for (int j=0; j<4; j++){
				result[i][j] = rows[i][j] - target[i][j];
			}
		}
		return result;
	};
	
	void decomposeTSR(CXXMatrix &T, CXXMatrix &S, CXXMatrix &R) {
		//Here I am going to analyse the current matrix and break it down as though 
		//it was a matrix that applies Rotation R, followed by Scale S, followed by translation T
		//i/e/ this matrix M = T*S*R
		
		CXXMatrix identityMatrix;
		T = identityMatrix;
		T[3] = rows[3];
		
		CXXMatrix SxR(*this);
		for (int i=0; i<3; i++) SxR[3][i] = 0.;
		SxR[3][3] = 1.;
        
		double scale = SxR.scaleComponent();
		S = identityMatrix * scale;	
        S[3][3] = 1.;
		
		R = SxR / scale;
		R[3][3] = 1.;		
	};
	
	void decomposeRST(CXXMatrix &R, CXXMatrix &S, CXXMatrix &T) {
		//Here I am going to analyse the current matrix and break it down as though 
		//it was a matrix that applies translation T, followed by Scale S, followed by Rotation R
		//i/e/ this matrix M = R*S*T
		
		CXXMatrix SxR(*this);
		for (int i=0; i<3; i++) SxR[3][i] = 0.;
		SxR[3][3] = 1.;
        
		double scale = SxR.scaleComponent();

        CXXMatrix identityMatrix;
		S = identityMatrix * scale;	
        S[3][3] = 1.;
		
		R = SxR / scale;
		R[3][3] = 1.;		

        CXXMatrix Sinv = identityMatrix / scale;
        Sinv[3][3] = 1.;
        
		CXXMatrix Rinv = R.transpose();
        CXXMatrix SinvxRinv = Sinv * Rinv;
        SinvxRinv[3][3] = 1.;
        
        T = SinvxRinv * (*this);
        //CXXMatrix RST = (R*(S*T));
	};
	
	void gLOperator(float *matrix){
		int i, k;
		for (i=0, k=0; i<4; i++){
			for (int j=0; j<4; j++){
				matrix[k++] = rows[i][j];
			}
		}
	};
	
	double det3x3(){
		double a, b, c, d, e, f, g, h, i;
		a = rows[0][0];
		b = rows[0][1];
		c = rows[0][2];
		d = rows[1][0];
		e = rows[1][1];
		f = rows[1][2];
		g = rows[2][0];
		h = rows[2][1];
		i = rows[2][2];
		double result = a*e*i + b*f*g + c*d*h - a*f*h - b*d*i - c*e*g;
        
		return result;
	};
	
	friend std::ostream &operator << ( std::ostream &out, const CXXMatrix &c ){
		out << "(";
		for (int i=0; i<4; i++) out << " " << c[i] << (i==3?")\n":"\n");
		return out;
	};
	
    FCXXCoord toAxisAngle() {
        CXXMatrix &m(*this);
        
        double angle,x,y,z; // variables for result
        double epsilon = 0.001; // margin to allow for rounding errors
        double epsilon2 = 0.01; // margin to distinguish between 0 and 180 degrees
        // optional check that input is pure rotation, 'isRotationMatrix' is defined at:
        // http://www.euclideanspace.com/maths/algebra/matrix/orthogonal/rotation/
        if ((fabs(m[0][1]-m[1][0])< epsilon)
            && (fabs(m[0][2]-m[2][0])< epsilon)
            && (fabs(m[1][2]-m[2][1])< epsilon)) {
            // singularity found
            // first check for identity matrix which must have +1 for all terms
            //  in leading diagonaland zero in other terms
            if ((fabs(m[0][1]+m[1][0]) < epsilon2)
                && (fabs(m[0][2]+m[2][0]) < epsilon2)
                && (fabs(m[1][2]+m[2][1]) < epsilon2)
                && (fabs(m[0][0]+m[1][1]+m[2][2]-3) < epsilon2)) {
                // this singularity is identity matrix so angle = 0
                return FCXXCoord (0.,1.,0.,0.);
            }
            // otherwise this singularity is angle = 180
            angle = M_PI;
            double xx = (m[0][0]+1)/2;
            double yy = (m[1][1]+1)/2;
            double zz = (m[2][2]+1)/2;
            double xy = (m[0][1]+m[1][0])/4;
            double xz = (m[0][2]+m[2][0])/4;
            double yz = (m[1][2]+m[2][1])/4;
            if ((xx > yy) && (xx > zz)) { // m[0][0] is the largest diagonal term
                if (xx< epsilon) {
                    x = 0;
                    y = 0.7071;
                    z = 0.7071;
                } else {
                    x = sqrt(xx);
                    y = xy/x;
                    z = xz/x;
                }
            } else if (yy > zz) { // m[1][1] is the largest diagonal term
                if (yy< epsilon) {
                    x = 0.7071;
                    y = 0;
                    z = 0.7071;
                } else {
                    y = sqrt(yy);
                    x = xy/y;
                    z = yz/y;
                }	
            } else { // m[2][2] is the largest diagonal term so base result on this
                if (zz< epsilon) {
                    x = 0.7071;
                    y = 0.7071;
                    z = 0;
                } else {
                    z = sqrt(zz);
                    x = xz/z;
                    y = yz/z;
                }
            }
			FCXXCoord aCoord(x, y, z, 0.);
			aCoord.normalise();
			x = aCoord[0];
			y = aCoord[1];
			z = aCoord[2];
            return FCXXCoord ((180./M_PI)*angle, x, y, z);        
        }
        // as we have reached here there are no singularities so we can handle normally
        double s = sqrt((m[2][1] - m[1][2])*(m[2][1] - m[1][2])
                        +(m[0][2] - m[2][0])*(m[0][2] - m[2][0])
                        +(m[1][0] - m[0][1])*(m[1][0] - m[0][1])); // used to normalise
        if (fabs(s) < 0.001) s=1;
		// prevent divide by zero, should not happen if matrix is orthogonal and should be
		// caught by singularity test above, but I've left it in just in case
        angle = acos(( m[0][0] + m[1][1] + m[2][2] - 1)/2);
        x = (m[2][1] - m[1][2])/s;
        y = (m[0][2] - m[2][0])/s;
        z = (m[1][0] - m[0][1])/s;

		FCXXCoord aCoord(x, y, z, 0.);
		aCoord.normalise();
		x = aCoord[0];
		y = aCoord[1];
		z = aCoord[2];

        return FCXXCoord ((180./M_PI)*angle,x,y,z);
    };
};

#endif
