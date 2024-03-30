#include <iostream>
#include "CXXCoord.h"

using namespace std;

class CXXMatrix4x4{
private:
	CXXCoord<CXXCoord_ftype>columns[4];
	
public:
	CXXMatrix4x4(){
		columns[0] = CXXCoord<CXXCoord_ftype>(1.,0.,0.,0.);
		columns[1] = CXXCoord<CXXCoord_ftype>(0.,1.,0.,0.);
		columns[2] = CXXCoord<CXXCoord_ftype>(0.,0.,1.,0.);
		columns[3] = CXXCoord<CXXCoord_ftype>(0.,0.,0.,1.);
	};
	
	CXXMatrix4x4(const CXXCoord<CXXCoord_ftype>&c0, const CXXCoord<CXXCoord_ftype>&c1, const CXXCoord<CXXCoord_ftype>&c2, const CXXCoord<CXXCoord_ftype>&c3) {
		columns[0] = c0;
		columns[1] = c1;
		columns[2] = c2;
		columns[3] = c3;
	};
	
	const CXXCoord<CXXCoord_ftype>&operator [] (const int iColumn) const {
		return columns[iColumn];
	};
	
	CXXCoord<CXXCoord_ftype>operator * (const CXXCoord<CXXCoord_ftype>&coord) const {
		double result[] = {0., 0., 0., 0.};
		for (int i=0; i<4; i++){
			for (int j=0; j<4; j++){
				result[i] += (columns[j][i] * coord[j]);
			}
		}
		return (CXXCoord<CXXCoord_ftype>(result));
	};
	
	CXXMatrix4x4 operator * (const CXXMatrix4x4 &matrix) const {
		double result[][4] = {{0., 0., 0., 0.},{0., 0., 0., 0.},{0., 0., 0., 0.},{0., 0., 0., 0.}};
		for (int i=0; i<4; i++){
			for (int j=0; j<4; j++){
				for (int k=0; k<4; k++){
					result[j][i] += columns[k][i] * matrix[j][k];
				}
			}
		}
		return CXXMatrix4x4 (CXXCoord<CXXCoord_ftype>(result[0]),
							 CXXCoord<CXXCoord_ftype>(result[1]),
							 CXXCoord<CXXCoord_ftype>(result[2]),
							 CXXCoord<CXXCoord_ftype>(result[3]));
	};
	
	inline CXXMatrix4x4  invert() {
		float detx;
		float xx[9],zz[9];
		//Treat as 3x3 matrix
		double result[][4] = {{0., 0., 0., 0.},{0., 0., 0., 0.},{0., 0., 0., 0.},{0., 0., 0., 1.}};
		int i, j, k;
		
		k=0;
		for (i=0; i<3; i++)
			for (j=0; j<3; j++)
				xx[k++] = columns[i][j];
		
		zz[0] = xx[4]*xx[8] - xx[5]*xx[7];
		zz[1] = xx[2]*xx[7] - xx[1]*xx[8];
		zz[2] = xx[1]*xx[5] - xx[2]*xx[4];
		zz[3] = xx[5]*xx[6] - xx[3]*xx[8];
		zz[4] = xx[0]*xx[8] - xx[2]*xx[6];
		zz[5] = xx[2]*xx[3] - xx[0]*xx[5];
		zz[6] = xx[3]*xx[7] - xx[4]*xx[6];
		zz[7] = xx[1]*xx[6] - xx[0]*xx[7];
		zz[8] = xx[0]*xx[4] - xx[1]*xx[3];
		detx = zz[0]*xx[0] + zz[1]*xx[3] + zz[2]*xx[6];
		
		if (detx == 0.) {
			cout << "Singular matrix :Returning identity\n";
			return CXXMatrix4x4();
		}
		
		k = 0;
		for (i=0; i<3; i++)
			for (j=0; j<3; j++)
				result[i][j] = zz[k++]/detx;
		
		return CXXMatrix4x4(CXXCoord<CXXCoord_ftype>(result[0]),
							CXXCoord<CXXCoord_ftype>(result[1]),
							CXXCoord<CXXCoord_ftype>(result[2]),
							CXXCoord<CXXCoord_ftype>(result[3]));
	};
	
	CXXMatrix4x4 transpose(void) const{
		double result[][4] = {{0., 0., 0., 0.},{0., 0., 0., 0.},{0., 0., 0., 0.},{0., 0., 0., 0.}};
		for (int i=0; i<4; i++){
			for (int j=0; j<4; j++){
				result[i][j] = columns[j][i];
			}
		}
		return CXXMatrix4x4(CXXCoord<CXXCoord_ftype>(result[0]),
							CXXCoord<CXXCoord_ftype>(result[1]),
							CXXCoord<CXXCoord_ftype>(result[2]),
							CXXCoord<CXXCoord_ftype>(result[3]));
	};
	
	void printf(){
		for (int i=0; i<4; i++){
			for (int j=0; j<4; j++){
				cout <<columns[j][i] <<" ";
			}
			cout << "\n";
		}
	};
};


