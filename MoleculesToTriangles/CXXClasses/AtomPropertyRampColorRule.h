/*
 * MoleculesToTriangles/CXXClasses/AtomPropertyRampColorRule.h
 *
 * Copyright 2009 by Martin Noble, University of Oxford
 * Author: Martin Noble
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */
#ifndef AtomPropertyRampColorRule_h
#define AtomPropertyRampColorRule_h

#include "ColorRule.h"
#include "NRStuff.h"
#include "CompoundSelection.h"
#include <memory>

class AtomPropertyRampColorRule : public ColorRule {
private:
//Modelled
        FCXXCoord startHSV, middleHSV, endHSV;
        float startValue, endValue;
        int rampType;
//For internal use
        CoordSpline spline;
        int ramp_points_size;
public:	
   enum RampType {
      BFactor, ResidueNumber
   };

   AtomPropertyRampColorRule() : startValue(1.), endValue(1000.0) {
      rank = 1.0;
      ramp_points_size = 6;
      type = AtomPropertyRamp;
      rampType = ResidueNumber;
      compoundSelection = std::shared_ptr<CompoundSelection>(new CompoundSelection("/*/*/*.*/*"));
      startHSV  = FCXXCoord (240., 0.4, 0.6, 1.);
      middleHSV = FCXXCoord (120., 0.4, 0.6, 1.);
      endHSV    = FCXXCoord (0.,   0.4, 0.6, 1.);
      updateSpline();
   };
        float getStartValue() const {
           return startValue;
        };
	void setStartValue(const float value){
	   startValue = value;
	   updateSpline();
	};
        float getEndValue() const {
           return endValue;
        };
	void setEndValue(float value){
	   endValue = value;
	   updateSpline();
	};
        FCXXCoord getStartRGB() {
           return hsvToRGB(startHSV);
        };
        void setStartRGB(FCXXCoord color){
	   startHSV = rgbToHSV(color);
           while (middleHSV[0] > startHSV[0]) middleHSV[0] -= 360.;
           while (endHSV[0] > middleHSV[0]) endHSV[0] -= 360.;
		updateSpline();
	};
        FCXXCoord getMiddleRGB() {
           return hsvToRGB(middleHSV);
        };
	void setMiddleRGB(FCXXCoord color){
	    middleHSV = rgbToHSV(color);
            //std::cout << "MiddleHSV" << middleHSV << std::endl;
            while (middleHSV[0] > startHSV[0]) middleHSV[0] -= 360.;
            while (endHSV[0] > middleHSV[0]) endHSV[0] -= 360.;
		updateSpline();
	};
        FCXXCoord getEndRGB() {
           return hsvToRGB(endHSV);
        };
	void setEndRGB(FCXXCoord color){
	   endHSV = rgbToHSV(color);
           while (middleHSV[0] > startHSV[0]) middleHSV[0] -= 360.;
           while (endHSV[0] > middleHSV[0]) endHSV[0] -= 360.;
		updateSpline();
	};
	int getRampType(){
		return rampType;
	};
	void setRampType(int type){
		rampType = type;
	};
        //virtual void prepareForSelectionInMMDB(int handle, mmdb::Manager *mmdb);
 
        void setNumberOfRampPoints(int _accu){
            ramp_points_size = _accu;
        }

        void updateSpline(){
		spline.clearSpline();
                for(int i=0;i<ramp_points_size;i++){
                    float f1 = (1.0-float(i)/ramp_points_size);
                    float f2 = (float(i)/ramp_points_size);
                    FCXXCoord hsv = middleHSV;
                    if(float(i)/ramp_points_size<0.5){
                        hsv[0] = f1*startHSV[0] + f2*middleHSV[0];
                    } else {
                        hsv[0] = f1*middleHSV[0] + f2*endHSV[0];
                    }
                    spline.addPair(f1*startValue+f2*endValue, hsv);
                }
		spline.calculateYDoublePrimes(1e30, 1e30);
	};

	FCXXCoord hsvToRGB(FCXXCoord hsv){
		FCXXCoord rgb;
                rgb[3] = hsv[3];
		int i;
		float f, p, q, t;
		if( hsv[1] == 0. ) {
			// achromatic (grey)
			rgb[0] = rgb[1] = rgb[2] = hsv[2];
			return rgb;
		}
		hsv[0] /= 60.f;			// sector 0 to 5
                //int floorhsv=floor(hsv[0]);
		i = fmod(floor( hsv[0] ), 6.f);
		f = hsv[0] - i;			// factorial part of h
		p = hsv[2] * ( 1 - hsv[1] );
		q = hsv[2] * ( 1 - hsv[1] * f );
		t = hsv[2] * ( 1 - (hsv[1] * ( 1 - f )) );
		switch( i ) {
			case 0:
				rgb[0] = hsv[2];
				rgb[1] = t;
				rgb[2] = p;
				break;
			case 1:
				rgb[0] = q;
				rgb[1]= hsv[2];
				rgb[2] = p;
				break;
			case 2:
				rgb[0] = p;
				rgb[1] = hsv[2];
				rgb[2] = t;
				break;
			case 3:
				rgb[0] = p;
				rgb[1] = q;
				rgb[2] = hsv[2];
				break;
			case 4:
				rgb[0] = t;
				rgb[1] = p;
				rgb[2] = hsv[2];
				break;
			default:		// case 5:
				rgb[0] = hsv[2];
				rgb[1] = p;
				rgb[2] = q;
				break;
		}
		return rgb;
	};
	FCXXCoord rgbToHSV(FCXXCoord rgb){
		FCXXCoord hsv;
        hsv[3] = rgb[3];
		float r = rgb[0];
		float g = rgb[1];
		float b = rgb[2];
		
		float min, max, delta;
		min = ((r<g?r:g)<b?(r<g?r:g):b);
		max = ((r>g?r:g)>b?(r>g?r:g):b);
		hsv[2] = max;				// v
		delta = max - min;
		if( max != 0 )
			hsv[1] = delta / max;		// s
		else {
			// r = g = b = 0		// s = 0, v is undefined
			hsv[1] = 0;
			hsv[0] = -1;
			return hsv;
		}
        if (fabs(delta) < 1e-6)
            hsv[0] = 0.;
		else if( r == max )
			hsv[0] = ( g - b ) / delta;		// between yellow & magenta
		else if( g == max )
			hsv[0] = 2 + ( b - r ) / delta;	// between cyan & yellow
		else
			hsv[0] = 4 + ( r - g ) / delta;	// between magenta & cyan
		hsv[0] *= 60.;				// degrees
		if( hsv[0] < 0 )
			hsv[0] += 360.;	
		return hsv;
	};
	
    static AtomPropertyRampColorRule *defaultRampRule();
	
	virtual FCXXCoord colorForAtom(const mmdb::Atom* atom){
		float value;
		if (rampType == ResidueNumber){
            const mmdb::Residue*residue = const_cast<mmdb::Atom*>(atom)->GetResidue();
			value = const_cast<mmdb::Residue*>(residue)->GetSeqNum();
			if (value < startValue) value = startValue;
			if (value > endValue) value = endValue;
		}
		else if (rampType == BFactor){
			value= atom->tempFactor;
			if (value < startValue) value = startValue;
			if (value > endValue) value = endValue;
		} else {
		  value = 0.5;
		}
        FCXXCoord rgb = colorForValue(value);
        
        return rgb;
    };
    
    FCXXCoord colorForValue(float value){
        if (value < startValue) value = startValue;
        if (value > endValue) value = endValue;
		FCXXCoord hsv = spline.coordForXEquals(value);
        
		while (hsv[0]<  0.)hsv[0] += 360.;
		while (hsv[0]>360.)hsv[0] -= 360.;
		hsv[1] = (hsv[1]<0.?0.:hsv[1]);
		hsv[1] = (hsv[1]>1.?1.:hsv[1]);
		hsv[2] = (hsv[2]<0.?0.:hsv[2]);
		hsv[2] = (hsv[2]>1.?1.:hsv[2]);
        hsv[3] = (hsv[3]<0.?0.:hsv[3]);
        hsv[3] = (hsv[3]>1.?1.:hsv[3]);
		FCXXCoord rgb = hsvToRGB(hsv);
        
        return rgb;
    };

};

#endif
