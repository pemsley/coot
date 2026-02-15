/*
 * MoleculesToTriangles/CXXClasses/DiscreteSegment.cpp
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

#include <string>
#include "DiscreteSegment.h"

extern "C" {
#include "lfit.h"
}

void DiscreteSegment::evaluateNormals()
{
    normalOnes.resize(calphas.size());
    normalTwos.resize(calphas.size());
    FCXXCoord xAxis(1.,0.,0.);
    FCXXCoord yAxis(0.,1.,0.);
    normalOnes[0] = xAxis;
    normalOnes[normalOnes.size()-1] = xAxis;
    normalTwos[0] = yAxis;
    normalTwos[normalOnes.size()-1] = yAxis;
    if (normalOnes.size()<3) return;
    
    for (int i=1; i<(calphas.size()-1); i++){
        FCXXCoord a (calphaCoords[i-1]);
        FCXXCoord b (calphaCoords[i]);
        FCXXCoord c (calphaCoords[i+1]);
        FCXXCoord ab = b-a;
        FCXXCoord cb = b-c;
        FCXXCoord ac = c-a;
        FCXXCoord normal = ab^cb;
        if (normal*normalOnes[i-1]<0.){
            normal *= -1.;
        }
        normal.normalise();
        normalOnes[i] = normal;
        FCXXCoord normal2 = normal^ac;
        normal2.normalise();
        normalTwos[i] = normal2;
    }
    normalOnes[0] = normalOnes[1];
    normalOnes[normalOnes.size()-1] = normalOnes[normalOnes.size()-2];
    normalTwos[0] = normalTwos[1];
    normalTwos[normalTwos.size()-1] = normalTwos[normalTwos.size()-2];
}

void DiscreteSegment::smoothBetas()
{
    float *xes = new float[calphas.size()];
    
    float *yes[3];
    float *n1s[3];
    float *n2s[3];
    
    for (int i=0; i<3; i++){
        yes[i] = new float[calphas.size()];
        n1s[i] = new float[calphas.size()];
        n2s[i] = new float[calphas.size()];
    }

    int *fixeds = new int[calphas.size()];
    
    for (int iCalpha=0; iCalpha<calphas.size(); iCalpha++){
        mmdb::Atom* atom = calphas[iCalpha];
        int nInStrand = 0;
        while ((atom->GetResidue()->SSE == mmdb::SSE_Strand ||
                atom->GetResidue()->SSE == mmdb::SSE_Bulge)
                && iCalpha+nInStrand < calphas.size()){
            fixeds[nInStrand] = 0;
            xes[nInStrand] = nInStrand;
            
            for (int i=0; i<3; i++){
                yes[i][nInStrand] = calphaCoords[nInStrand + iCalpha][i];
                n1s[i][nInStrand] = normalOnes[nInStrand + iCalpha][i];
                n2s[i][nInStrand] = normalTwos[nInStrand + iCalpha][i];
            }
            nInStrand++;
            if (iCalpha+nInStrand < calphas.size()){
                atom = calphas[iCalpha+nInStrand];
            }
        }
        if (nInStrand >= 2){
            fixeds[0] = fixeds[nInStrand-1] = 1;
            for (int i=0; i<3; i++){
                int order = nInStrand - 1;
                if (order >3) order = 3; 
                ForcePolynomial(order, xes, yes[i], fixeds, nInStrand, 0, nInStrand-1);
                ForcePolynomial(order, xes, n1s[i], fixeds, nInStrand, 0, nInStrand-1);
                ForcePolynomial(order, xes, n2s[i], fixeds, nInStrand, 0, nInStrand-1);
            }
            for (int i=0; i<nInStrand;i++){
                for (int j=0; j<3; j++){
                    calphaCoords[iCalpha + i][j] = yes[j][i];
                    normalOnes[iCalpha + i][j] = n1s[j][i];
                    normalTwos[iCalpha + i][j] = n2s[j][i];
                }
            }
        }
        iCalpha+= nInStrand;
    }    

    for (int i=0; i<3; i++){
        delete [] yes[i];    
        delete [] n1s[i];    
        delete [] n2s[i];    
    }
    
    delete [] xes;
    delete [] fixeds;
}

void DiscreteSegment::evaluateSplines(int samplesPerSegment)
{
   // std::cout << "evaluateSpines() size " << calphas.size() << std::endl;
    for (int i=0; i<calphas.size(); i++){
       // std::cout << "evaluateSpines() i = " << i << std::endl;
        coordSpline.addPair((float)i, calphaCoords[i]);
        normalOnesSpline.addPair((float)i, normalOnes[i]);
        normalTwosSpline.addPair((float)i, normalTwos[i]);
    }
    // std::cout << "evaluateSpines() A" << std::endl;
    coordSpline.calculateYDoublePrimes(1e30f, 1e30f, samplesPerSegment);
    // std::cout << "evaluateSpines() B" << std::endl;
    normalOnesSpline.calculateYDoublePrimes(1e30f, 1e30f, samplesPerSegment);
    // std::cout << "evaluateSpines() C" << std::endl;
    normalTwosSpline.calculateYDoublePrimes(1e30f, 1e30f, samplesPerSegment);
    // std::cout << "evaluateSpines() D" << std::endl;
}

FCXXCoord DiscreteSegment::coordFor(float xVal)
{
    return coordSpline.coordForXEquals(xVal);
}
FCXXCoord DiscreteSegment::normalOneFor(float xVal)
{
    return normalOnesSpline.coordForXEquals(xVal);
}
FCXXCoord DiscreteSegment::normalTwoFor(float xVal)
{
    return normalTwosSpline.coordForXEquals(xVal);
}
