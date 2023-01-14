/*
 *  DiscreetSegment.h
 *  MMDBRibbons
 *
 *  Created by Martin Noble on 17/07/2008.
 *  Copyright 2008 Dept. of Biochemistry, Oxford University. All rights reserved.
 *
 */

#include "mmdb2/mmdb_manager.h"
#include "NRStuff.h"

#include <vector>

class DiscreteSegment {
private:
    std::vector<mmdb::Atom*> calphas;
    std::vector<FCXXCoord  >calphaCoords;
    std::vector<FCXXCoord  >normalOnes;
    std::vector<FCXXCoord  >normalTwos;
    CoordSpline coordSpline;
    CoordSpline normalOnesSpline;
    CoordSpline normalTwosSpline;
public:
    ~DiscreteSegment(){
        calphas.clear();
        normalOnes.clear();
        normalTwos.clear();
    }
    void addCalpha(mmdb::Atom* calpha){ 
        calphas.push_back(calpha);
        calphaCoords.push_back(FCXXCoord (calpha->x, calpha->y, calpha->z));
    }
    FCXXCoord operator [] (int i) {
        return FCXXCoord (calphas[i]->x, calphas[i]->y, calphas[i]->z);
    }
    int nCalphas() {
        return int(calphas.size());
    }
    mmdb::Atom* calpha(int iCalpha) {
        return calphas[iCalpha];
    }    
    void evaluateNormals();
    void smoothBetas();
    void evaluateSplines();
    FCXXCoord coordFor(float xVal);
    FCXXCoord normalOneFor(float xVal);
    FCXXCoord normalTwoFor(float xVal);
};
