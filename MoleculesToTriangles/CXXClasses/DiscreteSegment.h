/*
 * MoleculesToTriangles/CXXClasses/DiscreteSegment.h
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

#include "mmdb2/mmdb_manager.h"
#include "NRStuff.h"

#include <vector>
#include <tuple>
#include <map>
#include <memory>

class ColorScheme;
class ColorRule;

class DiscreteSegment {
private:
    std::vector<mmdb::Atom*> calphas;
    std::vector<FCXXCoord  >calphaCoords;
    std::vector<FCXXCoord  >normalOnes;
    std::vector<FCXXCoord  >normalTwos;
    std::vector<std::tuple<float, float, float>> anisoValues;
    CoordSpline coordSpline;
    CoordSpline normalOnesSpline;
    CoordSpline normalTwosSpline;
    CoordSpline anisoSpline;
    CoordSpline colorSpline;
public:
    ~DiscreteSegment(){
        calphas.clear();
        normalOnes.clear();
        normalTwos.clear();
        anisoValues.clear();
    }
    void addCalpha(mmdb::Atom* calpha){
        calphas.push_back(calpha);
        calphaCoords.push_back(FCXXCoord (calpha->x, calpha->y, calpha->z, 1.0f));
        anisoValues.push_back(std::make_tuple(1.0f, 1.0f, 1.0f));
    }
    void addCalpha(mmdb::Atom* calpha, float radius){
        calphas.push_back(calpha);
        calphaCoords.push_back(FCXXCoord (calpha->x, calpha->y, calpha->z, radius));
        anisoValues.push_back(std::make_tuple(1.0f, 1.0f, 1.0f));
    }
    void addCalpha(mmdb::Atom* calpha, float ax, float ay, float az){
        calphas.push_back(calpha);
        calphaCoords.push_back(FCXXCoord (calpha->x, calpha->y, calpha->z, 1.0f));
        anisoValues.push_back(std::make_tuple(ax, ay, az));
    }
    std::tuple<float, float, float> anisoFor(float xVal) {
        FCXXCoord aniso = anisoSpline.coordForXEquals(xVal);
        return std::make_tuple(aniso.x(), aniso.y(), aniso.z());
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
    void evaluateColors(std::shared_ptr<ColorScheme> colorScheme, std::map<std::shared_ptr<ColorRule>, int> &handles);
    FCXXCoord coordFor(float xVal);
    FCXXCoord normalOneFor(float xVal);
    FCXXCoord normalTwoFor(float xVal);
    FCXXCoord colorFor(float xVal);
};
