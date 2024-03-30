/*
 * MoleculesToTriangles/CXXClasses/SecondaryColorScheme.h
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

#ifndef SecondaryColorScheme_h
#define SecondaryColorScheme_h

#include <vector>
#include "ColorScheme.h"

class SecondaryColorPair {
private:
    int secondary;
    FCXXCoord color;
public:
    SecondaryColorPair(){
    };
    SecondaryColorPair(int _secondary, const FCXXCoord &_color) :
    secondary(_secondary), color(_color) {};
    int getSecondary(){
        return secondary;
    };
    FCXXCoord getColor(){
        return color;
    };
};

class SecondaryColorScheme : public ColorScheme {
private:
    std::vector<SecondaryColorPair> pairs;
public:
    virtual FCXXCoord colorForAtom(mmdb::Atom* atom){
        FCXXCoord result = FCXXCoord (1.,1.,1.,0.);
        std::vector<SecondaryColorPair>::iterator pair = pairs.begin();
        while (pair != pairs.end()){
            if (atom->residue->SSE == pair->getSecondary()){
                result = pair->getColor();
            }
            pair++;
        }
        return result;
    };
    void addPair(SecondaryColorPair _pair){
        pairs.push_back(_pair);
    };
    static SecondaryColorScheme *defaultSecondaryScheme();
};

#endif
