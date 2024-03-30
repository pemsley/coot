/*
 *  SecondaryColorScheme.h
 *  Aesop
 *
 *  Created by Martin Noble on 23/02/2009.
 *  Copyright 2009 Dept. of Biochemistry, Oxford University. All rights reserved.
 *
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
