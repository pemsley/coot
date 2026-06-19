/*
 * MoleculesToTriangles/CXXClasses/DishyBase-gemmi.cc
 *
 * gemmi-native twin of DishyBase.cpp. Base atom names are UNPADDED to match
 * gemmi::Atom::name (the original used mmdb-padded names like " N1 ").
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#include "DishyBase-gemmi.hh"

void
coot::m2t::DishyBaseContainer_t::init() {
   cytidine_base_names = {"N1","C2","N3","C4","C5","C6","O2","N4"};
   uracil_base_names   = {"N1","C2","N3","C4","C5","C6","O2","O4"};
   adenine_base_names  = {"N9","C8","N7","C5","C4","N1","C2","N3","C6","N6"};
   guanine_base_names  = {"N9","C8","N7","C5","C4","N1","C2","N3","C6","O6","N2"};
   thymine_base_names  = {"N1","C2","N3","C4","C5","C6","O2","O4","C5M"};
}

std::vector<std::pair<int, int> > coot::m2t::DishyBase_t::bondingPattern = {
    std::pair<int,int>(0,1),
    std::pair<int,int>(1,2),
    std::pair<int,int>(2,3),
    std::pair<int,int>(3,4),
    std::pair<int,int>(4,0),
};
