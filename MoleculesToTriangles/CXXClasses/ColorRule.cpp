/*
 *  ColorRule.mm
 *  Aesop
 *
 *  Created by Martin Noble on 19/02/2009.
 *  Copyright 2009 Dept. of Biochemistry, Oxford University. All rights reserved.
 *
 */

#include "ColorRule.h"

#include "mmdb2/mmdb_manager.h"

bool ColorRule::compareRank(const std::shared_ptr<ColorRule> rule1, const std::shared_ptr<ColorRule> rule2) {
    return rule1->getRank() < rule2->getRank();
}


