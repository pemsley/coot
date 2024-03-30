/*
 *  SecondaryColorScheme.mm
 *  Aesop
 *
 *  Created by Martin Noble on 23/02/2009.
 *  Copyright 2009 Dept. of Biochemistry, Oxford University. All rights reserved.
 *
 */

#include "SecondaryColorScheme.h"
#include "mmdb2/mmdb_manager.h"

SecondaryColorScheme *SecondaryColorScheme::defaultSecondaryScheme()
{
    SecondaryColorScheme *result = new SecondaryColorScheme();
    result->addPair(SecondaryColorPair(mmdb::SSE_None,FCXXCoord (1.,1.,1.,0.)));
    result->addPair(SecondaryColorPair(mmdb::SSE_Strand,FCXXCoord (1.,1.,0.,0.)));
    result->addPair(SecondaryColorPair(mmdb::SSE_Helix,FCXXCoord (1.,0.,1.,0.)));
    return result;
}
