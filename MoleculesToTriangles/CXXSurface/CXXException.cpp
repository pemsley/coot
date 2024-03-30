/*
 *  CXXException.cpp
 *  lpbSolver
 *
 *  Created by gruber on Sun Jul 11 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include "CXXException.h"


CXXException::CXXException(string m ) {
	
	message = m;
}
	
int CXXException::Report() {
	
	cout << message << "\n";
	return 0;
}
