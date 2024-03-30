/*
 *  CXXException.h
 *  lpbSolver
 *
 *  Created by gruber on Sun Jul 11 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>

#ifndef  __STRING_H
#include <string>
#endif

#define __CXXException__

using namespace std;

class CXXException {
	
private:
	
	string message;	
	
public:
	
	CXXException( string m );
	int Report();
	
};
