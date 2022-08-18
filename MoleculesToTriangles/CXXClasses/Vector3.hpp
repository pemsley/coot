//---------------------------------------------------------------------------
//
//	File: Vector3.hpp
//
//  Abstract: C++ templates for 3-Vector operators and methods
// 			 
//  Disclaimer: IMPORTANT:  This Apple software is supplied to you by
//  Inc. ("Apple") in consideration of your agreement to the following terms, 
//  and your use, installation, modification or redistribution of this Apple 
//  software constitutes acceptance of these terms.  If you do not agree with 
//  these terms, please do not use, install, modify or redistribute this 
//  Apple software.
//  
//  In consideration of your agreement to abide by the following terms, and
//  subject to these terms, Apple grants you a personal, non-exclusive
//  license, under Apple's copyrights in this original Apple software (the
//  "Apple Software"), to use, reproduce, modify and redistribute the Apple
//  Software, with or without modifications, in source and/or binary forms;
//  provided that if you redistribute the Apple Software in its entirety and
//  without modifications, you must retain this notice and the following
//  text and disclaimers in all such redistributions of the Apple Software. 
//  Neither the name, trademarks, service marks or logos of Apple Inc. may 
//  be used to endorse or promote products derived from the Apple Software 
//  without specific prior written permission from Apple.  Except as 
//  expressly stated in this notice, no other rights or licenses, express
//  or implied, are granted by Apple herein, including but not limited to
//  any patent rights that may be infringed by your derivative works or by
//  other works in which the Apple Software may be incorporated.
//  
//  The Apple Software is provided by Apple on an "AS IS" basis.  APPLE
//  MAKES NO WARRANTIES, EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION
//  THE IMPLIED WARRANTIES OF NON-INFRINGEMENT, MERCHANTABILITY AND FITNESS
//  FOR A PARTICULAR PURPOSE, REGARDING THE APPLE SOFTWARE OR ITS USE AND
//  OPERATION ALONE OR IN COMBINATION WITH YOUR PRODUCTS.
//  
//  IN NO EVENT SHALL APPLE BE LIABLE FOR ANY SPECIAL, INDIRECT, INCIDENTAL
//  OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
//  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
//  INTERRUPTION) ARISING IN ANY WAY OUT OF THE USE, REPRODUCTION,
//  MODIFICATION AND/OR DISTRIBUTION OF THE APPLE SOFTWARE, HOWEVER CAUSED
//  AND WHETHER UNDER THEORY OF CONTRACT, TORT (INCLUDING NEGLIGENCE),
//  STRICT LIABILITY OR OTHERWISE, EVEN IF APPLE HAS BEEN ADVISED OF THE
//  POSSIBILITY OF SUCH DAMAGE.
// 
//  Copyright (c) 2007-2008 Apple Inc., All rights reserved.
//
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------

#ifndef _VECTOR3_HPP_
#define _VECTOR3_HPP_

#ifdef __cplusplus

#include <cmath>
#include <iostream>

//------------------------------------------------------------------------

#include "Vector2.hpp"

//------------------------------------------------------------------------

//------------------------------------------------------------------------

template <typename Type>
	struct Position3 
	{
		union 
		{
			Type V[3];
			
			struct 
			{
				Type x;
				Type y;
				Type z;
			}; // union
		}; // struct
	};

//------------------------------------------------------------------------

template <typename Type>
	class Vector3
	{
		public:
		
			Vector3();
			Vector3(const Type X, const Type Y);
			Vector3(const Type X, const Type Y, const Type Z);
			Vector3(const Type *v);
			Vector3(const Position2<Type> &p);
			Vector3(const Position2<Type> *p);
			Vector3(const Position3<Type> &p);
			Vector3(const Position3<Type> *p);
			Vector3(const Vector2<Type> &v);
			Vector3(const Vector2<Type> *v);

			Vector3(const Vector3 &v);

			Vector3 operator-(Vector3 &v);
			Vector3 operator+(Vector3 &v);
			Type    operator*(Vector3 &v);  // Interior dot product
			Vector3 operator^(Vector3 &v);  // Exterior cross product

			Vector3 operator-(Type t);
			Vector3 operator+(Type t);
			Vector3 operator*(Type s);
			Vector3 operator/(Type s);
			
			Vector3 operator+=(Type t);
			Vector3 operator-=(Type t);
			Vector3 operator*=(Type s);
			Vector3 operator/=(Type s);
			
			Vector3 operator+();
			Vector3 operator-();

			Type magnitude();
			Type cos(Vector3 &v);
			Type sin(Vector3 &v);
		
			Vector3 normalize();
			Vector3 normals(Vector3 &v, Vector3 &w);
			Vector3 normals(Position3<Type> &q, Position3<Type> &r);
			
			Position3<Type> diff(Vector3 &v);
			Position3<Type> position();
			
			bool renormalize(const Type e);
		
			Type IsInsideSphere(const Type radius);
		
		public:
			
			Type x;
			Type y;
			Type z;
	}; // class Vector3

//------------------------------------------------------------------------

//------------------------------------------------------------------------

typedef Position3<float>  FPosition3;
typedef Position3<double> DPosition3;

typedef Vector3<float>   FVector3;
typedef Vector3<double>  DVector3;

//------------------------------------------------------------------------

//------------------------------------------------------------------------


#endif

#endif