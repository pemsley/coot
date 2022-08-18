//---------------------------------------------------------------------------
//
//	File: Vector2.hpp
//
//  Abstract: C++ templates for 2-Vector operators and methods
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

#ifndef _VECTOR2_HPP_
#define _VECTOR2_HPP_

#ifdef __cplusplus

#include <cmath>
#include <iostream>

//------------------------------------------------------------------------

//------------------------------------------------------------------------

template <typename Type>
	struct Position2 
	{
		union 
		{
			Type V[2];
			
			struct 
			{
				union { Type x, u, s; };
				union { Type y, v, t; };
			}; // union
		}; // struct
	};

//------------------------------------------------------------------------

template <typename Type>
	class Vector2
	{
		public:
		
			Vector2();
			Vector2(const Type X, const Type Y);
			Vector2(const Type *p);
			Vector2(const Position2<Type> *p);

			Vector2(const Vector2 &r);

			Vector2 operator-(Vector2 &p);
			Vector2 operator+(Vector2 &p);
			Type    operator*(Vector2 &p);  // Interior dot product
			
			Vector2 operator+(Type k);
			Vector2 operator-(Type k);
			Vector2 operator*(Type k);
			Vector2 operator/(Type k);

			Vector2 operator+=(Type k);
			Vector2 operator-=(Type k);
			Vector2 operator*=(Type k);
			Vector2 operator/=(Type k);
			
			Vector2 operator+();
			Vector2 operator-();

			Type magnitude();
		
			Type cos(Vector2 &p);
			
			Vector2 normalize();
			
			bool renormalize(const Type e);
			
			void swap();

			Position2<Type> diff(Vector2 &v);
			Position2<Type> position( );
			
		public:
		
			union { Type x, u, s; };
			union { Type y, v, t; };
	}; // class Vector2

//------------------------------------------------------------------------

//------------------------------------------------------------------------

typedef Position2<float>  FPosition2;
typedef Position2<double> DPosition2;

typedef Vector2<float> FVector2;
typedef Vector2<double> DVector2;

//------------------------------------------------------------------------

//------------------------------------------------------------------------


#endif

#endif
