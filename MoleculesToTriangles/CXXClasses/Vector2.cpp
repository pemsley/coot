//---------------------------------------------------------------------------
//
//	File: Vector2.cpp
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

#include "Vector2.hpp"

//------------------------------------------------------------------------

//------------------------------------------------------------------------

template <typename Type>
	inline Vector2<Type>::Vector2()
	{
		x = 0;
		y = 0;
	} // Default Constructor

//------------------------------------------------------------------------

template <typename Type>
	inline Vector2<Type>::Vector2(const Type X, const Type Y)
	{
		x = X;
		y = Y;
	}// Constructor

//------------------------------------------------------------------------

template <typename Type>
	inline Vector2<Type>::Vector2(const Type *p)
	{
		if( p != NULL )
		{
			x = p[0];
			y = p[1];
		} // if
		else
		{
			x = 0;
			y = 0;
		} // else
	}// Constructor

//------------------------------------------------------------------------

template <typename Type>
	inline Vector2<Type>::Vector2(const Position2<Type> *p)
	{
		if( p != NULL )
		{
			x = p->x;
			y = p->y;
		} // if
		else
		{
			x = 0;
			y = 0;
		} // else
	}// Constructor

//------------------------------------------------------------------------

template <typename Type>
	inline Vector2<Type>::Vector2(const Vector2 &p)
	{
		x = p.x;
		y = p.y;
	}// Copy Constructor

//------------------------------------------------------------------------

template <typename Type>
	inline Vector2<Type> Vector2<Type>::operator+() 
	{ 
		return *this; 
	} // Vector2::operator+()

//------------------------------------------------------------------------

template <typename Type>
	inline Vector2<Type> Vector2<Type>::operator-() 
	{ 
		Vector2 p;

		p.x = -x;
		p.y = -y;

		return p;
	} // Vector2::operator-()

//------------------------------------------------------------------------

template <typename Type>
	inline Vector2<Type> Vector2<Type>::operator-(Type k)
	{
		Vector2 p;

		p.x = x - k;
		p.y = y - k;

		return p;
	} // Vector2::operator-

//------------------------------------------------------------------------

template <typename Type>
	inline Vector2<Type> Vector2<Type>::operator-(Vector2 &p)
	{
		Vector2 q;

		q.x = x - p.x;
		q.y = y - p.y;

		return q;
	} // Vector2::operator-

//------------------------------------------------------------------------

template <typename Type>
	inline Vector2<Type> Vector2<Type>::operator+(Type k)
	{
		Vector2 p;

		p.x = x + k;
		p.y = y + k;

		return p;
	} // Vector2::operator+

//------------------------------------------------------------------------

template <typename Type>
	inline Vector2<Type> Vector2<Type>::operator+(Vector2 &p)
	{
		Vector2 q;

		q.x = x + p.x;
		q.y = y + p.y;

		return q;
	} // Vector2::operator+

//------------------------------------------------------------------------

template <typename Type>
	inline Vector2<Type> Vector2<Type>::operator*(Type k)
	{
		Vector2 p;

		p.x = x * k;
		p.y = y * k;

		return p;
	} // Vector3::operator*

//------------------------------------------------------------------------

template <typename Type>
	inline Vector2<Type> Vector2<Type>::operator/(Type k)
	{
		Vector2 p;

		p.x = x / k;
		p.y = y / k;

		return p;
	} // Vector3::operator/

//------------------------------------------------------------------------

template <typename Type>
	inline Type Vector2<Type>::operator*(Vector2 &p)
	{
		Type m;

		m = x * p.x + y * p.y;

		return m;
	} // Vector2::operator*

//------------------------------------------------------------------------

template <typename Type>
	inline Vector2<Type> Vector2<Type>::operator-=(Type k)
	{
		x -= k;
		y -= k;

		return *this;
	} // Vector2::operator*=

//------------------------------------------------------------------------

template <typename Type>
	inline Vector2<Type> Vector2<Type>::operator+=(Type k)
	{
		x += k;
		y += k;

		return *this;
	} // Vector2::operator*=

//------------------------------------------------------------------------

template <typename Type>
	inline Vector2<Type> Vector2<Type>::operator*=(Type k)
	{
		x *= k;
		y *= k;

		return *this;
	} // Vector2::operator*=

//------------------------------------------------------------------------

template <typename Type>
	inline Vector2<Type> Vector2<Type>::operator/=(Type k)
	{
		if(s != 0)
		{
			x /= k;
			y /= k;
		} // if
		else 
		{
			x = 0;
			y = 0;
		} // else

		return *this;
	} // Vector2::operator/=

//------------------------------------------------------------------------

template <typename Type>
	inline Type Vector2<Type>::magnitude()
	{
		Type l = std::sqrt(x * x + y * y);

		return l;
	} // Vector2::magnitude

//------------------------------------------------------------------------

template <typename Type>
	inline Vector2<Type> Vector2<Type>::normalize()
	{
		Vector2<Type> w;
		
		Type L = this->magnitude();
		
		L = 1/L;
		
		w.x = L * x;
		w.y = L * y;
		
		return w;
	} // Vector2::normalize

//------------------------------------------------------------------------

template <typename Type>
	inline bool Vector2<Type>::renormalize(const Type e)
	{
		bool bRenormalized = false;
		
		Type L = this->magnitude();
		
		if( std::abs( L - 1 ) > e )
		{
			L = 1/L;
			
			x *= L;
			y *= L;
			
			bRenormalized = true;
		} // if
		
		return bRenormalized;
	} // Vector2::renormalize

//------------------------------------------------------------------------

template <typename Type>
	inline Type Vector2<Type>::cos(Vector2<Type> &p)
	{
		Type lp = p.magnitude();
		Type lq = this->magnitude();
		
		// interior scalar product c
		
		Type c = *this * p;
		
		lp = 1 / lp;
		lq = 1 / lq;
		
		// A = cos(a) = ( p * q ) / ( ||p|| ||q|| ) = c / ( lp * lq )
		
		Type A = lp * lq * c;
		
		return A;
	} // Vector2::cos

//------------------------------------------------------------------------

template <typename Type>
	inline Position2<Type> Vector2<Type>::position( )
	{
		Position2<Type> p;
		
		p.x = x;
		p.y = y;
		
		return p;
	} // Vector2::position

//------------------------------------------------------------------------

template <typename Type>
	inline Position2<Type> Vector2<Type>::diff(Vector2<Type> &v)
	{
		Position2<Type> p;
		
		p.x = std::abs(v.x - x);
		p.y = std::abs(v.y - y);
		
		return p;
	} // Vector2::diff

//------------------------------------------------------------------------

template <typename Type>
	inline void Vector2<Type>::swap() 
	{ 
		Type temp = x; 
		
		x = y; 
		y = temp; 
	} // Vector2<Type>::swap()

//------------------------------------------------------------------------

//------------------------------------------------------------------------

template class Position2<float>;
template class Position2<double>;

template class Vector2<float>;
template class Vector2<double>;

//------------------------------------------------------------------------

//------------------------------------------------------------------------
