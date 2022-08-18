//---------------------------------------------------------------------------
//
//	File: Vector3.cpp
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

#include "Vector2.hpp"
#include "Vector3.hpp"

//------------------------------------------------------------------------

//------------------------------------------------------------------------

template <typename Type>
	inline Vector3<Type>::Vector3()
	{
		x = 0;
		y = 0;
		z = 0;
	} // Default Constructor

//------------------------------------------------------------------------

template <typename Type>
	inline Vector3<Type>::Vector3(const Type X, const Type Y)
	{
		x = X;
		y = Y;
		z = 0;
	}// Constructor

//------------------------------------------------------------------------

template <typename Type>
	inline Vector3<Type>::Vector3(const Type X, const Type Y, const Type Z)
	{
		x = X;
		y = Y;
		z = Z;
	}// Constructor

//------------------------------------------------------------------------

template <typename Type>
	inline Vector3<Type>::Vector3(const Vector2<Type> &v)
	{
		x = v.x;
		y = v.y;
		z = 0;
	}// Constructor

//------------------------------------------------------------------------

template <typename Type>
	inline Vector3<Type>::Vector3(const Vector2<Type> *v)
	{
		if( v != NULL )
		{
			x = v->x;
			y = v->y;
		} // if
		else
		{
			x = 0;
			y = 0;
		} // else
		
		z = 0;
	}// Constructor

//------------------------------------------------------------------------

template <typename Type>
	inline Vector3<Type>::Vector3(const Type *v)
	{
		if( v != NULL )
		{
			x = v[0];
			y = v[1];
			z = v[2];
		} // if
		else
		{
			x = 0;
			y = 0;
			z = 0;
		} // else
	}// Constructor

//------------------------------------------------------------------------

template <typename Type>
	inline Vector3<Type>::Vector3(const Position2<Type> &p)
	{
		x = p.x;
		y = p.y;
		z = 0;
	}// Constructor

//------------------------------------------------------------------------

template <typename Type>
	inline Vector3<Type>::Vector3(const Position2<Type> *p)
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

		z = 0;
	}// Constructor

//------------------------------------------------------------------------

template <typename Type>
	inline Vector3<Type>::Vector3(const Position3<Type> &p)
	{
		x = p.x;
		y = p.y;
		z = p.z;
	}// Constructor

//------------------------------------------------------------------------

template <typename Type>
	inline Vector3<Type>::Vector3(const Position3<Type> *p)
	{
		if( p != NULL )
		{
			x = p->x;
			y = p->y;
			z = p->z;
		} // if
		else
		{
			x = 0;
			y = 0;
			z = 0;
		} // else
	}// Constructor

//------------------------------------------------------------------------

template <typename Type>
	inline Vector3<Type>::Vector3(const Vector3 &v)
	{
		x = v.x;
		y = v.y;
		z = v.z;
	}// Copy Constructor

//------------------------------------------------------------------------

template <typename Type>
	inline Vector3<Type> Vector3<Type>::operator+() 
	{ 
		return *this; 
	} // Vector2::operator+()

//------------------------------------------------------------------------

template <typename Type>
	inline Vector3<Type> Vector3<Type>::operator-() 
	{ 
		Vector3 p;

		p.x = -x;
		p.y = -y;
		p.z = -z;

		return p;
	} // Vector2::operator-()

//------------------------------------------------------------------------

template <typename Type>
	inline Vector3<Type> Vector3<Type>::operator-(Type t)
	{
		Vector3 w;

		w.x = x - t;
		w.y = y - t;
		w.z = z - t;

		return w;
	} // Vector3::operator-

//------------------------------------------------------------------------

template <typename Type>
	inline Vector3<Type> Vector3<Type>::operator-(Vector3 &v)
	{
		Vector3 w;

		w.x = x - v.x;
		w.y = y - v.y;
		w.z = z - v.z;

		return w;
	} // Vector3::operator-

//------------------------------------------------------------------------

template <typename Type>
	inline Vector3<Type> Vector3<Type>::operator+(Type t)
	{
		Vector3 w;

		w.x = x + t;
		w.y = y + t;
		w.z = z + t;

		return w;
	} // Vector3::operator+

//------------------------------------------------------------------------

template <typename Type>
	inline Vector3<Type> Vector3<Type>::operator+(Vector3 &v)
	{
		Vector3 w;

		w.x = x + v.x;
		w.y = y + v.y;
		w.z = z + v.z;

		return w;
	} // Vector3::operator+

//------------------------------------------------------------------------

template <typename Type>
	inline Vector3<Type> Vector3<Type>::operator*(Type s)
	{
		Vector3 w;

		w.x = x * s;
		w.y = y * s;
		w.z = z * s;

		return w;
	} // Vector3::operator*

//------------------------------------------------------------------------

template <typename Type>
	inline Vector3<Type> Vector3<Type>::operator/(Type s)
	{
		Vector3 w;

		w.x = x / s;
		w.y = y / s;
		w.z = z / s;

		return w;
	} // Vector3::operator/

//------------------------------------------------------------------------

template <typename Type>
	inline Type Vector3<Type>::operator*(Vector3 &v)
	{
		Type m = x * v.x + y * v.y + z * v.z;

		return m;
	} // Vector3::operator*

//------------------------------------------------------------------------
//
// Here, using the notation "^" for denoting cross products, comes from
// diffferential forms.  In differential forms, cross product is an
// exterior product denoted by "^".
//
//------------------------------------------------------------------------

template <typename Type>
	inline Vector3<Type> Vector3<Type>::operator^(Vector3<Type> &v)
	{
		Vector3 w;

		w.x = y * v.z - z * v.y;
		w.y = z * v.x - x * v.z;
		w.z = x * v.y - y * v.x;

		return w;
	} // Vector3::operator^

//------------------------------------------------------------------------

template <typename Type>
	inline Vector3<Type> Vector3<Type>::operator-=(Type t)
	{
		x -= t;
		y -= t;
		z -= t;

		return *this;
	} // Vector3::operator*=

//------------------------------------------------------------------------

template <typename Type>
	inline Vector3<Type> Vector3<Type>::operator+=(Type t)
	{
		x += t;
		y += t;
		z += t;

		return *this;
	} // Vector3::operator*=

//------------------------------------------------------------------------

template <typename Type>
	inline Vector3<Type> Vector3<Type>::operator*=(Type s)
	{
		x *= s;
		y *= s;
		z *= s;

		return *this;
	} // Vector3::operator*=

//------------------------------------------------------------------------

template <typename Type>
	inline Vector3<Type> Vector3<Type>::operator/=(Type s)
	{
		if(s != 0)
		{
			x /= s;
			y /= s;
			z /= s;
		} // if
		else 
		{
			x = 0;
			y = 0;
			z = 0;
		} // else

		return *this;
	} // Vector3::operator/=

//------------------------------------------------------------------------

template <typename Type>
	inline Type Vector3<Type>::magnitude()
	{
		Type l = std::sqrt(x * x + y * y + z * z);

		return l;
	} // Vector3::magnitude

//------------------------------------------------------------------------

template <typename Type>
	inline Type Vector3<Type>::IsInsideSphere( const Type radius )
	{
		Type L = 0;
		Type R = 0;
		Type Z = 0;
		
		L = x * x + y * y;
		R = radius * radius;
		
		if( L > R ) 
		{
			// On or outside the sphere.
			
			Z = 0;
		} 
		else
		{
			// Inside the sphere.

			Z = std::sqrt( R - L );
		} // else
		
		return Z;
	} // Vector3::isInsideSphere

//------------------------------------------------------------------------

template <typename Type>
	inline Type Vector3<Type>::cos(Vector3<Type> &v)
	{
		Type lu = this->magnitude();
		Type lv = v.magnitude();
		
		// interior scalar product c
		
		Type c = *this * v;
		
		lu = 1 / lu;
		lv = 1 / lv;
		
		// A = cos(a) = ( u * v ) / ( ||u|| ||v|| ) = c / ( lu * lv )
		
		Type A = lu * lv * c;
		
		return A;
	} // Vector3::cos

//------------------------------------------------------------------------

template <typename Type>
	inline Type Vector3<Type>::sin(Vector3<Type> &v)
	{
		Type lu = this->magnitude();
		Type lv = v.magnitude();
		
		// exterior vector product w
		
		Vector3 w = *this ^ v;
		
		// lw = || u ^ v ||
		
		Type lw = w.magnitude();
		
		lu = 1 / lu;
		lv = 1 / lv;
		
		// A = sin(a) = || u ^ v || / ( ||u|| ||v|| ) = lw / ( lu * lv )
		
		Type A = lu * lv * lw;
		
		return A;
	} // Vector3::sin

//------------------------------------------------------------------------

template <typename Type>
	inline Vector3<Type> Vector3<Type>::normalize()
	{
		Vector3<Type> w;
		
		Type L = this->magnitude();
		
		L = 1/L;
		
		w.x = L * x;
		w.y = L * y;
		w.z = L * z;
		
		return w;
	} // Vector3::normalize

//------------------------------------------------------------------------

template <typename Type>
	inline bool Vector3<Type>::renormalize(const Type e)
	{
		bool bRenormalized = false;
		
		Type L = this->magnitude();
		
		if( std::abs( L - 1 ) > e )
		{
			L = 1/L;
			
			x *= L;
			y *= L;
			z *= L;
			
			bRenormalized = true;
		} // if
		
		return bRenormalized;
	} // Vector3::renormalize

//------------------------------------------------------------------------

template <typename Type>
	inline Position3<Type> Vector3<Type>::position( )
	{
		Position3<Type> p;
		
		p.x = x;
		p.y = y;
		p.z = z;
		
		return p;
	} // Vector3::position

//------------------------------------------------------------------------

template <typename Type>
	inline Position3<Type> Vector3<Type>::diff(Vector3<Type> &v)
	{
		Position3<Type> p;
		
		p.x = std::abs(v.x - x);
		p.y = std::abs(v.y - y);
		p.z = std::abs(v.z - z);
		
		return p;
	} // Vector3::diff

//------------------------------------------------------------------------

template <typename Type>
	inline Vector3<Type> Vector3<Type>::normals(Vector3<Type> &v, Vector3<Type> &w)
	{
		Vector3<Type> n;
		Vector3<Type> dv;
		Vector3<Type> dw;

		dv = v - *this;
		dv = dv.normalize();
		
		dw = w - *this;
		dw = dw.normalize();
		
		n = dv ^ dw;
		n = n.normalize();

		return n;
	} // normals

//------------------------------------------------------------------------

template <typename Type>
	inline Vector3<Type> Vector3<Type>::normals(Position3<Type> &q, Position3<Type> &r)
	{
		Vector3<Type> v(q);
		Vector3<Type> w(r);
		
		Vector3<Type> n;
		Vector3<Type> dv;
		Vector3<Type> dw;

		dv = v - *this;
		dv = dv.normalize();
		
		dw = w - *this;
		dw = dw.normalize();
		
		n = dv ^ dw;
		n = n.normalize();

		return n;
	} // normals

//------------------------------------------------------------------------

//------------------------------------------------------------------------

template class Position3<float>;
template class Position3<double>;

template class Vector3<float>;
template class Vector3<double>;

//------------------------------------------------------------------------

//------------------------------------------------------------------------
