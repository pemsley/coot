//---------------------------------------------------------------------------
//
//	File: Quaternion.cpp
//
//  Abstract: C++ templates for Quaternion operators and methods
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
//  Copyright (c) 2008 Apple Inc., All rights reserved.
//
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------

#include "Quaternion.hpp"

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------

static const double kTwiceRadians2Degrees = 360.0/M_PI;
static const double kHalfDegrees2Radians  = M_PI/360.0;

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------

template <typename Type>
	inline Quaternion<Type>::Quaternion()
	{
		t = 0;
		x = 0;
		y = 0;
		z = 0;
	} // Default Constructor

//---------------------------------------------------------------------------

template <typename Type>
	inline Quaternion<Type>::Quaternion(const Type T, 
										const Type X, 
										const Type Y,
										const Type Z)
	{
		t = T;
		x = X;
		y = Y;
		z = Z;
	}// Constructor

//---------------------------------------------------------------------------

template <typename Type>
	inline Quaternion<Type>::Quaternion(const Type T, 
										const Type *v)
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
		
		t = T;
	}// Constructor

//---------------------------------------------------------------------------

template <typename Type>
	inline Quaternion<Type>::Quaternion(const Type T, 
										const Vector3<Type> &v)
	{
		x = v.x;
		y = v.y;
		z = v.z;
		t = T;
	}// Constructor

//---------------------------------------------------------------------------

template <typename Type>
	inline Quaternion<Type>::Quaternion(const Type T, 
										const Vector3<Type> *v)
	{
		if( v != NULL )
		{
			x = v->x;
			y = v->y;
			z = v->z;
		} // if
		else
		{
			x = 0;
			y = 0;
			z = 0;
		} // else
		
		t = T;
	}// Constructor

//---------------------------------------------------------------------------

template <typename Type>
	inline Quaternion<Type>::Quaternion(const Type T,
										const Position3<Type> &p)
	{
		x = p.x;
		y = p.y;
		z = p.z;
		t = T;
	}// Constructor

//---------------------------------------------------------------------------

template <typename Type>
	inline Quaternion<Type>::Quaternion(const Type T,
										const Position3<Type> *p)
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
		
		t = T;
	}// Constructor

//---------------------------------------------------------------------------

template <typename Type>
	inline Quaternion<Type>::Quaternion(const Quaternion &v)
	{
		x = v.x;
		y = v.y;
		z = v.z;
		t = v.t;
	}// Copy Constructor

//---------------------------------------------------------------------------

template <typename Type>
	inline Quaternion<Type> Quaternion<Type>::operator+() 
	{ 
		return *this; 
	} // Quaternion::operator+()

//---------------------------------------------------------------------------
//
// Quaternion conjugate
//
//---------------------------------------------------------------------------

template <typename Type>
	inline Quaternion<Type> Quaternion<Type>::operator-() 
	{ 
		Quaternion p;

		p.x = -x;
		p.y = -y;
		p.z = -z;
		p.t =  t;

		return p;
	} // Vector2::operator-()

//---------------------------------------------------------------------------

template <typename Type>
	inline Quaternion<Type> Quaternion<Type>::operator-(Type r)
	{
		Quaternion w;

		w.x = x - r;
		w.y = y - r;
		w.z = z - r;
		w.t = t - r;

		return w;
	} // Quaternion::operator-

//---------------------------------------------------------------------------

template <typename Type>
	inline Quaternion<Type> Quaternion<Type>::operator-(Quaternion &v)
	{
		Quaternion w;

		w.x = x - v.x;
		w.y = y - v.y;
		w.z = z - v.z;
		w.t = t - v.t;

		return w;
	} // Quaternion::operator-

//---------------------------------------------------------------------------

template <typename Type>
	inline Quaternion<Type> Quaternion<Type>::operator+(Type r)
	{
		Quaternion w;

		w.x = x + r;
		w.y = y + r;
		w.z = z + r;
		w.t = t + r;

		return w;
	} // Quaternion::operator+

//---------------------------------------------------------------------------

template <typename Type>
	inline Quaternion<Type> Quaternion<Type>::operator+(Quaternion &v)
	{
		Quaternion w;

		w.x = x + v.x;
		w.y = y + v.y;
		w.z = z + v.z;
		w.t = t + v.t;

		return w;
	} // Quaternion::operator+

//---------------------------------------------------------------------------

template <typename Type>
	inline Quaternion<Type> Quaternion<Type>::operator*(Type s)
	{
		Quaternion w;

		w.x = x * s;
		w.y = y * s;
		w.z = z * s;
		w.t = t * s;

		return w;
	} // Quaternion::operator*

//---------------------------------------------------------------------------

template <typename Type>
	inline Quaternion<Type> Quaternion<Type>::operator/(Type s)
	{
		Quaternion w;

		Type S = 1 / s;
		
		w.x = S * x;
		w.y = S * y;
		w.z = S * z;
		w.t = S * t;

		return w;
	} // Quaternion::operator/

//---------------------------------------------------------------------------

template <typename Type>
	inline Type Quaternion<Type>::operator*(Quaternion &v)
	{
		Type m = t * v.t + x * v.x + y * v.y + z * v.z;

		return m;
	} // Quaternion::operator*

//---------------------------------------------------------------------------
//
// This is better known as the Hamilton product
//
//---------------------------------------------------------------------------

template <typename Type>
	inline Quaternion<Type> Quaternion<Type>::operator^(Quaternion<Type> &v)
	{
		Quaternion w;
		
		w.x = v.y * z - v.z * y + v.t * x + v.x * t;
		w.y = v.z * x - v.x * z + v.t * y + v.y * t;
		w.z = v.x * y - v.y * x + v.t * z + v.z * t;
		w.t = v.t * t - v.x * x - v.y * y - v.z * z;
		
		return w;
	} // Quaternion::operator^

//---------------------------------------------------------------------------

template <typename Type>
	inline Quaternion<Type> Quaternion<Type>::operator-=(Type r)
	{
		x -= r;
		y -= r;
		z -= r;
		t -= r;
		
		return *this;
	} // Quaternion::operator*=

//---------------------------------------------------------------------------

template <typename Type>
	inline Quaternion<Type> Quaternion<Type>::operator+=(Type r)
	{
		x += r;
		y += r;
		z += r;
		t += r;

		return *this;
	} // Quaternion::operator*=

//---------------------------------------------------------------------------

template <typename Type>
	inline Quaternion<Type> Quaternion<Type>::operator*=(Type s)
	{
		x *= s;
		y *= s;
		z *= s;
		t *= s;

		return *this;
	} // Quaternion::operator*=

//---------------------------------------------------------------------------

template <typename Type>
	inline Quaternion<Type> Quaternion<Type>::operator/=(Type s)
	{
		if(s != 0)
		{
			x /= s;
			y /= s;
			z /= s;
			t /= s;
		} // if
		else 
		{
			x = 0;
			y = 0;
			z = 0;
			t = 0;
		} // else

		return *this;
	} // Quaternion::operator/=

//---------------------------------------------------------------------------

template <typename Type>
	inline Type Quaternion<Type>::magnitude()
	{
		Type l = std::sqrt( t * t + x * x + y * y + z * z );

		return l;
	} // Quaternion::magnitude

//---------------------------------------------------------------------------

template <typename Type>
	inline Quaternion<Type> Quaternion<Type>::normalize()
	{
		Quaternion<Type> w;
		
		Type L = this->magnitude();
		
		L = 1/L;
		
		w.x = L * x;
		w.y = L * y;
		w.z = L * z;
		w.t = L * t;
		
		return w;
	} // Quaternion::normalize

//---------------------------------------------------------------------------

template <typename Type>
	inline bool Quaternion<Type>::renormalize(const Type e)
	{
		bool bRenormalized = false;
		
		Type L = this->magnitude();

		if( std::abs( L - 1 ) > e )
		{
			L = 1/L;
			
			x *= L;
			y *= L;
			z *= L;
			t *= L;
			
			bRenormalized = true;
		} // if
		
		return bRenormalized;
	} // Quaternion::renormalize

//---------------------------------------------------------------------------
//
// An identity rotation is expressed as rotation by 0 about any axis.
// The "angle" term in a quaternion is really the cosine of the half-angle.
// So, if the cosine of the half-angle is one (or, 1.0 within our tolerance),
// then you have an identity rotation.
//
//---------------------------------------------------------------------------

template <typename Type>
	inline bool Quaternion<Type>::IsIdentityRotation(const Type e, Type *R)
	{
		bool bIdentityRotation = false;
		
        if( std::abs( std::abs(t) - 1 ) < e )
		{
			// Identity rotation.
			
			R[0] = 0;
			R[1] = 1;
			R[2] = 0;
			R[3] = 0;
			
			bIdentityRotation = true;
		} // if
		
		return bIdentityRotation;
	} // Quaternion::IsIdentityRotation

//---------------------------------------------------------------------------

template <typename Type>
	inline Position3<Type> Quaternion<Type>::position()
	{
		Position3<Type> p;
		
		p.x = x;
		p.y = y;
		p.z = z;
		
		return p;
	} // Quaternion::position

//---------------------------------------------------------------------------

template <typename Type>
	inline Vector3<Type> Quaternion<Type>::vector()
	{
		Vector3<Type> u(x,y,z);
		
		return u;
	} // Quaternion::vector

//---------------------------------------------------------------------------
//
// Convert from degrees to radians, get the half-angle.
//
//---------------------------------------------------------------------------
//
// Convert a GL-style rotation to a quaternion.  The GL rotation 
// looks like this:
//
//		{angle, x, y, z}, 
//
// the corresponding quaternion looks like this:
// 
//		{{v}, cos(angle/2)}, 
//
// where {v} is {x, y, z} / sin(angle/2).
//
//---------------------------------------------------------------------------

template <typename Type>
	inline void Quaternion<Type>::RotationToQuaternion( const Type *R )
	{		
		Type theta = kHalfDegrees2Radians * R[0]; // The half-angle
		
		Type A = std::sin( theta );
		
		x = A * R[1]; 
		y = A * R[2]; 
		z = A * R[3];
		
		t = std::cos( theta );
	} // Quaternion::RotationToQuaternion

//---------------------------------------------------------------------------
//
// Turn the quaternion back into an {angle, {axis}} rotation.
//
//---------------------------------------------------------------------------

template <typename Type>
	inline void Quaternion<Type>::QuaternionToRotation( Type *R )
	{
		Type theta = std::acos( t );
		
		Type A = 1/std::sin( theta );
		
		R[0] = kTwiceRadians2Degrees * theta;
		
		R[1] = A * x;
		R[2] = A * y;
		R[3] = A * z;
	} // Quaternion::QuaternionToRotation

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------

template class Quaternion<float>;
template class Quaternion<double>;

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
