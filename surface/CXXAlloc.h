/*
 *  CXXAlloc.h
 *  CXXSurface
 *
 *  Created by Martin Noble on 21/05/2009.
 *  Copyright 2009 LMB, Oxford University. All rights reserved.
 *
 */

#ifndef CXX_mot_CXXAlloc_included
#define CXX_mot_CXXAlloc_included

#include <stddef.h>
#include <stdlib.h>

#ifdef USENEDMALLOC
#include "nedmalloc.h"
#define cxxmalloc nedalloc::nedmalloc
#define cxxfree nedalloc::nedfree
#elif defined  USETCMALLOC
#include "google/tcmalloc.h"
#define cxxmalloc tc_malloc
#define cxxfree tc_free
#else
#define cxxmalloc malloc
#define cxxfree free
#endif


namespace CXX_old {
	inline void destruct(char*) {}
	inline void destruct(wchar_t*) {}
	template <typename T> inline void destruct(T* t) { t->~T(); }
    
	template <class T> class CXXAlloc;
	
    // specialize for void:
    template <> class CXXAlloc<void> {
    public: 
		typedef void*       pointer;
		typedef const void* const_pointer;
		// reference to void members are impossible.
		typedef void value_type;
		template <class U> struct rebind { typedef CXXAlloc<U>
			other; };
    };
	
    template <class T> class CXXAlloc {
    public:
		typedef size_t    size_type;
		typedef ptrdiff_t difference_type;
		typedef T*        pointer;
		typedef const T*  const_pointer;
		typedef T&        reference;
		typedef const T&  const_reference;
		typedef T         value_type;
		template <class U> struct rebind { 
            typedef CXXAlloc<U> other; 
        };
		
		CXXAlloc() throw() {};
		CXXAlloc(const CXXAlloc& other) throw(){
		};
		template <class U> CXXAlloc(const CXXAlloc<U>&) throw(){};
		~CXXAlloc() throw(){
		};
		
		pointer address(reference x) const {return &x;}
		const_pointer address(const_reference x) const{return &x;};
		
		pointer allocate(size_type size, CXXAlloc<void>::const_pointer hint = 0) 
		{
			return static_cast<pointer>(::cxxmalloc(size*sizeof(T)));
		};
		void deallocate(pointer p, size_type n) {
			if (p) ::cxxfree(p);
		};
		size_type max_size() const throw(){
			return size_t(-1) / sizeof(value_type);
		}
		void construct(pointer p, const T& val)
		{
			::new(static_cast<void*>(p)) T(val);
		}
		void destroy(pointer p) {
			CXX_old::destruct(p);
		};
    };

    template <typename T, typename U>
    inline bool operator==(const CXXAlloc<T>&, const CXXAlloc<U>){return true;}
    
    template <typename T, typename U>
    inline bool operator!=(const CXXAlloc<T>&, const CXXAlloc<U>){return false;}
    
};


#endif
