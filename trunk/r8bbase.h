//$ nobt
//$ nocpp

/**
 * \file r8bbase.h
 * \brief The "base" inclusion file with basic classes and functions.
 *
 * This is the "base" inclusion file for the "r8brain-free-src" sample rate
 * converter. This inclusion file contains implementations of several small
 * utility classes and functions used by the library.
 *
 * r8brain-free-src Copyright (c) 2013 Aleksey Vaneev
 * See the "License.txt" file for license.
 *
 * \mainpage
 * \section intro_sec Introduction
 *
 * Open source high-quality professional audio sample rate converter (SRC)
 * library. Features routines for realtime SRC (non-realtime conversion is
 * also possible), both up- and down-sampling, from/to any sample rate,
 * including non-integer sample rates. SRC routines were implemented in
 * multi-platform C++ code, and have a high level of optimality. The user can
 * select the transition band/steepness of the low-pass (reconstruction)
 * filter, expressed as a percentage of the full spectral bandwidth, and the
 * desired stop-band attenuation in decibel.
 *
 * The structure of this library's objects is such that they can be frequently
 * created and destroyed in large applications with minimal performance impact
 * due to a high level of reusability of its most "initialization-expensive"
 * objects: the fast Fourier transform and filter response objects.
 *
 * The algorithm at first produces 2X to 3X oversampled signal and then
 * performs fractional delay interpolation using a bank of very short (8 taps)
 * cross-faded FIR filters. This puts the algorithm into the league of the
 * fastest among the most precise SRC algorithms. The more precise alternative
 * being only the whole number-factored SRC, which is slower.
 *
 * \section license License
 *
 * The MIT License (MIT)
 * 
 * r8brain-free-src Copyright (c) 2013 Aleksey Vaneev
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 * 
 * Please credit the creator of this library in your documentation in the
 * following way: "Sample rate converter designed by Aleksey Vaneev"
 *
 * @version 0.1
 */

#ifndef R8BBASE_INCLUDED
#define R8BBASE_INCLUDED

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "r8bconf.h"

#if defined( R8B_WIN )
	#include <windows.h>
#elif defined( R8B_MAC )
	#include <libkern/OSAtomic.h>
#elif defined( R8B_LNX )
	#include <pthread.h>
#endif

/**
 * \brief The "r8brain-free-src" library namespace.
 *
 * The "r8brain-free-src" sample rate converter library namespace.
 */

namespace r8b {

#if !defined( M_PI )
	/**
	 * The macro equals to "pi" constant, fits 53-bit floating point mantissa.
	 */

	#define M_PI 3.14159265358979324
#endif // M_PI

#if !defined( M_2PI )
	/**
	 * The M_2PI macro equals to "2 * pi" constant, fits 53-bit floating point
	 * mantissa.
	 */

	#define M_2PI 6.28318530717958648
#endif // M_2PI

/**
 * A special macro that defines empty copy-constructor and copy operator with
 * the "private:" prefix. This macro should be used in classes that cannot be
 * copied in a standard C++ way.
 *
 * This macro do not need to be defined in classes derived from a class where
 * such macro was already used.
 *
 * @param ClassName The name of the class which uses this macro.
 */

#define R8BNOCTOR( ClassName ) \
	private: \
		ClassName( const ClassName& ) { } \
		ClassName& operator = ( const ClassName& ) { return( *this ); }

/**
 * \brief The default base class for objects created on heap.
 *
 * Class that implements "new" and "delete" operators that use standard
 * malloc() and free() functions.
 */

class CStdClassAllocator
{
public:
	/**
	 * @param n The size of the object, in bytes.
	 * @param p Pointer to object's pre-allocated memory block.
	 * @return Pointer to object.
	 */

	void* operator new( size_t n, void* p )
	{
		return( p );
	}

	/**
	 * @param n The size of the object, in bytes.
	 * @return Pointer to the allocated memory block for the object.
	 */

	void* operator new( size_t n )
	{
		return( :: malloc( n ));
	}

	/**
	 * @param n The size of the object, in bytes.
	 * @return Pointer to the allocated memory block for the object.
	 */

	void* operator new[]( size_t n )
	{
		return( :: malloc( n ));
	}

	/**
	 * Operator frees a previously allocated memory block for the object.
	 *
	 * @param p Pointer to the allocated memory block for the object.
	 */

	void operator delete( void* p )
	{
		:: free( p );
	}

	/**
	 * Operator frees a previously allocated memory block for the object.
	 *
	 * @param p Pointer to the allocated memory block for the object.
	 */

	void operator delete[]( void* p )
	{
		:: free( p );
	}
};

/**
 * \brief The default base class for objects that allocate blocks of memory.
 *
 * Memory buffer allocator that uses "stdlib" standard memory functions.
 */

class CStdMemAllocator : public CStdClassAllocator
{
public:
	/**
	 * Function allocates memory block.
	 *
	 * @param Size The size of the block, in bytes.
	 * @result The pointer to the allocated block.
	 */

	static void* allocmem( const size_t Size )
	{
		return( :: malloc( Size ));
	}

	/**
	 * Function reallocates a previously allocated memory block.
	 *
	 * @param p Pointer to the allocated block, can be NULL.
	 * @param Size The new size of the block, in bytes.
	 * @result The pointer to the (re)allocated block.
	 */

	static void* reallocmem( void* p, const size_t Size )
	{
		return( :: realloc( p, Size ));
	}

	/**
	 * Function frees a previously allocated memory block.
	 *
	 * @param p Pointer to the allocated block, can be NULL.
	 */

	static void freemem( void* p )
	{
		:: free( p );
	}
};

/**
 * \brief Templated memory buffer class for element buffers of fixed size.
 *
 * Fixed memory buffer object. Supports allocation of a fixed amount of
 * memory. Does not store buffer's capacity - the user should know the actual
 * size of the buffer. Does not feature "internal" storage, memory is always
 * allocated via the R8B_MEMALLOCCLASS class's functions.
 *
 * This class manages memory space only - it does not perform element class
 * construction operations.
 *
 * @param T The class of the stored elements (e.g. "double").
 */

template< class T >
class CFixedBuffer : public R8B_MEMALLOCCLASS
{
	R8BNOCTOR( CFixedBuffer );

public:
	CFixedBuffer()
		: Data( NULL )
	{
	}

	/**
	 * Constructor allocates memory so that the specified number of elements
	 * of type T can be stored in *this buffer object.
	 *
	 * @param Capacity Storage for this number of elements to allocate.
	 */

	CFixedBuffer( const int Capacity )
	{
		R8BASSERT( Capacity > 0 || Capacity == 0 );

		Data = (T*) allocmem( Capacity * sizeof( T ));

		R8BASSERT( Data != NULL || Capacity == 0 );
	}

	~CFixedBuffer()
	{
		freemem( Data );
	}

	/**
	 * Function allocates memory so that the specified number of elements of
	 * type T can be stored in *this buffer object.
	 *
	 * @param Capacity Storage for this number of elements to allocate.
	 */

	void alloc( const int Capacity )
	{
		R8BASSERT( Capacity > 0 || Capacity == 0 );

		freemem( Data );
		Data = (T*) allocmem( Capacity * sizeof( T ));

		R8BASSERT( Data != NULL || Capacity == 0 );
	}

	/**
	 * Function deallocates a previously allocated buffer.
	 */

	void free()
	{
		freemem( Data );
		Data = NULL;
	}

	/**
	 * @return Pointer to the first element of the allocated buffer, NULL if
	 * not allocated.
	 */

	T* getPtr() const
	{
		return( Data );
	}

	/**
	 * @return Pointer to the first element of the allocated buffer, NULL if
	 * not allocated.
	 */

	operator T* () const
	{
		return( Data );
	}

private:
	T* Data; ///< Element buffer pointer.
};

/**
 * \brief Spin lock-based thread synchonization class.
 *
 * Object for spin-lock multi-threaded synchronization. This object is
 * non-reentrant, it is invalid to enter it in the same thread twice: but the
 * effect will be evident: a dead-lock.
 *
 * This object is non-relocatable, it should not be copied nor moved after its
 * construction.
 *
 * Windows spin-lock implementation uses critical section since it's usually
 * only 2 times less efficient than a "manually coded" spin-lock while works
 * stably on both single and multi-processor systems.
 */

struct CSyncSpinLock
{
	R8BNOCTOR( CSyncSpinLock );

public:
	CSyncSpinLock()
	{
	#if defined( R8B_WIN )
		InitializeCriticalSectionAndSpinCount( &CritSec, 4000 );
	#elif defined( R8B_MAC )
		SpinLock = OS_SPINLOCK_INIT;
	#elif defined( R8B_LNX )
		pthread_spin_init( &SpinLock, PTHREAD_PROCESS_PRIVATE );
	#endif
	}

	~CSyncSpinLock()
	{
	#if defined( R8B_WIN )
		DeleteCriticalSection( &CritSec );
	#elif defined( R8B_LNX )
		pthread_spin_destroy( &SpinLock );
	#endif
	}

	/**
	 * Function acquires *this spin lock. Function does not return until the
	 * lock is obtained.
	 */

	void enter()
	{
	#if defined( R8B_WIN )
		EnterCriticalSection( &CritSec );
	#elif defined( R8B_MAC )
		OSSpinLockLock( &SpinLock );
	#elif defined( R8B_LNX )
		pthread_spin_lock( &SpinLock );
	#endif
	}

	/**
	 * Function releases *this spin lock.
	 */

	void leave()
	{
	#if defined( R8B_WIN )
		LeaveCriticalSection( &CritSec );
	#elif defined( R8B_MAC )
		OSSpinLockUnlock( &SpinLock );
	#elif defined( R8B_LNX )
		pthread_spin_unlock( &SpinLock );
	#endif
	}

private:
	#if defined( R8B_WIN )
		CRITICAL_SECTION CritSec; ///< Standard Windows critical section
			/// structure.
	#elif defined( R8B_MAC )
		volatile OSSpinLock SpinLock; ///< Spin-lock variable.
	#elif defined( R8B_LNX )
		pthread_spinlock_t SpinLock; ///< Spin-lock variable.
	#endif
};

/**
 * \brief A "keeper" class for spin lock-based synchronization.
 *
 * Sync keeper object for CSyncSpinLock objects.
 */

struct CSyncKeeperSpin
{
	R8BNOCTOR( CSyncKeeperSpin );

public:
	CSyncKeeperSpin()
		: SyncObj( NULL )
	{
	}

	/**
	 * @param aSyncObj Pointer to the sync object which should be used for
	 * sync'ing. Can be equal to NULL.
	 */

	CSyncKeeperSpin( CSyncSpinLock* const aSyncObj )
		: SyncObj( aSyncObj )
	{
		if( SyncObj != NULL )
		{
			SyncObj -> enter();
		}
	}

	/**
	 * @param aSyncObj Reference to the sync object which should be used for
	 * sync'ing.
	 */

	CSyncKeeperSpin( CSyncSpinLock& aSyncObj )
		: SyncObj( &aSyncObj )
	{
		SyncObj -> enter();
	}

	~CSyncKeeperSpin()
	{
		if( SyncObj != NULL )
		{
			SyncObj -> leave();
		}
	}

private:
	CSyncSpinLock* SyncObj; ///< Sync object which should be used (can be
		/// NULL).
};

/**
 * The synchronization macro. The R8BSYNCSPIN( obj ) macro should be put
 * before sections of the code that may be potentially changed asynchronously
 * by several threads simultaneously. The R8BSYNCSPIN( obj ) macro "enters"
 * the synchronization section that will block execution of other threads that
 * use same R8BSYNCSPIN( obj ) macro. The section begins with the
 * R8BSYNCSPIN( obj ) macro and finishes at the end of the current block. When
 * used, also works as the "memory barrier" instruction forcing caches of
 * processor cores to be synchronized.
 *
 * @param SyncObject An object of the CSyncSpinLock type that is used for
 * synchronization.
 */

#define R8BSYNCSPIN( SyncObject ) R8BSYNCSPIN_( SyncObject, __LINE__ )
#define R8BSYNCSPIN_( SyncObject, id ) R8BSYNCSPIN__( SyncObject, id )
#define R8BSYNCSPIN__( SyncObject, id ) \
	CSyncKeeperSpin SyncKeeper##id( SyncObject )

/**
 * \brief Sine signal generator class.
 *
 * Class implements sine signal generator with optional biasing.
 */

struct CSineGen
{
	/**
	 * Function initializes *this sine signal generator.
	 *
	 * @param si Sine function increment, in radians.
	 * @param g Overall gain. If biasing is not planned to be used, this value
	 * should be twice as high.
	 * @param ph Starting phase, in radians. Add 0.5 * M_PI for cosine
	 * function.
	 */

	void init( const double si, const double g,
		const double ph = -0.5 * M_PI )
	{
		bias = 0.5 * g;
		svalue1 = sin( ph ) * bias;
		svalue2 = sin( ph - si ) * bias;
		sincr = 2.0 * cos( si );
	}

	/**
	 * @return Next value of the sine function, with biasing.
	 */

	double genBias()
	{
		const double r = svalue1 + bias;

		const double tmp = svalue1;
		svalue1 = sincr * svalue1 - svalue2;
		svalue2 = tmp;

		return( r );
	}

	/**
	 * @return Next value of the sine function, without biasing.
	 */

	double gen()
	{
		const double r = svalue1;

		const double tmp = svalue1;
		svalue1 = sincr * svalue1 - svalue2;
		svalue2 = tmp;

		return( r );
	}

private:
	double bias; ///< Value bias.
	double svalue1; ///< Current sine value.
	double svalue2; ///< Previous sine value.
	double sincr; ///< Sine value increment.
};

/**
 * Function normalizes FIR filter so that its frequency response at DC is
 * equal to DCGain.
 *
 * @param[in,out] p Filter coefficients.
 * @param l Filter length.
 * @param DCGain Filter's gain at DC (linear, non-decibel value).
 * @param pstep "p" array step.
 */

template< class T >
inline void normalizeFIRFilter( T* const p, const int l, const double DCGain,
	const int pstep = 1 )
{
	double s = 0.0;
	T* pp = p;
	int i = l;

	while( i > 0 )
	{
		s += *pp;
		pp += pstep;
		i--;
	}

	s = DCGain / s;
	pp = p;
	i = l;

	while( i > 0 )
	{
		*pp *= s;
		pp += pstep;
		i--;
	}
}

/**
 * Function calculates frequency response of the specified FIR filter at the
 * specified circular frequency.
 *
 * @param flt FIR filter's coefficients.
 * @param fltlen Number of coefficients (taps) in the filter.
 * @param th Circular frequency [0; pi].
 * @param[out] re0 Resulting real part of the complex frequency response.
 * @param[out] im0 Resulting imaginary part of the complex frequency response.
 */

template< class T >
inline void calcFIRFilterResponse( const T* flt, int fltlen, const double th,
	double& re0, double& im0 )
{
	double re = 0.0;
	double im = 0.0;

	double svalue1 = 0.0;
	double svalue2 = sin( -th );
	double sincr = 2.0 * cos( th );
	double cvalue1 = 1.0;
	double cvalue2 = sin( M_PI * 0.5 - th );

	while( fltlen > 0 )
	{
		re += svalue1 * flt[ 0 ];
		im += cvalue1 * flt[ 0 ];
		flt++;
		fltlen--;

		double tmp = svalue1;
		svalue1 = sincr * svalue1 - svalue2;
		svalue2 = tmp;

		tmp = cvalue1;
		cvalue1 = sincr * cvalue1 - cvalue2;
		cvalue2 = tmp;
	}

	re0 = re;
	im0 = im;
}

/**
 * @param v Input value.
 * @return Calculated bit occupancy of the specified input value. Bit
 * occupancy means how many significant lower bits are necessary to store a
 * specified value. Function treats the input value as unsigned.
 */

inline int getBitOccupancy( const int v )
{
	static const int OccupancyTable[] =
	{
		1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4,
		5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
		6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
		6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
		8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
		8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
		8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
		8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
		8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
		8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
		8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
		8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8
	};

	int t, tt;
	tt = v >> 16;

	if( tt != 0 )
	{
		return(( t = v >> 24 ) ? 24 + OccupancyTable[ t ] :
			16 + OccupancyTable[ tt & 0xFF ]);
	}
	else
	{
		return(( t = v >> 8 ) ? 8 + OccupancyTable[ t ] :
			OccupancyTable[ v ]);
	}
}

#if !defined( min )

/**
 * @param v1 Value 1.
 * @param v2 Value 2.
 * @return The minimum of 2 values.
 */

template< class T >
inline T min( const T& v1, const T& v2 )
{
	return( v1 < v2 ? v1 : v2 );
}

#endif // min

#if !defined( max )

/**
 * @param v1 Value 1.
 * @param v2 Value 2.
 * @return The maximum of 2 values.
 */

template< class T >
inline T max( const T& v1, const T& v2 )
{
	return( v1 > v2 ? v1 : v2 );
}

#endif // max

/**
 * @param x Value to square.
 * @return Squared value of the argument.
 */

template< class T >
inline T sqr( const T x )
{
	return( x * x );
}

/**
 * @param v Input value.
 * @param p Power factor.
 * @return Returns pow() function's value with input value's sign check.
 */

template< class T >
inline T pows( const T v, const T p )
{
	return( v < 0.0 ? -pow( -v, p ) : pow( v, p ));
}

/**
 * @param v Input value.
 * @return Calculated single-argument Gaussian function of the input value.
 */

template< class T >
inline T gauss( const T v )
{
	return( exp( -( v * v )));
}

/**
 * @param v Input value.
 * @return Calculated inverse hyperbolic sine of the input value.
 */

template< class T >
inline T asinh( const T v )
{
	return( log( v + sqrt( v * v + 1.0 )));
}

} // namespace r8b

#endif // R8BBASE_INCLUDED
