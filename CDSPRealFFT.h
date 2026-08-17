//$ nobt
//$ nocpp

/**
 * @file CDSPRealFFT.h
 *
 * @brief Real-valued FFT transform class.
 *
 * This file includes FFT object implementation. All created FFT objects are
 * kept in a global list after use, for a future reusal. Such approach
 * minimizes time necessary to initialize the FFT object of the required
 * length.
 *
 * r8brain-free-src Copyright (c) 2013-2026 Aleksey Vaneev
 *
 * See the "LICENSE" file for license.
 */

#ifndef R8B_CDSPREALFFT_INCLUDED
#define R8B_CDSPREALFFT_INCLUDED

#include "r8bbase.h"

#if R8B_PFFFT_DOUBLE
	#include "fft/pffft_double.h"
#elif R8B_PFFFT
	#include "fft/pffft.h"
#elif !R8B_IPP
	#include "fft/fft4g.h"
#endif // !R8B_IPP

namespace r8b {

#if R8B_PFFFT
	typedef float realfft_t; ///< Forward transform's native type.
#else // R8B_PFFFT
	typedef double realfft_t; ///< Forward transform's native type.
#endif // R8B_PFFFT

/**
 * @brief Real-valued FFT setup class.
 *
 * This class stores pre-calculated data for FFT functions of a specific
 * length. Only a single object of this class is required for all FFT
 * transform objects utilized concurrently.
 */

class CDSPRealFFTSetup : public R8B_BASECLASS
{
	R8BNOCTOR( CDSPRealFFTSetup )

	friend class CPtrKeeper< CDSPRealFFTSetup >;
	friend class CDSPRealFFT;
	friend class CDSPRealFFTKeeper;

public:
	/**
	 * @brief Return a multiplication constant that should be used after
	 * inverse transform to obtain a correct value scale.
	 */

	double getInvMulConst() const
	{
		return( InvMulConst );
	}

	/**
	 * @brief Returns the length (the number of real values in a transform) of
	 * *this* FFT setup object, expressed as Nth power of 2.
	 */

	int getLenBits() const
	{
		return( LenBits );
	}

	/**
	 * @brief Returns the length (the number of real values in a transform) of
	 * *this* FFT setup object.
	 */

	int getLen() const
	{
		return( Len );
	}

private:
	int LenBits; ///< Length of FFT block (expressed as Nth power of 2).
	int Len; ///< Length of FFT block (number of real values).
	double InvMulConst; ///< Inverse FFT multiply constant.

	#if R8B_IPP
		IppsFFTSpec_R_64f* SpecPtr; ///< Pointer to the initialized
			///< specification buffer to be passed to IPP's FFT functions.
		CFixedBuffer< unsigned char > SpecBuffer; ///< Working buffer.
		int WorkBufferSize; ///< Size of the working buffer.
	#elif R8B_PFFFT
		PFFFT_Setup* setup; ///< PFFFT setup object.
	#elif R8B_PFFFT_DOUBLE
		PFFFTD_Setup* setup; ///< PFFFTD setup object.
	#else // R8B_PFFFT_DOUBLE
		CFixedBuffer< int > wi; ///< Coefficients buffer (ints).
		CFixedBuffer< double > wd; ///< Coefficients buffer (doubles).
	#endif // R8B_PFFFT_DOUBLE

	CDSPRealFFTSetup()
	{
	}

	/**
	 * @brief Initializes *this* FFT setup object.
	 *
	 * @param aLenBits The length of FFT block (Nth power of 2), specifies the
	 * number of real values in a block. Values from 1 to 30 inclusive are
	 * supported.
	 */

	CDSPRealFFTSetup( const int aLenBits )
		: LenBits( aLenBits )
		, Len( 1 << aLenBits )
	#if R8B_OOURA
		, InvMulConst( 2.0 / Len )
	#else // R8B_OOURA
		, InvMulConst( 1.0 / Len )
	#endif // R8B_OOURA
	{
	#if R8B_IPP

		int SpecBufferSize;
		int InitBufferSize;

		ippsFFTGetSize_R_64f( LenBits, IPP_FFT_NODIV_BY_ANY,
			ippAlgHintFast, &SpecBufferSize, &InitBufferSize,
			&WorkBufferSize );

		CFixedBuffer< unsigned char > InitBuffer( InitBufferSize );
		SpecBuffer.alloc( SpecBufferSize );

		ippsFFTInit_R_64f( &SpecPtr, LenBits, IPP_FFT_NODIV_BY_ANY,
			ippAlgHintFast, SpecBuffer, InitBuffer );

	#elif R8B_PFFFT

		setup = pffft_new_setup( Len, PFFFT_REAL );

	#elif R8B_PFFFT_DOUBLE

		setup = pffftd_new_setup( Len, PFFFT_REAL );

	#else // R8B_PFFFT_DOUBLE

		wi.alloc( (int) ceil( 2.0 + sqrt( (double) ( Len >> 1 ))));
		wd.alloc( Len >> 1 );

		ooura_fft :: rdftinit( Len, wi, wd );

	#endif // R8B_PFFFT_DOUBLE
	}

	~CDSPRealFFTSetup()
	{
		#if R8B_PFFFT
			pffft_destroy_setup( setup );
		#elif R8B_PFFFT_DOUBLE
			pffftd_destroy_setup( setup );
		#endif // R8B_PFFFT_DOUBLE
	}
};

/**
 * @brief Real-valued FFT transform class.
 *
 * Class implements a wrapper for real-valued discrete fast Fourier transform
 * functions. The object of this class can only be obtained via the
 * CDSPRealFFTKeeper class.
 *
 * Uses functions from the FFT package by: Copyright(C) 1996-2001 Takuya OOURA
 * http://www.kurims.kyoto-u.ac.jp/~ooura/fft.html
 *
 * Also uses Intel IPP library functions if available (if the R8B_IPP=1 macro
 * was defined). Note that IPP library's FFT functions are 2-3 times more
 * efficient on the modern Intel Core i7-3770K processor than Ooura's
 * functions. It may be worthwhile investing in IPP. Note, that FFT functions
 * take less than 20% of the overall sample rate conversion time. However,
 * when the "power of 2" resampling is used the performance of FFT functions
 * becomes "everything".
 */

class CDSPRealFFT : public R8B_BASECLASS
{
	R8BNOCTOR( CDSPRealFFT )

	friend class CPtrKeeper< CDSPRealFFT >;
	friend class CDSPRealFFTKeeper;

public:
	/**
	 * @brief Return a multiplication constant that should be used after
	 * inverse transform to obtain a correct value scale.
	 */

	double getInvMulConst() const
	{
		return( s -> InvMulConst );
	}

	/**
	 * @brief Returns the length (the number of real values in a transform) of
	 * *this* FFT object, expressed as Nth power of 2.
	 */

	int getLenBits() const
	{
		return( s -> LenBits );
	}

	/**
	 * @brief Returns the length (the number of real values in a transform) of
	 * *this* FFT object.
	 */

	int getLen() const
	{
		return( Len );
	}

	/**
	 * @brief Performs in-place forward FFT.
	 *
	 * If FFT's native type differs from `double`, the result will be placed
	 * to `p` using values of the native type.
	 *
	 * @param[in,out] p Pointer to data block to transform, length should be
	 * equal to *this* object's getLen().
	 * @return The pointer `p` casted to transform's native type.
	 */

	realfft_t* forward( double* const p ) const
	{
	#if R8B_IPP

		ippsFFTFwd_RToPerm_64f( p, p, s -> SpecPtr, WorkBuffer );

		return( p );

	#elif R8B_PFFFT

		float* const tp = &work[ Len ];
		int i;

		for( i = 0; i < Len; i++ )
		{
			tp[ i ] = (float) p[ i ];
		}

		float* const op = constructptr< float >( (void*) p, (size_t) Len );

		pffft_transform_ordered( s -> setup, tp, op, work, PFFFT_FORWARD );

		return( op );

	#elif R8B_PFFFT_DOUBLE

		pffftd_transform_ordered( s -> setup, p, p, work, PFFFT_FORWARD );

		return( p );

	#else // R8B_PFFFT_DOUBLE

		ooura_fft :: rdft( Len, 1, p, s -> wi, s -> wd );

		return( p );

	#endif // R8B_PFFFT_DOUBLE
	}

	/**
	 * @brief Performs in-place inverse FFT.
	 *
	 * If FFT's native type differs from `double`, the result will be placed
	 * to `p` using `double` values.
	 *
	 * @param[in,out] p Pointer to data block to transform, length should be
	 * equal to *this* object's getLen().
	 */

	void inverse( realfft_t* const p ) const
	{
	#if R8B_IPP

		ippsFFTInv_PermToR_64f( p, p, s -> SpecPtr, WorkBuffer );

	#elif R8B_PFFFT

		float* const tp = &work[ Len ];

		pffft_transform_ordered( s -> setup, p, tp, work, PFFFT_BACKWARD );

		double* const op = constructptr< double >( (void*) p, (size_t) Len );
		int i;

		for( i = 0; i < Len; i++ )
		{
			op[ i ] = tp[ i ];
		}

	#elif R8B_PFFFT_DOUBLE

		pffftd_transform_ordered( s -> setup, p, p, work, PFFFT_BACKWARD );

	#else // R8B_PFFFT_DOUBLE

		ooura_fft :: rdft( Len, -1, p, s -> wi, s -> wd );

	#endif // R8B_PFFFT_DOUBLE
	}

	/**
	 * @brief Multiplies two complex-valued data blocks and places result in
	 * a new data block.
	 *
	 * Length of all data blocks should be equal to *this* object's block
	 * length. Input blocks should have been produced with the forward()
	 * function of *this* object.
	 *
	 * @param ip1 Input data block 1.
	 * @param ip2 Input data block 2.
	 * @param[out] op Output data block, should not be equal to aip1 nor
	 * aip2.
	 */

	void multiplyBlocks( const realfft_t* const ip1,
		const realfft_t* const ip2, realfft_t* const op ) const
	{
	#if R8B_IPP

		ippsMulPerm_64f( (Ipp64f*) ip1, (Ipp64f*) ip2, (Ipp64f*) op, Len );

	#else // R8B_IPP

		op[ 0 ] = ip1[ 0 ] * ip2[ 0 ];
		op[ 1 ] = ip1[ 1 ] * ip2[ 1 ];

		int i = 2;

		while( i < Len )
		{
			op[ i ] = ip1[ i ] * ip2[ i ] - ip1[ i + 1 ] * ip2[ i + 1 ];
			op[ i + 1 ] = ip1[ i ] * ip2[ i + 1 ] + ip1[ i + 1 ] * ip2[ i ];
			i += 2;
		}

	#endif // R8B_IPP
	}

	/**
	 * @brief Multiplies two complex-valued data blocks in-place.
	 *
	 * Length of both data blocks should be equal to *this* object's block
	 * length. Blocks should have been produced with the forward() function of
	 * *this* object.
	 *
	 * @param ip Input data block 1.
	 * @param[in,out] op Output/input data block 2.
	 */

	void multiplyBlocks( const realfft_t* const ip, realfft_t* const op ) const
	{
	#if R8B_IPP

		ippsMulPerm_64f( (Ipp64f*) op, (Ipp64f*) ip, (Ipp64f*) op, Len );

	#else // R8B_IPP

		op[ 0 ] *= ip[ 0 ];
		op[ 1 ] *= ip[ 1 ];

		realfft_t t;
		int i = 2;

		while( i < Len )
		{
			t = op[ i ] * ip[ i ] - op[ i + 1 ] * ip[ i + 1 ];
			op[ i + 1 ] = op[ i ] * ip[ i + 1 ] + op[ i + 1 ] * ip[ i ];
			op[ i ] = t;
			i += 2;
		}

	#endif // R8B_IPP
	}

	/**
	 * @brief Multiplies two complex-valued data blocks in-place, considering
	 * that the `aip` block contains "zero-phase" response.
	 *
	 * Length of both data blocks should be equal to *this* object's block
	 * length. Blocks should have been produced with the forward() function of
	 * *this* object.
	 *
	 * @param ip Input data block 1, "zero-phase" response. This block should
	 * be first transformed via the convertToZP() function.
	 * @param[in,out] op Output/input data block 2.
	 */

	void multiplyBlocksZP( const realfft_t* ip, realfft_t* op ) const
	{
	// SIMD implementations assume that pointers are address-aligned.

	#if !R8B_PFFFT && defined( R8B_SSE2 )

		int c8 = Len >> 3;

		while( c8 != 0 )
		{
			const __m128d iv1 = _mm_load_pd( ip );
			const __m128d iv2 = _mm_load_pd( ip + 2 );
			const __m128d ov1 = _mm_load_pd( op );
			const __m128d ov2 = _mm_load_pd( op + 2 );
			_mm_store_pd( op, _mm_mul_pd( iv1, ov1 ));
			_mm_store_pd( op + 2, _mm_mul_pd( iv2, ov2 ));

			const __m128d iv3 = _mm_load_pd( ip + 4 );
			const __m128d ov3 = _mm_load_pd( op + 4 );
			const __m128d iv4 = _mm_load_pd( ip + 6 );
			const __m128d ov4 = _mm_load_pd( op + 6 );
			_mm_store_pd( op + 4, _mm_mul_pd( iv3, ov3 ));
			_mm_store_pd( op + 6, _mm_mul_pd( iv4, ov4 ));

			ip += 8;
			op += 8;
			c8--;
		}

		int c = Len & 7;

		while( c != 0 )
		{
			*op *= *ip;
			ip++;
			op++;
			c--;
		}

	#elif !R8B_PFFFT && defined( R8B_NEON )

		int c8 = Len >> 3;

		while( c8 != 0 )
		{
			const float64x2_t iv1 = vld1q_f64( ip );
			const float64x2_t iv2 = vld1q_f64( ip + 2 );
			const float64x2_t ov1 = vld1q_f64( op );
			const float64x2_t ov2 = vld1q_f64( op + 2 );
			vst1q_f64( op, vmulq_f64( iv1, ov1 ));
			vst1q_f64( op + 2, vmulq_f64( iv2, ov2 ));

			const float64x2_t iv3 = vld1q_f64( ip + 4 );
			const float64x2_t iv4 = vld1q_f64( ip + 6 );
			const float64x2_t ov3 = vld1q_f64( op + 4 );
			const float64x2_t ov4 = vld1q_f64( op + 6 );
			vst1q_f64( op + 4, vmulq_f64( iv3, ov3 ));
			vst1q_f64( op + 6, vmulq_f64( iv4, ov4 ));

			ip += 8;
			op += 8;
			c8--;
		}

		int c = Len & 7;

		while( c != 0 )
		{
			*op *= *ip;
			ip++;
			op++;
			c--;
		}

	#else // SIMD

		int i;

		for( i = 0; i < Len; i++ )
		{
			op[ i ] *= ip[ i ];
		}

	#endif // SIMD
	}

	/**
	 * @brief Converts the specified forward-transformed block into
	 * "zero-phase" form, suitable for use with the multiplyBlocksZP()
	 * function.
	 *
	 * @param[in,out] p Block to transform.
	 */

	void convertToZP( realfft_t* const p ) const
	{
		int i = 2;

		while( i < Len )
		{
			p[ i + 1 ] = p[ i ];
			i += 2;
		}
	}

private:
	CDSPRealFFTSetup* s; ///< Pointer to static FFT transform setup object.
	int Len; ///< Length of FFT block (number of real values).
	CDSPRealFFT* Next; ///< Next object in a singly-linked list.

	#if R8B_IPP
		CFixedBuffer< unsigned char > WorkBuffer; ///< Working buffer.
	#elif R8B_PFFFT
		CFixedBuffer< float > work; ///< Working buffer, includes temporary
			///< buffer for double-to-float conversions.
	#elif R8B_PFFFT_DOUBLE
		CFixedBuffer< double > work; ///< Working buffer.
	#endif // R8B_PFFFT_DOUBLE

	/**
	 * Constructor initializes FFT object.
	 *
	 * @param s0 Pointer to FFT setup object.
	 */

	CDSPRealFFT( CDSPRealFFTSetup* const s0 )
		: s( s0 )
		, Len( s -> Len )
		, Next( R8B_NULL )
	{
	#if R8B_IPP

		WorkBuffer.alloc( s -> WorkBufferSize );

	#elif R8B_PFFFT

		work.alloc( Len * 2 );

	#elif R8B_PFFFT_DOUBLE

		work.alloc( Len );

	#endif // R8B_PFFFT_DOUBLE
	}

	~CDSPRealFFT()
	{
		while( Next != R8B_NULL )
		{
			CDSPRealFFT* const nn = Next -> Next;
			Next -> Next = R8B_NULL;
			delete Next;
			Next = nn;
		}
	}
};

/**
 * @brief A "keeper" class for real-valued FFT transform objects.
 *
 * Class implements "keeper" functionality for handling CDSPRealFFT objects.
 * The allocated FFT objects are placed on the global static list of objects
 * for future reuse instead of deallocation.
 */

class CDSPRealFFTKeeper : public R8B_BASECLASS
{
	R8BNOCTOR( CDSPRealFFTKeeper )

public:
	CDSPRealFFTKeeper()
		: Object( R8B_NULL )
	{
	}

	/**
	 * @brief Acquires FFT object with the specified block length.
	 *
	 * @param LenBits The length of FFT block (Nth power of 2), in the range
	 * [1; 30] inclusive, specifies the number of real values in a FFT block.
	 */

	CDSPRealFFTKeeper( const int LenBits )
	{
		Object = acquire( LenBits );
	}

	~CDSPRealFFTKeeper()
	{
		if( Object != R8B_NULL )
		{
			release( Object );
		}
	}

	/**
	 * @brief Returns pointer to the acquired FFT object.
	 */

	const CDSPRealFFT* operator -> () const
	{
		R8BASSERT( Object != R8B_NULL );

		return( Object );
	}

	/**
	 * @brief Acquires FFT object with the specified block length. This
	 * function can be called any number of times.
	 *
	 * @param LenBits The length of FFT block (Nth power of 2), in the range
	 * [1; 30] inclusive, specifies the number of real values in a FFT block.
	 */

	void init( const int LenBits )
	{
		if( Object != R8B_NULL )
		{
			if( Object -> getLenBits() == LenBits )
			{
				return;
			}

			release( Object );
		}

		Object = acquire( LenBits );
	}

	/**
	 * @brief Releases a previously acquired FFT object.
	 */

	void reset()
	{
		if( Object != R8B_NULL )
		{
			release( Object );
			Object = R8B_NULL;
		}
	}

private:
	CDSPRealFFT* Object; ///< FFT object.

	/**
	 * @brief Acquires FFT object from the global pool.
	 *
	 * @param LenBits FFT block length (expressed as Nth power of 2).
	 */

	CDSPRealFFT* acquire( const int LenBits )
	{
		R8BASSERT( LenBits > 0 && LenBits <= 30 );

		R8BSYNC( getStateSync() );

		CDSPRealFFTSetup* s = getFFTSetupObjects()[ LenBits ];

		if( s == R8B_NULL )
		{
			s = new CDSPRealFFTSetup( LenBits );
			getFFTSetupObjects()[ LenBits ] = s;
		}

		#if R8B_OOURA

			if( getFFTObjects()[ LenBits ] == R8B_NULL )
			{
				getFFTObjects()[ LenBits ] = new CDSPRealFFT( s );
			}

			return( getFFTObjects()[ LenBits ]);

		#else // R8B_OOURA

			if( getFFTObjects()[ LenBits ] == R8B_NULL )
			{
				return( new CDSPRealFFT( s ));
			}

			CDSPRealFFT* const ffto = getFFTObjects()[ LenBits ].unkeep();
			getFFTObjects()[ LenBits ] = ffto -> Next;
			ffto -> Next = R8B_NULL;

			return( ffto );

		#endif // R8B_OOURA
	}

	/**
	 * @brief Releases a previously acquired FFT object.
	 *
	 * @param ffto FFT object to release.
	 */

	void release( CDSPRealFFT* const ffto )
	{
	#if R8B_OOURA

		(void) ffto;

	#else // R8B_OOURA

		R8BSYNC( getStateSync() );

		ffto -> Next = getFFTObjects()[ ffto -> getLenBits() ].unkeep();
		getFFTObjects()[ ffto -> getLenBits() ] = ffto;

	#endif // R8B_OOURA
	}

	/**
	 * @brief Returns pointer to a pool of FFT setup objects.
	 */

	static CPtrKeeper< CDSPRealFFTSetup >* getFFTSetupObjects()
	{
		R8B_EXITDTOR static CPtrKeeper< CDSPRealFFTSetup > FFTSetups[ 31 ];

		return( FFTSetups );
	}

	/**
	 * @brief Returns pointer to a pool of free FFT objects.
	 */

	static CPtrKeeper< CDSPRealFFT >* getFFTObjects()
	{
		R8B_EXITDTOR static CPtrKeeper< CDSPRealFFT > FFTObjects[ 31 ];

		return( FFTObjects );
	}

	/**
	 * @brief Returns reference to FFT object pool sync object.
	 */

	static CSyncObject& getStateSync()
	{
		R8B_EXITDTOR static CSyncObject StateSync;

		return( StateSync );
	}
};

/**
 * @brief Calculates the minimum-phase transform of the filter kernel, using
 * a discrete Hilbert transform in cepstrum domain.
 *
 * For more details, see part III.B of
 * http://www.hpl.hp.com/personal/Niranjan_Damera-Venkata/files/ComplexMinPhase.pdf
 *
 * @param[in,out] Kernel Filter kernel buffer.
 * @param KernelLen Filter kernel's length, in samples.
 * @param LenMult Kernel length multiplier. Used as a coefficient of
 * oversampling in the frequency domain. Such oversampling is needed to
 * improve the precision of the minimum-phase transform. If the filter's
 * attenuation is high, this multiplier should be increased or otherwise the
 * required attenuation will not be reached due to "smoothing" effect of this
 * transform.
 * @param DoFinalMul `true` if the final multiplication after transform
 * should be performed. Such multiplication returns the gain of the signal to
 * its original value. This parameter can be set to `false` if normalization
 * of the resulting filter kernel is planned to be used.
 * @param[out] DCGroupDelay If not `nullptr`, this variable receives group
 * delay at DC offset, in samples (can be a non-integer value).
 */

inline void calcMinPhaseTransform( double* const Kernel, const int KernelLen,
	const int LenMult = 2, const bool DoFinalMul = true,
	double* const DCGroupDelay = R8B_NULL )
{
	R8BASSERT( KernelLen > 0 );
	R8BASSERT( LenMult >= 2 );

	const int LenBits = getBitOccupancy(( KernelLen * LenMult ) - 1 );
	const int Len = 1 << LenBits;
	const int Len2 = Len >> 1;
	int i;

	CFixedBuffer< double > ip( Len );
	CFixedBuffer< realfft_t > ip2( Len2 + 1 );

	memcpy( &ip[ 0 ], Kernel, (size_t) KernelLen * sizeof( ip[ 0 ]));
	memset( &ip[ KernelLen ], 0,
		(size_t) ( Len - KernelLen ) * sizeof( ip[ 0 ]));

	CDSPRealFFTKeeper ffto( LenBits );

	realfft_t* aip = ffto -> forward( ip );
	realfft_t* const aip2 = &ip2[ 0 ];

	const realfft_t nzbias = ( sizeof( realfft_t ) < sizeof( double ) ?
		1e-35 : 1e-300 );

	// Create the "log |c|" spectrum while saving the original power spectrum
	// in the "ip2" buffer.

	aip2[ 0 ] = aip[ 0 ];
	aip[ 0 ] = log( fabs( aip[ 0 ]) + nzbias );
	aip2[ Len2 ] = aip[ 1 ];
	aip[ 1 ] = log( fabs( aip[ 1 ]) + nzbias );

	for( i = 1; i < Len2; i++ )
	{
		aip2[ i ] = sqrt( aip[ i * 2 ] * aip[ i * 2 ] +
			aip[ i * 2 + 1 ] * aip[ i * 2 + 1 ]);

		aip[ i * 2 ] = log( aip2[ i ] + nzbias );
		aip[ i * 2 + 1 ] = 0.0;
	}

	// Convert to cepstrum and apply discrete Hilbert transform.

	ffto -> inverse( aip );

	const double m1 = ffto -> getInvMulConst();
	const double m2 = -m1;

	ip[ 0 ] = 0.0;

	for( i = 1; i < Len2; i++ )
	{
		ip[ i ] *= m1;
	}

	ip[ Len2 ] = 0.0;

	for( i = Len2 + 1; i < Len; i++ )
	{
		ip[ i ] *= m2;
	}

	// Convert Hilbert-transformed cepstrum back to the "log |c|" spectrum and
	// perform its exponentiation, multiplied by the power spectrum previously
	// saved in the "ip2" buffer.

	aip = ffto -> forward( ip );

	aip[ 0 ] = aip2[ 0 ];
	aip[ 1 ] = aip2[ Len2 ];

	for( i = 1; i < Len2; i++ )
	{
		aip[ i * 2 + 0 ] = cos( aip[ i * 2 + 1 ]) * aip2[ i ];
		aip[ i * 2 + 1 ] = sin( aip[ i * 2 + 1 ]) * aip2[ i ];
	}

	ffto -> inverse( aip );

	if( DoFinalMul )
	{
		for( i = 0; i < KernelLen; i++ )
		{
			Kernel[ i ] = ip[ i ] * m1;
		}
	}
	else
	{
		memcpy( &Kernel[ 0 ], &ip[ 0 ],
			(size_t) KernelLen * sizeof( Kernel[ 0 ]));
	}

	if( DCGroupDelay != R8B_NULL )
	{
		*DCGroupDelay = calcFIRFilterGroupDelay( Kernel, KernelLen, 0.0 );
	}
}

} // namespace r8b

#endif // R8B_CDSPREALFFT_INCLUDED
