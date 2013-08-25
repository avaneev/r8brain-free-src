//$ nocpp

/**
 * \file CDSPRealFFT.h
 *
 * \brief Real-valued FFT transform class.
 *
 * This file includes FFT object implementation. All created FFT objects are
 * kept in a global list after use for future reusal. Such approach minimizes
 * time necessary to initialize the FFT object of the required size.
 *
 * r8brain-free-src Copyright (c) 2013 Aleksey Vaneev
 * See the "License.txt" file for license.
 */

#ifndef R8B_CDSPREALFFT_INCLUDED
#define R8B_CDSPREALFFT_INCLUDED

#include "r8bbase.h"

#if !R8B_IPP
	#include "fft4g.h"
#endif // !R8B_IPP

namespace r8b {

/**
 * \brief Real-valued FFT transform class.
 *
 * Class implements a wrapper for real-valued discrete fast Fourier transform
 * functions. The object of this class can only be obtained via the
 * CDSPRealFFTKeeper class.
 *
 * Uses functions from the FFT package by: Copyright(C) 1996-2001 Takuya OOURA
 * http://www.kurims.kyoto-u.ac.jp/~ooura/fft.html
 *
 * Also uses Intel IPP library functions if available (the R8B_IPP=1 macro was
 * defined). Note that IPP library's FFT functions are 2-3 times more
 * efficient on the modern Intel Core i7-3770K processor than Ooura's
 * functions. It may be worthwhile investing in IPP. Note, that FFT functions
 * take less than 20% of the overall sample rate conversion time. However,
 * when the "power of 2" resampling is used the performance of FFT functions
 * becomes "everything".
 */

class CDSPRealFFT : public R8B_BASECLASS
{
	R8BNOCTOR( CDSPRealFFT );

	friend class CDSPRealFFTKeeper;

public:
	/**
	 * @return A multiplication constant that should be used after inverse
	 * transform to obtain a correct value scale.
	 */

	double getInvMulConst() const
	{
		return( InvMulConst );
	}

	/**
	 * @return The size (the number of real values in a transform) of *this
	 * FFT object, expressed as Nth power of 2.
	 */

	int getSizeBits() const
	{
		return( SizeBits );
	}

	/**
	 * @return The size (the number of real values in a transform) of *this
	 * FFT object.
	 */

	int getSize() const
	{
		return( Size );
	}

	/**
	 * Function performs in-place forward FFT.
	 *
	 * @param[in,out] p Pointer to data block to transform, size should be
	 * equal to *this object's Size.
	 */

	void forward( double* const p ) const
	{
	#if R8B_IPP

		ippsFFTFwd_RToPerm_64f( p, p, SPtr, WorkBuffer );

	#else // R8B_IPP

		ooura_fft :: rdft( Size, 1, p, ip.getPtr(), w.getPtr() );

	#endif // R8B_IPP
	}

	/**
	 * Function performs in-place inverse FFT.
	 *
	 * @param[in,out] p Pointer to data block to transform, size should be
	 * equal to *this object's Size.
	 */

	void inverse( double* const p ) const
	{
	#if R8B_IPP

		ippsFFTInv_PermToR_64f( p, p, SPtr, WorkBuffer );

	#else // R8B_IPP

		ooura_fft :: rdft( Size, -1, p, ip.getPtr(), w.getPtr() );

	#endif // R8B_IPP
	}

	/**
	 * Function multiplies two complex-valued data blocks and places result in
	 * a new data block. Length of all data blocks should be equal to *this
	 * object's block size. Input blocks should have been produced with the
	 * forward() function of *this object.
	 *
	 * @param ip1 Input data block 1.
	 * @param ip2 Input data block 2.
	 * @param[out] op Output data block, should not be equal to ip1 nor ip2.
	 */

	void multiplyBlocks( const double* const ip1, const double* const ip2,
		double* const op ) const
	{
		#if R8B_IPP

			ippsMulPerm_64f( (Ipp64f*) ip1, (Ipp64f*) ip2, (Ipp64f*) op,
				Size );

		#else // R8B_IPP

			op[ 0 ] = ip1[ 0 ] * ip2[ 0 ];
			op[ 1 ] = ip1[ 1 ] * ip2[ 1 ];

			int i = 2;

			while( i < Size )
			{
				op[ i ] = ip1[ i ] * ip2[ i ] - ip1[ i + 1 ] * ip2[ i + 1 ];
				op[ i + 1 ] =
					ip1[ i ] * ip2[ i + 1 ] + ip1[ i + 1 ] * ip2[ i ];

				i += 2;
			}

		#endif // R8B_IPP
	}

	/**
	 * Function is similar to the multiplyBlocks() function, but instead of
	 * replacing data in the output buffer, the data is summed with the output
	 * buffer.
	 *
	 * @param ip1 Input data block 1.
	 * @param ip2 Input data block 2.
	 * @param[out] op Output data block, should not be equal to ip1 nor ip2.
	 */

	void multiplyBlocksAdd( const double* const ip1, const double* const ip2,
		double* const op ) const
	{
		op[ 0 ] += ip1[ 0 ] * ip2[ 0 ];
		op[ 1 ] += ip1[ 1 ] * ip2[ 1 ];

		#if R8B_IPP

			ippsAddProduct_64fc( (const Ipp64fc*) ( ip1 + 2 ),
				(const Ipp64fc*) ( ip2 + 2 ), (Ipp64fc*) ( op + 2 ),
				( Size >> 1 ) - 1 );

		#else // R8B_IPP

			int i = 2;

			while( i < Size )
			{
				op[ i ] += ip1[ i ] * ip2[ i ] - ip1[ i + 1 ] * ip2[ i + 1 ];
				op[ i + 1 ] +=
					ip1[ i ] * ip2[ i + 1 ] + ip1[ i + 1 ] * ip2[ i ];

				i += 2;
			}

		#endif // R8B_IPP
	}

	/**
	 * Function multiplies two complex-valued data blocks in-place. Length of
	 * all data blocks should be equal to *this object's block size. Input
	 * blocks should have been produced with the forward() function of *this
	 * object.
	 *
	 * @param ip Input data block 1.
	 * @param[in,out] op Output/input data block 2.
	 */

	void multiplyBlocks( const double* const ip, double* const op ) const
	{
	#if R8B_IPP

		ippsMulPerm_64f( (Ipp64f*) op, (Ipp64f*) ip, (Ipp64f*) op, Size );

	#else // R8B_IPP

		op[ 0 ] *= ip[ 0 ];
		op[ 1 ] *= ip[ 1 ];

		int i = 2;

		while( i < Size )
		{
			const double t = op[ i ] * ip[ i ] - op[ i + 1 ] * ip[ i + 1 ];
			op[ i + 1 ] = op[ i ] * ip[ i + 1 ] + op[ i + 1 ] * ip[ i ];
			op[ i ] = t;
			i += 2;
		}

	#endif // R8B_IPP
	}

	/**
	 * Function performs in-place spectrum squaring. May cause aliasing
	 * if the filter was not zero-padded before the forward() function call.
	 *
	 * @param[in,out] p Pointer to data block to square, size should be equal
	 * to *this object's Size. This data block should contain complex spectrum
	 * data, previously obtained via the forward() function.
	 */

	void sqr( double* const p ) const
	{
		p[ 0 ] *= p[ 0 ];
		p[ 1 ] *= p[ 1 ];

	#if R8B_IPP

		ippsSqr_64fc( (Ipp64fc*) ( p + 2 ), (Ipp64fc*) ( p + 2 ),
			( Size >> 1 ) - 1 );

	#else // R8B_IPP

		int i = 2;

		while( i < Size )
		{
			const double r = p[ i ] * p[ i ] - p[ i + 1 ] * p[ i + 1 ];
			p[ i + 1 ] = p[ i ] * ( p[ i + 1 ] + p[ i + 1 ]);
			p[ i ] = r;
			i += 2;
		}

	#endif // R8B_IPP
	}

private:
	int SizeBits; ///< Size of FFT block (expressed as Nth power of 2).
		///<
	int Size; ///< Size of FFT block (number of real values).
		///<
	double InvMulConst; ///< Inverse FFT multiply constant.
		///<
	CDSPRealFFT* Next; ///< Next object in a singly-linked list.
		///<

	#if R8B_IPP
		IppsFFTSpec_R_64f* SPtr; ///< Pointer to initialized data buffer
			///< to be passed to IPP's FFT functions.
			///<
		CFixedBuffer< unsigned char > SpecBuffer; ///< Working buffer.
			///<
		CFixedBuffer< unsigned char > WorkBuffer; ///< Working buffer.
			///<
	#else // R8B_IPP
		CFixedBuffer< int > ip; ///< Working buffer (ints).
			///<
		CFixedBuffer< double > w; ///< Working buffer (doubles).
			///<
	#endif // R8B_IPP

	/**
	 * A simple class that keeps the pointer to the object and deletes it
	 * automatically.
	 */

	class CObjKeeper
	{
		R8BNOCTOR( CObjKeeper );

	public:
		CObjKeeper()
			: Object( NULL )
		{
		}

		~CObjKeeper()
		{
			delete Object;
		}

		CObjKeeper& operator = ( CDSPRealFFT* const aObject )
		{
			Object = aObject;
			return( *this );
		}

		operator CDSPRealFFT* () const
		{
			return( Object );
		}

	private:
		CDSPRealFFT* Object; ///< FFT object being kept.
			///<
	};

	CDSPRealFFT()
	{
	}

	/**
	 * Constructor initializes FFT object.
	 *
	 * @param SizeBits Size of FFT block (Nth power of 2), specifies the
	 * number of real values in a block. Values from 1 to 30 inclusive are
	 * supported.
	 */

	CDSPRealFFT( const int aSizeBits )
		: SizeBits( aSizeBits )
		, Size( 1 << aSizeBits )
	#if R8B_IPP
		, InvMulConst( 1.0 / (double) Size )
	#else // R8B_IPP
		, InvMulConst( 2.0 / (double) Size )
	#endif // R8B_IPP
	{
	#if R8B_IPP

		int SpecSize;
		int SpecBufferSize;
		int BufferSize;

		ippsFFTGetSize_R_64f( SizeBits, IPP_FFT_NODIV_BY_ANY,
			ippAlgHintFast, &SpecSize, &SpecBufferSize, &BufferSize );

		CFixedBuffer< unsigned char > InitBuffer( SpecBufferSize );
		SpecBuffer.alloc( SpecSize );
		WorkBuffer.alloc( BufferSize );

		ippsFFTInit_R_64f( &SPtr, SizeBits, IPP_FFT_NODIV_BY_ANY,
			ippAlgHintFast, SpecBuffer, InitBuffer );

	#else // R8B_IPP

		ip.alloc( (int) ceil( 2.0 + sqrt( (double) ( Size >> 1 ))));
		ip[ 0 ] = 0;
		w.alloc( Size >> 1 );

	#endif // R8B_IPP
	}

	~CDSPRealFFT()
	{
		delete Next;
	}
};

/**
 * \brief A "keeper" class for real-valued FFT transform objects.
 *
 * Class implements "keeper" functionality for handling CDSPRealFFT objects.
 * The allocated FFT objects are placed on the global static list of objects
 * for future reuse instead of deallocation.
 */

class CDSPRealFFTKeeper : public R8B_BASECLASS
{
	R8BNOCTOR( CDSPRealFFTKeeper );

public:
	CDSPRealFFTKeeper()
		: Object( NULL )
	{
	}

	/**
	 * Function acquires FFT object with the specified block size.
	 *
	 * @param SizeBits Size of FFT block (Nth power of 2), in the range
	 * [1; 30] inclusive, specifies the number of real values in a FFT block.
	 */

	CDSPRealFFTKeeper( const int SizeBits )
	{
		Object = acquire( SizeBits );
	}

	~CDSPRealFFTKeeper()
	{
		if( Object != NULL )
		{
			release( Object );
		}
	}

	/**
	 * @return Pointer to the acquired FFT object.
	 */

	const CDSPRealFFT* operator -> () const
	{
		R8BASSERT( Object != NULL );

		return( Object );
	}

	/**
	 * Function acquires FFT object with the specified block size. This
	 * function can be called any number of times.
	 *
	 * @param SizeBits Size of FFT block (Nth power of 2), in the range
	 * [1; 30] inclusive, specifies the number of real values in a FFT block.
	 */

	void init( const int SizeBits )
	{
		if( Object != NULL )
		{
			if( Object -> SizeBits == SizeBits )
			{
				return;
			}

			release( Object );
		}

		Object = acquire( SizeBits );
	}

	/**
	 * Function releases a previously acquired FFT object.
	 */

	void reset()
	{
		if( Object != NULL )
		{
			release( Object );
			Object = NULL;
		}
	}

private:
	CDSPRealFFT* Object; ///< FFT object.
		///<

	static CSyncObject StateSync; ///< FFTObjects synchronizer.
		///<
	static CDSPRealFFT :: CObjKeeper FFTObjects[]; ///< Pool of FFT objects of
		///< various sizes.
		///<

	/**
	 * Function acquires FFT object from the global pool.
	 *
	 * @param SizeBits FFT block size (expressed as Nth power of 2).
	 */

	CDSPRealFFT* acquire( const int SizeBits )
	{
		R8BASSERT( SizeBits > 0 && SizeBits <= 30 );

		R8BSYNC( StateSync );

		if( FFTObjects[ SizeBits ] == NULL )
		{
			return( new CDSPRealFFT( SizeBits ));
		}

		CDSPRealFFT* ffto = FFTObjects[ SizeBits ];
		FFTObjects[ SizeBits ] = ffto -> Next;

		return( ffto );
	}

	/**
	 * Function releases a previously acquired FFT object.
	 *
	 * @param ffto FFT object to release.
	 */

	void release( CDSPRealFFT* const ffto )
	{
		R8BSYNC( StateSync );

		ffto -> Next = FFTObjects[ ffto -> SizeBits ];
		FFTObjects[ ffto -> SizeBits ] = ffto;
	}
};

/**
 * Function calculates the minimum-phase transform of the filter kernel, using
 * a discrete Hilbert transform in cepstrum domain.
 *
 * For more details, see part III.B of
 * http://www.hpl.hp.com/personal/Niranjan_Damera-Venkata/files/ComplexMinPhase.pdf
 *
 * @param[in,out] Kernel Filter kernel buffer.
 * @param KernelLen Filter kernel's length, in samples.
 * @param LenMult Kernel length multiplier. Used as a coefficient of the
 * "oversampling" in the frequency domain. Such oversampling is needed to
 * improve the precision of the minimum-phase transform. If the filter's
 * attenuation is high, this multiplier should be increased or otherwise the
 * required attenuation will not be reached due to "smoothing" effect of this
 * transform.
 * @param DoFinalMul "True" if the final multiplication after transform should
 * be performed or not. Such multiplication returns the gain of the signal to
 * its original value. This parameter can be set to "false" if normalization
 * of the resulting filter kernel is planned to be used.
 * @param[out] DCGroupDelay If not NULL, this variable receives group delay
 * at DC offset, in samples (can be a non-integer value).
 */

inline void calcMinPhaseTransform( double* const Kernel, const int KernelLen,
	const int LenMult = 2, const bool DoFinalMul = true,
	double* const DCGroupDelay = NULL )
{
	R8BASSERT( KernelLen > 0 );
	R8BASSERT( LenMult >= 2 );

	const int SizeBits = getBitOccupancy(( KernelLen * LenMult ) - 1 );
	const int Size = 1 << SizeBits;
	const int Size2 = Size >> 1;
	int i;

	CFixedBuffer< double > ip( Size );
	CFixedBuffer< double > ip2( Size2 + 1 );

	memcpy( &ip[ 0 ], Kernel, KernelLen * sizeof( double ));
	memset( &ip[ KernelLen ], 0, ( Size - KernelLen ) * sizeof( double ));

	CDSPRealFFTKeeper ffto( SizeBits );
	ffto -> forward( ip );

	// Create the "log |c|" spectrum while saving the original power spectrum
	// in the "ip2" buffer.

	ip2[ 0 ] = ip[ 0 ];
	ip[ 0 ] = log( fabs( ip[ 0 ]));
	ip2[ Size2 ] = ip[ 1 ];
	ip[ 1 ] = log( fabs( ip[ 1 ]));

	for( i = 1; i < Size2; i++ )
	{
		ip2[ i ] = sqrt( ip[ i * 2 ] * ip[ i * 2 ] +
			ip[ i * 2 + 1 ] * ip[ i * 2 + 1 ]);

		ip[ i * 2 ] = log( ip2[ i ]);
		ip[ i * 2 + 1 ] = 0.0;
	}

	// Convert to cepstrum and apply discrete Hilbert transform.

	ffto -> inverse( ip );

	ip[ 0 ] = 0.0;

	for( i = 1; i < Size2; i++ )
	{
		ip[ i ] *= ffto -> getInvMulConst();
	}

	ip[ Size2 ] = 0.0;

	for( i = Size2 + 1; i < Size; i++ )
	{
		ip[ i ] *= -ffto -> getInvMulConst();
	}

	// Convert Hilbert-transformed cepstrum back to the "log |c|" spectrum and
	// perform its exponentiation, multiplied by the power spectrum previously
	// saved in the "ip2" buffer.

	ffto -> forward( ip );

	ip[ 0 ] = ip2[ 0 ];
	ip[ 1 ] = ip2[ Size2 ];

	for( i = 1; i < Size2; i++ )
	{
		const double p = ip2[ i ];
		ip[ i * 2 + 0 ] = cos( ip[ i * 2 + 1 ]) * p;
		ip[ i * 2 + 1 ] = sin( ip[ i * 2 + 1 ]) * p;
	}

	ffto -> inverse( ip );

	if( DoFinalMul )
	{
		for( i = 0; i < KernelLen; i++ )
		{
			Kernel[ i ] = ip[ i ] * ffto -> getInvMulConst();
		}
	}
	else
	{
		memcpy( &Kernel[ 0 ], &ip[ 0 ], KernelLen * sizeof( double ));
	}

	if( DCGroupDelay != NULL )
	{
		double tmp;

		calcFIRFilterResponseAndGroupDelay( Kernel, KernelLen, 0.0,
			tmp, tmp, *DCGroupDelay );
	}
/*
	const int thcount = 10000;
	const double g = 10.0 / log( 10.0 );

	for( i = 0; i < thcount; i++ )
	{
		const double th = M_PI * i / ( thcount - 1 );
		double re;
		double im;
		double gd;

		calcFIRFilterResponseAndGroupDelay( Kernel, KernelLen, th,
			re, im, gd );

		printf( "%.3f\t%.10f\n", g * log( re * re + im * im ), gd );
//		printf( "%.3f\n", g * log( re * re + im * im ));
	}*/
}

} // namespace r8b

#endif // VOX_CDSPREALFFT_INCLUDED
