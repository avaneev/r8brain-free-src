//$ nocpp

/**
 * \file CDSPFracInterpolator.h
 *
 * \brief Fractional delay interpolator and filter bank classes.
 *
 * This file includes fractional delay interpolator class.
 *
 * r8brain-free-src Copyright (c) 2013 Aleksey Vaneev
 * See the "License.txt" file for license.
 */

#ifndef R8B_CDSPFRACINTERPOLATOR_INCLUDED
#define R8B_CDSPFRACINTERPOLATOR_INCLUDED

#include "CDSPSincFilterGen.h"

namespace r8b {

/**
 * \brief Sinc function-based fractional delay filter bank class.
 *
 * Class implements storage and initialization of a bank of sinc-based
 * fractional delay filters, expressed as 4-point Hermite spline coefficients
 * of 3rd order. The filters are windowed by one of the "Vaneev" windowing
 * functions. The FilterLen and FilterFracs parameters can be varied freely
 * without breaking the resampler.
 *
 * @param FilterLen Specifies the number of samples (taps) each fractional
 * delay filter should have. This must be an even value, the minimal value for
 * FilterLen is 4. To achieve a higher resampling precision, the oversampling
 * should be used in the first place instead of using a higher FilterLen
 * value. The lower this value is the lower the signal-to-noise performance of
 * the interpolator will be.
 * @param FilterFracs The number of fractional delay positions to sample. This
 * do not have to be a high value since spline interpolation provides a good
 * match. 1024 covers the most demanding precision requirement. For image
 * resizing this value can be set to a much lower value (e.g. FilterLen=10,
 * FilterFracs=64).
 */

template< int FilterLen, int FilterFracs >
class CDSPFracDelayFilterBank : public R8B_BASECLASS
{
public:
	CDSPFracDelayFilterBank()
	{
		R8BASSERT( FilterLen >= 4 );
		R8BASSERT(( FilterLen & 1 ) == 0 );
		R8BASSERT( FilterFracs > 0 );

		CDSPSincFilterGen :: CWindowFunc WinFunc;

		if( FilterLen >= 38 )
		{
			WinFunc = &CDSPSincFilterGen :: calcWindowVaneev38;
		}
		else
		if( FilterLen >= 32 )
		{
			WinFunc = &CDSPSincFilterGen :: calcWindowVaneev32;
		}
		else
		if( FilterLen >= 24 )
		{
			WinFunc = &CDSPSincFilterGen :: calcWindowVaneev24;
		}
		else
		if( FilterLen >= 20 )
		{
			WinFunc = &CDSPSincFilterGen :: calcWindowVaneev20;
		}
		else
		if( FilterLen >= 14 )
		{
			WinFunc = &CDSPSincFilterGen :: calcWindowVaneev14;
		}
		else
		{
			WinFunc = &CDSPSincFilterGen :: calcWindowVaneev10;
		}

		CDSPSincFilterGen sinc;
		sinc.Len2 = FilterLen / 2;

		double* p = Filters;
		int i;

		for( i = -1; i <= FilterFracs + 2; i++ )
		{
			sinc.FracDelay = (double) ( FilterFracs - i ) / FilterFracs;
			sinc.initFracVaneev();
			sinc.generateFrac( p, 4, WinFunc );
			normalizeFIRFilter( p, FilterLen, 1.0, 4 );
			p += FilterLen4;
		}

		// Calculate Hermite spline interpolation coefficients.

		p = Filters;
		const int FuncPos2 = FilterLen4;
		const int FuncPos3 = FilterLen4 * 2;
		const int FuncPos4 = FilterLen4 * 3;

		for( i = 0; i <= FilterFracs; i++ )
		{
			int l = FilterLen;

			while( l > 0 )
			{
				calcHermiteCoeffs( p, p[ 0 ], p[ FuncPos2 ], p[ FuncPos3 ],
					p[ FuncPos4 ]);

				p += 4;
				l--;
			}
		}
	}

	/**
	 * @param i Filter index, in the range 0 to FilterFracs, inclusive.
	 * @return Reference to the filter.
	 */

	const double& operator []( const int i ) const
	{
		R8BASSERT( i >= 0 && i <= FilterFracs );

		return( Filters[ i * FilterLen4 ]);
	}

private:
	static const int FilterLen4 = FilterLen << 2; ///< = FilterLen * 4.
		///<
	double Filters[ FilterLen * 4 * ( FilterFracs + 4 )]; ///< Table of
		///< fractional delay filters for all discrete fractional x = 0..1
		///< sample positions, and 3rd order Hermite spline interpolation
		///< coefficients.
		///<
};

/**
 * \brief Fractional delay filter-based interpolator class.
 *
 * Class implements the fractional delay interpolator. This implementation at
 * first puts the input signal into a ring buffer and then performs
 * interpolation. The interpolation is performed using sinc-based fractional
 * delay filters. These filters are contained in a bank, and for higher
 * precision they are interpolated between adjacent filters via 4-point
 * Hermite spline interpolation of 3rd order.
 *
 * While this class is used by the CDSPResampler class, it is not the only
 * possible application of this class. For example, this class with the
 * FilterLen template parameter as low as 16 can be used for real-time "pitch"
 * changes and resampling in soft-synths given a suitable oversampling and
 * dynamic low-pass filtering is performed on a prior step.
 *
 * To increase sample timing precision, this class uses "resettable counter"
 * approach. This gives less than "1 per 100 billion" sample timing error when
 * converting 44100 to 48000 sample rate.
 *
 * VERY IMPORTANT: the interpolation step should not exceed FilterLen/2+1
 * samples or the algorithm in its current form will fail. However, this
 * condition can be easily met if the input signal is downsampled first before
 * the interpolation is performed.
 *
 * Please see the r8bbase.cpp file for the example of definition of the static
 * constant "filter bank" object for the given template parameters of this
 * class.
 *
 * @param FilterLen Specifies the number of samples (taps) each fractional
 * delay filter should have. See the r8b::CDSPFracDelayFilterBank class for
 * more details.
 * @param FilterFracs The number of fractional delay positions to sample. See
 * the r8b::CDSPFracDelayFilterBank class for more details.
 * @param BufLenBits The length of the ring buffer, expressed as Nth power of
 * 2. This value can be reduced if it is known that short input buffers will
 * be passed to the interpolator. The minimum value of this parameter is 5,
 * and 1<<BufLenBits should be at least 3 times larger than the FilterLen.
 */

template< int FilterLen, int FilterFracs, int BufLenBits >
class CDSPFracInterpolator : public R8B_BASECLASS
{
	R8BNOCTOR( CDSPFracInterpolator );

public:
	/**
	 * Constructor initalizes the interpolator. It is important to call the
	 * getMaxOutLen() function afterwards to obtain the optimal output buffer
	 * length.
	 *
	 * @param aSrcSampleRate Source sample rate.
	 * @param aDstSampleRate Destination sample rate.
	 * @param InitFracPos Initial fractional position, in samples, in the
	 * range [0; 1). A non-zero value can be specified to remove the
	 * fractional delay introduced by a minimum-phase filter.
	 */

	CDSPFracInterpolator( const double aSrcSampleRate,
		const double aDstSampleRate, const double InitFracPos )
		: SrcSampleRate( aSrcSampleRate )
		, DstSampleRate( aDstSampleRate )
	{
		R8BASSERT( SrcSampleRate > 0.0 );
		R8BASSERT( DstSampleRate > 0.0 );
		R8BASSERT( InitFracPos >= 0.0 && InitFracPos < 1.0 );
		R8BASSERT( BufLenBits >= 5 );
		R8BASSERT(( 1 << BufLenBits ) >= FilterLen * 3 );

		BufLeft = 0;
		WritePos = 0;
		ReadPos = BufLen - FilterLenD2Minus1; // Set "read" position to
			// account for filter's latency at zero fractional delay which
			// equals to FilterLenD2Minus1.

		memset( Buf + ReadPos, 0, FilterLenD2Minus1 * sizeof( double ));

		InCounter = 0;
		InPosInt = 0;
		InPosFrac = InitFracPos;
		InPosShift = InitFracPos;
	}

	/**
	 * @return The number of samples that should be passed to *this object
	 * before the actual output starts.
	 */

	int getInLenBeforeOutStart() const
	{
		return( FilterLenD2Plus1 );
	}

	/**
	 * @param MaxInLen The number of samples planned to process at once, at
	 * most.
	 * @return The minimal length of the output buffer required when
	 * processing the "MaxInLen" number of input samples.
	 */

	int getMaxOutLen( const int MaxInLen ) const
	{
		R8BASSERT( MaxInLen >= 0 );

		return( (int) ceil( MaxInLen * DstSampleRate / SrcSampleRate ) + 1 );
	}

	/**
	 * Function changes the destination sample rate "on the fly". Note that
	 * the getMaxOutLen() function may needed to be called after calling this
	 * function as the maximal number of output samples produced by the
	 * interpolator depends on the destination sample rate.
	 *
	 * It can be a useful approach to construct *this object passing the
	 * maximal possible destination sample rate to the constructor, obtaining
	 * the getMaxOutLen() value and then setting the destination sample rate
	 * to whatever value is needed.
	 *
	 * @param NewDstSampleRate New destination sample rate.
	 */

	void setDstSampleRate( const double NewDstSampleRate )
	{
		R8BASSERT( DstSampleRate > 0.0 );

		DstSampleRate = NewDstSampleRate;
		InCounter = 0;
		InPosInt = 0;
		InPosShift = InPosFrac;
	}

	/**
	 * Function performs input sample stream interpolation.
	 *
	 * @param ip Input sample buffer.
	 * @param[out] op0 Output sample buffer, the capacity of this buffer
	 * should be equal to the value returned by the getMaxOutLen() function
	 * for the given "l". This buffer can be equal to "ip" only if the
	 * getMaxOutLen( l ) function's returned value is lesser than "l".
	 * @param l The number of samples available in the input sample buffer.
	 * @return The number of output samples written to the "op" buffer. The
	 * output latency is equal to FilterLenD2Plus1 samples.
	 */

	int process( const double* ip, double* const op0, int l )
	{
		R8BASSERT( l >= 0 );
		R8BASSERT( ip != op || l == 0 );

		double* op = op0;

		while( l > 0 )
		{
			// Add new input samples to both halves of the ring buffer.

			const int b = min( min( l, BufLen - WritePos ),
				BufLeftMax - BufLeft );

			double* const wp1 = Buf + WritePos;
			double* const wp2 = wp1 + BufLen;
			int i;

			for( i = 0; i < b; i++ )
			{
				wp1[ i ] = ip[ i ];
				wp2[ i ] = ip[ i ];
			}

			ip += b;
			WritePos = ( WritePos + b ) & BufLenMask;
			l -= b;
			BufLeft += b;

			// Produce as many output samples as possible.

			while( BufLeft >= FilterLenD2Plus1 )
			{
				double x = InPosFrac * FilterFracs;
				const int fti = (int) x; // Function table index.
				x -= fti; // Coefficient for interpolation between adjacent
					// fractional delay filters.
				const double x2 = x * x;
				const double x3 = x2 * x;
				const double* const ftp = &FilterBank[ fti ];
				const double* const rp = Buf + ReadPos;
				double s = 0.0;

				for( i = 0; i < FilterLen; i++ )
				{
					const int ii = i * 4;
					s += ( ftp[ ii ] + ftp[ ii + 1 ] * x +
						ftp[ ii + 2 ] * x2 + ftp[ ii + 3 ] * x3 ) * rp[ i ];
				}

				*op = s;
				op++;

				InCounter++;
				const double NextInPos =
					InCounter * SrcSampleRate / DstSampleRate + InPosShift;

				const int NextInPosInt = (int) NextInPos;
				const int PosIncr = NextInPosInt - InPosInt;
				InPosInt = NextInPosInt;
				InPosFrac = NextInPos - NextInPosInt;

				ReadPos = ( ReadPos + PosIncr ) & BufLenMask;
				BufLeft -= PosIncr;
			}
		}

		if( InCounter > 1000 )
		{
			// Reset the interpolation position counter to achieve a higher
			// sample timing precision.

			InCounter = 0;
			InPosInt = 0;
			InPosShift = InPosFrac;
		}

		return( (int) ( op - op0 ));
	}

private:
	static const int FilterLenD2Minus1 = ( FilterLen >> 1 ) - 1; ///< =
		///< ( FilterLen >> 1 ) - 1. This value also equals to filter's
		///< latency in samples (taps).
		///<
	static const int FilterLenD2Plus1 = ( FilterLen >> 1 ) + 1; ///< =
		///< ( FilterLen >> 1 ) + 1.
		///<
	static const int BufLen = 1 << BufLenBits; ///< The length of the ring
		///< buffer. The actual size is twice as large to allow "beyond max"
		///< positioning.
		///<
	static const int BufLenMask = BufLen - 1; ///< Mask used for quick buffer
		///< position wrapping.
		///<
	static const int BufLeftMax = BufLen - FilterLenD2Minus1; ///< The number
		///< of new samples that the ring buffer can hold at most. The
		///< remaining FilterLenD2Minus1 samples hold "previous" input samples
		///< for the filter.
		///<
	double Buf[ BufLen * 2 ]; ///< The ring buffer.
		///<
	double SrcSampleRate; ///< Source sample rate.
		///<
	double DstSampleRate; ///< Destination sample rate.
		///<
	int BufLeft; ///< The number of samples left in the buffer to process.
		///< When this value is below FilterLenD2Plus1, the interpolation
		///< cycle ends.
		///<
	int WritePos; ///< The current buffer write position. Incremented together
		///< with the BufLeft variable.
		///<
	int ReadPos; ///< The current buffer read position.
		///<
	int InCounter; ///< Interpolation step counter.
		///<
	int InPosInt; ///< Interpolation position (integer part).
		///<
	double InPosFrac; ///< Interpolation position (fractional part).
		///<
	double InPosShift; ///< Interpolation position fractional shift.
		///<
	static const CDSPFracDelayFilterBank< FilterLen, FilterFracs >
		FilterBank; ///< Static filter bank object.
		///<

	CDSPFracInterpolator()
	{
	}
};

} // namespace r8b

#endif // R8B_CDSPFRACINTERPOLATOR_INCLUDED
