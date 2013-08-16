//$ nocpp

/**
 * \file CDSPFracInterpolator.h
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
 * fractional delay filters. The filters are windowed by the Blackman
 * windowing function.
 *
 * @param N Specifies the number of samples (taps) each filter
 * should have. This must be an even value. Instead of using a high N value
 * for higher precision the oversampling should be used in the first place.
 * @param Fracs The number of fractional delay positions to sample. Using
 * higher values does not increase the overall precision considerably.
 */

template< int N, int Fracs >
class CDSPFracDelayFilterBank : public R8B_BASECLASS
{
public:
	CDSPFracDelayFilterBank()
	{
		CDSPSincFilterGen sinc;
		sinc.Len2 = FuncLenD2;
		sinc.FracDelay = 1.0;

		double* p = &FuncTable[ 0 ];
		sinc.initFrac();
		sinc.generateFrac( p, 2 );
		normalizeFIRFilter( p, FuncLen, 1.0, 2 );
		p += FuncLen2;

		int i;

		for( i = 1; i <= FuncFrac; i++ )
		{
			sinc.FracDelay = (double) ( FuncFrac - i ) / FuncFrac;
			sinc.initFrac();
			sinc.generateFrac( p, 2 );
			normalizeFIRFilter( p, FuncLen, 1.0, 2 );

			int z;

			for( z = 0; z < FuncLen; z++ )
			{
				// Calculate delta between adjacent fractional delay filters.

				p[ 1 - FuncLen2 ] = p[ 0 ] - p[ -FuncLen2 ];
				p += 2;
			}
		}
	}

protected:
	static const int FuncLen = N; ///< Fractional delay function length in
		///< samples (taps).
		///<
	static const int FuncLen2 = FuncLen << 1; ///< = FuncLen * 2.
	static const int FuncLenD2 = FuncLen >> 1; ///< = FuncLen / 2.
	static const int FuncLenD2Minus1 = FuncLenD2 - 1; ///< = FuncLenD2 - 1.
		///< This value also equals to filter's latency in samples (taps).
		///<
	static const int FuncLenD2Plus1 = FuncLenD2 + 1; ///< = FuncLenD2 + 1.
	static const int FuncFrac = Fracs; ///< The number of fractional sample
		///< positions to use.
		///<
	double FuncTable[ FuncLen2 * ( FuncFrac + 1 )]; ///< Table of fractional
		///< delay functions for all discrete fractional x = 0..1 sample
		///< positions, and deltas between adjacent functions.
		///<
};

/**
 * \brief Fractional delay interpolator class.
 *
 * Class implements the fractional delay interpolator. This implementation at
 * first puts the input signal into a ring buffer and then performs
 * interpolation. The interpolation is performed using very short (8 taps)
 * sinc-based fractional delay filters. These filters are contained in a bank,
 * and for the sake of higher precision two close filters are crossfaded to
 * achieve an even higher precision.
 *
 * To increase sample timing precision, this class uses "resettable counter"
 * approach. This gives less than "1 per 100 billion" sample timing error when
 * converting 44100 to 48000 sample rate.
 *
 * VERY IMPORTANT: the interpolation step should not exceed FuncLenD2Plus1
 * samples or the algorithm in its current form will fail. However, this
 * condition can be easily met if the input signal gets downsampled first
 * before the interpolation is performed.
 */

class CDSPFracInterpolator : public CDSPFracDelayFilterBank< 8, 128 >
{
	R8BNOCTOR( CDSPFracInterpolator );

public:
	/**
	 * Constructor initalizes the interpolator. It is important to call the
	 * getOutputBufLen() function afterwards to obtain the optimal output
	 * buffer length.
	 *
	 * @param aSrcSampleRate Source signal sample rate.
	 * @param aDstSampleRate Destination signal sample rate.
	 */

	CDSPFracInterpolator( const double aSrcSampleRate,
		const double aDstSampleRate )
		: SrcSampleRate( aSrcSampleRate )
		, DstSampleRate( aDstSampleRate )
	{
		R8BASSERT( SrcSampleRate > 0.0 );
		R8BASSERT( DstSampleRate > 0.0 );

		BufLeft = 0;
		WritePos = 0;
		ReadPos = BufLen - FuncLenD2Minus1; // Set "read" position to account
			// for filter's latency at zero fractional delay which equals to
			// FuncLenD2Minus1.

		memset( Buf + ReadPos, 0, FuncLenD2Minus1 * sizeof( double ));

		InCounter = 0;
		InPosInt = 0;
		InPosFrac = 0.0;
		InPosShift = 0.0;
	}

	/**
	 * @param l The number of samples at most planned to process at once.
	 * @return The minimal length of the output buffer required when
	 * processing the "l" number of input samples.
	 */

	int getOutputBufLen( const int l ) const
	{
		R8BASSERT( l >= 0 );

		return( (int) ceil( l * DstSampleRate / SrcSampleRate ) + 1 );
	}

	/**
	 * Function performs input sample stream interpolation.
	 *
	 * @param ip Input sample buffer.
	 * @param[out] op Output sample buffer, the capacity of this buffer should
	 * be equal to the value returned by the getOutputBufLen() function for
	 * the given "l". This buffer should not be equal to "ip".
	 * @param l The number of samples available in the input sample buffer.
	 * @return The number of output samples written to the "op" buffer.
	 */

	int process( const double* ip, double* op, int l )
	{
		R8BASSERT( l >= 0 );
		R8BASSERT( ip != NULL || l == 0 );
		R8BASSERT( op != NULL );

		double* const op0 = op;

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

			while( BufLeft >= FuncLenD2Plus1 )
			{
				double x = InPosFrac * FuncFrac;
				const int fti = (int) x; // Function table index.
				x -= fti; // Coefficient for cross-fade between adjacent
					// fractional delay filters.
				const double* const ftp = &FuncTable[ fti * FuncLen2 ];
				const double* const rp = Buf + ReadPos;

/*				*op = ( ftp[ 0 ] + ftp[ 1 ] * x ) * rp[ 0 ] +
					( ftp[  2 ] + ftp[  3 ] * x ) * rp[ 1 ] +
					( ftp[  4 ] + ftp[  5 ] * x ) * rp[ 2 ] +
					( ftp[  6 ] + ftp[  7 ] * x ) * rp[ 3 ] +
					( ftp[  8 ] + ftp[  9 ] * x ) * rp[ 4 ] +
					( ftp[ 10 ] + ftp[ 11 ] * x ) * rp[ 5 ] +
					( ftp[ 12 ] + ftp[ 13 ] * x ) * rp[ 6 ] +
					( ftp[ 14 ] + ftp[ 15 ] * x ) * rp[ 7 ];*/

				double s = 0.0;

				for( i = 0; i < FuncLen; i++ )
				{
					s += ( ftp[ i * 2 ] + ftp[ i * 2 + 1 ] * x ) * rp[ i ];
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
	static const int BufLenBits = 9; ///< The length of the ring buffer,
		///< expressed as Nth power of 2.
		///<
	static const int BufLen = 1 << BufLenBits; ///< The length of the ring
		///< buffer. The actual size is twice as large to allow "beyond max"
		///< positioning.
		///<
	static const int BufLenMask = BufLen - 1; ///< Mask used for quick buffer
		///< position wrapping.
		///<
	static const int BufLeftMax = BufLen - FuncLenD2Minus1; ///< The number of
		///< new samples that the ring buffer can hold at most. The remaining
		///< FuncLenD2Minus1 samples hold "previous" input samples for the
		///< filter.
		///<
	double Buf[ BufLen * 2 ]; ///< The ring buffer.
	double SrcSampleRate; ///< Source sample rate.
	double DstSampleRate; ///< Destination sample rate.
	int BufLeft; ///< The number of samples left in the buffer to process.
		///< When this value is below FuncLenD2Plus1, the interpolation cycle
		///< ends.
		///<
	int WritePos; ///< The current buffer write position. Incremented together
		///< with the BufLeft variable.
		///<
	int ReadPos; ///< The current buffer read position.
	int InCounter; ///< Interpolation position counter.
	int InPosInt; ///< Interpolation position (integer part).
	double InPosFrac; ///< Interpolation position (fractional part).
	double InPosShift; ///< Interpolation position fractional shift.

	CDSPFracInterpolator()
	{
	}
};

} // namespace r8b

#endif // R8B_CDSPFRACINTERPOLATOR_INCLUDED
