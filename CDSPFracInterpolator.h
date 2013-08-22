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
 * of 3rd order. The filters are windowed by the Hann^6 windowing function.
 * The N and Fracs parameters can be varied freely without breaking the
 * resampler.
 *
 * @param N Specifies the number of samples (taps) each filter should have.
 * This must be an even value, the minimal value for N is 4. To achieve a
 * higher resampling precision, the oversampling should be used in the first
 * place instead of using a higher N value.
 * @param Fracs The number of fractional delay positions to sample. This do
 * not have to be a high value since spline interpolation provides a good
 * match. 1024 covers the most demanding precision requirement. For image
 * resizing this value can be set to a much lower value (e.g. N=10, Fracs=64).
 */

template< int N, int Fracs >
class CDSPFracDelayFilterBank : public R8B_BASECLASS
{
protected:
	static const int FuncLen = N; ///< Fractional delay function length in
		///< samples (taps).
		///<
	static const int FuncLen2 = FuncLen << 1; ///< = FuncLen * 2.
		///<
	static const int FuncLen4 = FuncLen << 2; ///< = FuncLen * 4.
		///<
	static const int FuncLenD2 = FuncLen >> 1; ///< = FuncLen / 2.
		///<
	static const int FuncLenD2Minus1 = FuncLenD2 - 1; ///< = FuncLenD2 - 1.
		///< This value also equals to filter's latency in samples (taps).
		///<
	static const int FuncLenD2Plus1 = FuncLenD2 + 1; ///< = FuncLenD2 + 1.
		///<
	static const int FuncFrac = Fracs; ///< The number of fractional sample
		///< positions to use.
		///<

	/**
	 * \brief Static fractional delay filter table class.
	 *
	 * Static fractional delay filter function table storage and
	 * initialization class.
	 */

	struct CFuncTable
	{
		double FuncTable[ FuncLen4 * ( FuncFrac + 4 )]; ///< Table of
			///< fractional delay functions for all discrete fractional x =
			///< 0..1 sample positions, and 3rd order Hermite spline
			///< interpolation coefficients.
			///<

		double& operator []( const int i )
		{
			return( FuncTable[ i ]);
		}

		const double& operator []( const int i ) const
		{
			return( FuncTable[ i ]);
		}

		CFuncTable()
		{
			R8BASSERT( N >= 4 );
			R8BASSERT( Fracs > 0 );

			CDSPSincFilterGen sinc;
			sinc.Len2 = FuncLenD2;

			double* p = (double*) FuncTable;
			int i;

			for( i = -1; i <= FuncFrac + 2; i++ )
			{
				sinc.FracDelay = (double) ( FuncFrac - i ) / FuncFrac;
				sinc.initFracVaneev();
				sinc.generateFrac( p, 4,
					&CDSPSincFilterGen :: calcWindowVaneev );

				normalizeFIRFilter( p, FuncLen, 1.0, 4 );
				p += FuncLen4;
			}

			// Calculate Hermite spline interpolation coefficients.

			p = FuncTable;
			const int FuncPos2 = FuncLen4;
			const int FuncPos3 = FuncLen4 * 2;
			const int FuncPos4 = FuncLen4 * 3;

			for( i = 0; i <= FuncFrac; i++ )
			{
				int l = FuncLen;

				while( l > 0 )
				{
					calcHermiteCoeffs( p, p[ 0 ], p[ FuncPos2 ],
						p[ FuncPos3 ], p[ FuncPos4 ]);

					p += 4;
					l--;
				}
			}
		}
	};

	static const CFuncTable FuncTable; ///< Static function table object.
		///<
};

/**
 * \brief Fractional delay filter-based interpolator class.
 *
 * Class implements the fractional delay interpolator. This implementation at
 * first puts the input signal into a ring buffer and then performs
 * interpolation. The interpolation is performed using short (40 taps)
 * sinc-based fractional delay filters. These filters are contained in a bank,
 * and for higher precision they are interpolated between adjacent filters via
 * 4-point Hermite spline interpolation of 3rd order.
 *
 * To increase sample timing precision, this class uses "resettable counter"
 * approach. This gives less than "1 per 100 billion" sample timing error when
 * converting 44100 to 48000 sample rate.
 *
 * VERY IMPORTANT: the interpolation step should not exceed FuncLenD2Plus1
 * samples or the algorithm in its current form will fail. However, this
 * condition can be easily met if the input signal is downsampled first before
 * the interpolation is performed.
 */

class CDSPFracInterpolator : public CDSPFracDelayFilterBank< 38, 1280 >
{
	R8BNOCTOR( CDSPFracInterpolator );

public:
	/**
	 * Constructor initalizes the interpolator. It is important to call the
	 * getMaxOutLen() function afterwards to obtain the optimal output buffer
	 * length.
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
	 * @return The number of samples that should be passed to *this object
	 * before the actual output starts.
	 */

	int getInLenBeforeOutStart() const
	{
		return( FuncLenD2Plus1 );
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
	 * Function performs input sample stream interpolation.
	 *
	 * @param ip Input sample buffer.
	 * @param[out] op Output sample buffer, the capacity of this buffer should
	 * be equal to the value returned by the getMaxOutLen() function for the
	 * given "l". This buffer can be equal to "ip" only if the
	 * getMaxOutLen( l ) function's returned value is lesser than "l".
	 * @param l The number of samples available in the input sample buffer.
	 * @return The number of output samples written to the "op" buffer. The
	 * output latency is equal to FuncLenD2Plus1 samples (5 samples with
	 * 8-sample filter).
	 */

	int process( const double* ip, double* op, int l )
	{
		R8BASSERT( l >= 0 );
		R8BASSERT( ip != op || l == 0 );

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
				const double x2 = x * x;
				const double x3 = x2 * x;
				const double* const ftp = &FuncTable[ fti * FuncLen4 ];
				const double* const rp = Buf + ReadPos;
				double s = 0.0;

				for( i = 0; i < FuncLen; i++ )
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
		///<
	double SrcSampleRate; ///< Source sample rate.
		///<
	double DstSampleRate; ///< Destination sample rate.
		///<
	int BufLeft; ///< The number of samples left in the buffer to process.
		///< When this value is below FuncLenD2Plus1, the interpolation cycle
		///< ends.
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

	CDSPFracInterpolator()
	{
	}
};

} // namespace r8b

#endif // R8B_CDSPFRACINTERPOLATOR_INCLUDED
