//$ nocpp

/**
 * \file CDSPResampler.h
 *
 * \brief The master sample rate converter (resampler) class.
 *
 * This file includes the master sample rate converter (resampler) class that
 * combines all elements of this library into a single front-end class.
 *
 * r8brain-free-src Copyright (c) 2013 Aleksey Vaneev
 * See the "License.txt" file for license.
 */

#ifndef R8B_CDSPRESAMPLER_INCLUDED
#define R8B_CDSPRESAMPLER_INCLUDED

#include "CDSPBlockConvolver.h"
#include "CDSPFracInterpolator.h"

namespace r8b {

/**
 * \brief The master sample rate converter (resampler) class.
 *
 * This class can be considered the "master" sample rate converter (resampler)
 * class since it combines all functionality of this library into a single
 * front-end class to perform sample rate conversion to/from any sample rate,
 * including non-integer sample rates.
 *
 * Note that objects of this class can be constructed on the stack as it has a
 * small member data size.
 */

class CDSPResampler : public R8B_BASECLASS
{
	R8BNOCTOR( CDSPResampler );

public:
	/**
	 * Constructor initalizes the resampler object.
	 *
	 * Note that increasing the transition band and decreasing attenuation
	 * reduces the filter length, this in turn reduces the "input before
	 * output" delay. However, the filter length has only a minor influence on
	 * the overall resampling speed.
	 *
	 * In most cases ReqTransBand=3 and ReqAtten=96 can be used with good
	 * results. When resampling impulse responses, the ReqAtten=70 can be
	 * used. For full dynamic range 24-bit audio (e.g. for storage purposes)
	 * the ReqAtten=144 may be used, but as a rule ReqAtten above 100 is
	 * beyond human perception. When upsampling 88200 or 96000 audio to higher
	 * sample rates the ReqTransBand can be considerably increased, up to 20.
	 *
	 * It should be noted that ReqAtten specifies the minimal difference
	 * between the loudest input signal component and the produced aliasing
	 * artifacts during resampling. For example, if ReqAtten=96 was specified
	 * when performing 2x upsampling, the further analysis of the resulting
	 * signal may display high-frequency components which are quieter than the
	 * loudest part of the input signal by only 96 decibel. The high-frequency
	 * part won't become "magically" completely silent lean after resampling.
	 * You have to specify higher ReqAtten values if you need a totally clean
	 * high-frequency content. On the other hand, it may not be reasonable to
	 * have a high-frequency content cleaner than the input signal itself: if
	 * the input signal is 16-bit, setting ReqAtten to 144 will make its
	 * high-frequency content 24-bit, but the original part of the signal will
	 * remain 16-bit.
	 *
	 * @param SrcSampleRate Source signal sample rate.
	 * @param DstSampleRate Destination signal sample rate. The "power of 2"
	 * ratios between the source and destination sample rates force resampler
	 * to use several fast "power of 2" resampling steps, without using
	 * fractional interpolation at all. Note that the "power of 2" upsampling
	 * (but not downsampling) requires a lot of buffer memory: e.g. upsampling
	 * by a factor of 16 requires an intermediate buffer MaxInLen*(16+8)
	 * samples long. So, when doing the "power of 2" upsampling it is highly
	 * recommended to do it in small steps, e.g. no more than 256 samples at
	 * once (also set MaxInLen to 256).
	 * @param ReqTransBand Required transition band, in percent of the
	 * spectral space of the input signal (or the output signal if
	 * downsampling is performed) between filter's -3 dB point and the Nyquist
	 * frequency. The range is from CDSPFIRFilter::getLPMinTransBand() to
	 * CDSPFIRFilter::getLPMaxTransBand(), inclusive. Transition bands
	 * between 25% and 45% should be considered "coarse" as they have a higher
	 * error in both the -3 dB point position and attenuation.
	 * @param ReqAtten Required stop-band attenuation in decibel, in the range
	 * CDSPFIRFilter::getLPMinAtten() to CDSPFIRFilter::getLPMaxAtten(),
	 * inclusive.
	 * @param MaxInLen The maximal planned length of the input buffer (in
	 * samples) that will be passed to the resampler. The resampler relies on
	 * this value as it allocates intermediate buffers. Input buffers longer
	 * than this value should never be supplied to the resampler. Note that
	 * the resampler may use the input buffer itself for intermediate sample
	 * data storage.
	 * @see CDSPFIRFilterCache::getLPFilter()
	 */

	CDSPResampler( const double SrcSampleRate, const double DstSampleRate,
		const double ReqTransBand, const double ReqAtten, const int MaxInLen )
	{
		R8BASSERT( SrcSampleRate > 0.0 );
		R8BASSERT( DstSampleRate > 0.0 );
		R8BASSERT( MaxInLen >= 0 );

		if( SrcSampleRate == DstSampleRate )
		{
			ConvCount = 0;
			return;
		}

		int SrcSRMult;
		int SrcSRDiv = 1;
		int MaxOutLen = MaxInLen;
		int ConvBufCapacities[ 2 ];

		if( DstSampleRate * 2 > SrcSampleRate )
		{
			// Only a single convolver with 2X upsampling is required.

			SrcSRMult = 2;
			const double NormFreq = ( DstSampleRate > SrcSampleRate ? 0.5 :
				0.5 * DstSampleRate / SrcSampleRate );

			Convs[ 0 ] = new CDSPBlockConvolver(
				CDSPFIRFilterCache :: getLPFilter( NormFreq, ReqTransBand,
				ReqAtten ), rsmUpsample2X );

			ConvCount = 1;
			MaxOutLen = Convs[ 0 ] -> getMaxOutLen( MaxOutLen );
			ConvBufCapacities[ 0 ] = MaxOutLen;

			// Find if the destination to source sample rate ratio is
			// a "power of 2" value.

			int UseConvCount = 1;

			while( true )
			{
				const double TestSR = SrcSampleRate * ( 1 << UseConvCount );

				if( TestSR > DstSampleRate )
				{
					UseConvCount = 0; // Power of 2 not found.
					break;
				}

				if( TestSR == DstSampleRate )
				{
					break; // Power of 2 found.
				}

				UseConvCount++;
			}

			if( UseConvCount > 0 )
			{
				R8BASSERT( UseConvCount <= ConvCountMax );

				ConvBufCapacities[ 1 ] = 0;
				ConvCount = UseConvCount;
				int i;

				for( i = 1; i < UseConvCount; i++ )
				{
					const double tb = ( i >= 2 ? 48 : 34 );

					Convs[ i ] = new CDSPBlockConvolver(
						CDSPFIRFilterCache :: getLPFilter( 0.5, tb,
						ReqAtten ), rsmUpsample2X );

					MaxOutLen = Convs[ i ] -> getMaxOutLen( MaxOutLen );
					ConvBufCapacities[ i & 1 ] = MaxOutLen;
				}

				ConvBufs[ 0 ].alloc( ConvBufCapacities[ 0 ]);

				if( ConvBufCapacities[ 1 ] > 0 )
				{
					ConvBufs[ 1 ].alloc( ConvBufCapacities[ 1 ]);
				}

				return; // No interpolator is needed.
			}

			ConvBufs[ 0 ].alloc( ConvBufCapacities[ 0 ]);
		}
		else
		{
			SrcSRMult = 1;
			ConvBufCapacities[ 0 ] = 0;
			ConvCount = 0;
			const double CheckSR = DstSampleRate * 4;

			while( CheckSR * SrcSRDiv <= SrcSampleRate )
			{
				SrcSRDiv *= 2;

				// If downsampling is even deeper, use a less steep filter at
				// this step.

				const double tb =
					( CheckSR * SrcSRDiv <= SrcSampleRate ? 48 : 34 );

				Convs[ ConvCount ] = new CDSPBlockConvolver(
					CDSPFIRFilterCache :: getLPFilter( 0.5, tb, ReqAtten ),
					rsmDownsample2X );

				MaxOutLen = Convs[ ConvCount ] -> getMaxOutLen( MaxOutLen );
				ConvCount++;

				R8BASSERT( ConvCount < ConvCountMax );
			}

			const double NormFreq = DstSampleRate * SrcSRDiv / SrcSampleRate;
			const CDSPResamplingMode rsm =
				( NormFreq == 0.5 ? rsmDownsample2X : rsmNone );

			Convs[ ConvCount ] = new CDSPBlockConvolver(
				CDSPFIRFilterCache :: getLPFilter( NormFreq, ReqTransBand,
				ReqAtten ), rsm );

			MaxOutLen = Convs[ ConvCount ] -> getMaxOutLen( MaxOutLen );
			ConvCount++;

			if( rsm == rsmDownsample2X )
			{
				// No further interpolation is necessary.

				return;
			}
		}

		Interp = new CDSPFracInterpolator(
			SrcSampleRate * SrcSRMult / SrcSRDiv, DstSampleRate );

		MaxOutLen = Interp -> getMaxOutLen( MaxOutLen );

		if( MaxOutLen <= ConvBufCapacities[ 0 ])
		{
			InterpBuf = ConvBufs[ 0 ];
		}
		else
		if( MaxOutLen <= MaxInLen )
		{
			InterpBuf = NULL;
		}
		else
		{
			TmpBuf.alloc( MaxOutLen );
			InterpBuf = TmpBuf;
		}
	}

	/**
	 * Function performs sample rate conversion.
	 *
	 * If the source and destination sample rates are equal, the resampler
	 * will do nothing and will simply return the input buffer unchanged.
	 *
	 * @param ip0 Input buffer. This buffer may be used as output buffer by
	 * this function.
	 * @param l The number of samples available in the input buffer.
	 * @param[out] op0 This variable receives the pointer to the resampled
	 * data. This pointer may point to the address within the "ip" input
	 * buffer, or to *this object's internal buffer.
	 * @return The number of samples available in the "op0" output buffer.
	 */

	int process( double* const ip0, int l, double*& op0 )
	{
		R8BASSERT( l >= 0 );

		if( ConvCount == 0 )
		{
			op0 = ip0;
			return( l );
		}

		double* ip = ip0;
		double* op;
		int i;

		for( i = 0; i < ConvCount; i++ )
		{
			op = ( ConvBufs[ i & 1 ] == NULL ? ip0 : ConvBufs[ i & 1 ]);
			l = Convs[ i ] -> process( ip, op, l );
			ip = op;
		}

		if( Interp == NULL )
		{
			op0 = op;
			return( l );
		}

		op = ( InterpBuf == NULL ? ip0 : InterpBuf );
		op0 = op;

		return( Interp -> process( ip, op, l ));
	}

private:
	static const int ConvCountMax = 8; ///< 8 convolvers with the
		///< built-in 2x up- or downsampling is enough for 256x up- or
		///< downsampling.
		///<
	CPtrKeeper< CDSPBlockConvolver* > Convs[ ConvCountMax ]; ///< Convolvers.
		///<
	int ConvCount; ///< The number of objects defined in the Convs[] array.
		///< Equals to 0 if sample rate conversion is not needed.
		///<
	CPtrKeeper< CDSPFracInterpolator* > Interp; ///< Fractional interpolator
		///< object. Equals NULL if no fractional interpolation is required
		///< meaning the "power of 2" resampling performed or no resampling is
		///< performed at all.
		///<
	CFixedBuffer< double > ConvBufs[ 2 ]; ///< Intermediate convolution
		///< buffers to use, used only when at least 2x upsampling is
		///< performed. These buffers are used in flip-flop manner. If NULL
		///< then the input buffer will be used instead.
		///<
	CFixedBuffer< double > TmpBuf; ///< Additional output buffer, can be
		///< addressed by the InterpBuf pointer.
		///<
	double* InterpBuf; ///< Final output interpolation buffer to use. If NULL
		///< then the input buffer will be used instead. Otherwise this
		///< pointer points to either one of ConvBufs or TmpBuf.
		///<

	CDSPResampler()
	{
	}
};

} // namespace r8b

#endif // R8B_CDSPRESAMPLER_INCLUDED
