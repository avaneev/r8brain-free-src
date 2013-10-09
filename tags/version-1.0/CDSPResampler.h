//$ nocpp

/**
 * @file CDSPResampler.h
 *
 * @brief The master sample rate converter (resampler) class.
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
 * @brief The master sample rate converter (resampler) class.
 *
 * This class can be considered the "master" sample rate converter (resampler)
 * class since it combines all functionality of this library into a single
 * front-end class to perform sample rate conversion to/from any sample rate,
 * including non-integer sample rates.
 *
 * Note that objects of this class can be constructed on the stack as it has a
 * small member data size. The default template parameters of this class are
 * suited for 32-bit fixed point resampling (not to be confused with 32-bit
 * floating point with 24-bit mantissa) and so the ReqAtten is set to 192.66
 * by default.
 *
 * Use the CDSPResampler16 class for 16-bit resampling.
 *
 * Use the CDSPResampler24 class for 24-bit resampling (including 32-bit
 * floating point resampling).
 *
 * @param CInterpClass Interpolator class that should be used by the
 * resampler. The desired interpolation quality can be defined via the
 * template parameters of the interpolator class. See
 * r8b::CDSPFracInterpolator and r8b::CDSPFracDelayFilterBank for description
 * of the template parameters.
 */

template< class CInterpClass =
	CDSPFracInterpolator< R8B_FLTLEN, R8B_FLTFRACS, 9 > >
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
	 * In most cases ReqTransBand=2 and ReqAtten=97 can be used with good
	 * results: this produces full dynamic range 16-bit audio. For full
	 * dynamic range 24-bit audio (e.g. for storage purposes) the ReqAtten=145
	 * may be used. For comparison purposes the ReqAtten=193 may be used.
	 * However, as a rule ReqAtten above 100 is beyond human perception. When
	 * upsampling 88200 or 96000 audio to higher sample rates the ReqTransBand
	 * can be considerably increased, up to 20.
	 *
	 * It should be noted that the ReqAtten specifies the minimal difference
	 * between the loudest input signal component and the produced aliasing
	 * artifacts during resampling. For example, if ReqAtten=100 was specified
	 * when performing 2x upsampling, the analysis of the resulting signal may
	 * display high-frequency components which are quieter than the loudest
	 * part of the input signal by only 100 decibel meaning the high-frequency
	 * part did not become "magically" completely silent after resampling. You
	 * have to specify a higher ReqAtten value if you need a totally clean
	 * high-frequency content. On the other hand, it may not be reasonable to
	 * have a high-frequency content cleaner than the input signal itself: if
	 * the input signal is 16-bit, setting ReqAtten to 150 will make its
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
	 * once (also set the MaxInLen to 256).
	 * @param MaxInLen The maximal planned length of the input buffer (in
	 * samples) that will be passed to the resampler. The resampler relies on
	 * this value as it allocates intermediate buffers. Input buffers longer
	 * than this value should never be supplied to the resampler. Note that
	 * the resampler may use the input buffer itself for intermediate sample
	 * data storage.
	 * @param ReqTransBand Required transition band, in percent of the
	 * spectral space of the input signal (or the output signal if
	 * downsampling is performed) between filter's -3 dB point and the Nyquist
	 * frequency. The range is from CDSPFIRFilter::getLPMinTransBand() to
	 * CDSPFIRFilter::getLPMaxTransBand(), inclusive.
	 * @param ReqAtten Required stop-band attenuation in decibel, in the range
	 * CDSPFIRFilter::getLPMinAtten() to CDSPFIRFilter::getLPMaxAtten(),
	 * inclusive. The actual attenuation may be 0.40-4.41 dB higher.
	 * @param ReqPhase Required filter's phase response. Note that this
	 * setting does not affect interpolator's phase response which is always
	 * linear-phase. Also note that if the "power of 2" resampling was engaged
	 * by the resampler together with the minimum-phase response, the audio
	 * stream may become fractionally delayed by up to 1 sample, depending on
	 * the minimum-phase filter's actual fractional delay. If the output
	 * stream should always start at "time zero" offset with minimum-phase
	 * filters the UsePower2 should be set to "false". Linear-phase filters
	 * do not have fractional delay.
	 * @param UsePower2 "True" if the "power of 2" resampling optimization
	 * should be used when possible.
	 * @see CDSPFIRFilterCache::getLPFilter()
	 */

	CDSPResampler( const double SrcSampleRate, const double DstSampleRate,
		const int MaxInLen, const double ReqTransBand = 2.0,
		const double ReqAtten = 192.66,
		const EDSPFilterPhaseResponse ReqPhase = fprLinearPhase,
		const bool UsePower2 = true )
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
		double PrevLatencyFrac = 0.0;

		if( DstSampleRate * 2 > SrcSampleRate )
		{
			// Only a single convolver with 2X upsampling is required.

			SrcSRMult = 2;
			const double NormFreq = ( DstSampleRate > SrcSampleRate ? 0.5 :
				0.5 * DstSampleRate / SrcSampleRate );

			Convs[ 0 ] = new CDSPBlockConvolver(
				CDSPFIRFilterCache :: getLPFilter( NormFreq, ReqTransBand,
				ReqAtten, ReqPhase ), rsmUpsample2X, 0.0 );

			ConvCount = 1;
			MaxOutLen = Convs[ 0 ] -> getMaxOutLen( MaxOutLen );
			ConvBufCapacities[ 0 ] = MaxOutLen;
			PrevLatencyFrac = Convs[ 0 ] -> getLatencyFrac();

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

			if( UsePower2 && UseConvCount > 0 )
			{
				R8BASSERT( UseConvCount <= ConvCountMax );

				ConvBufCapacities[ 1 ] = 0;
				ConvCount = UseConvCount;
				int i;

				for( i = 1; i < UseConvCount; i++ )
				{
					const double tb = ( i >= 2 ? 45.0 : 34.0 );

					Convs[ i ] = new CDSPBlockConvolver(
						CDSPFIRFilterCache :: getLPFilter( 0.5, tb, ReqAtten,
						ReqPhase ), rsmUpsample2X, PrevLatencyFrac );

					MaxOutLen = Convs[ i ] -> getMaxOutLen( MaxOutLen );
					ConvBufCapacities[ i & 1 ] = MaxOutLen;
					PrevLatencyFrac = Convs[ i ] -> getLatencyFrac();
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
					( CheckSR * SrcSRDiv <= SrcSampleRate ? 45.0 : 34.0 );

				Convs[ ConvCount ] = new CDSPBlockConvolver(
					CDSPFIRFilterCache :: getLPFilter( 0.5, tb, ReqAtten,
					ReqPhase ), rsmDownsample2X, PrevLatencyFrac );

				MaxOutLen = Convs[ ConvCount ] -> getMaxOutLen( MaxOutLen );
				PrevLatencyFrac = Convs[ ConvCount ] -> getLatencyFrac();
				ConvCount++;

				R8BASSERT( ConvCount < ConvCountMax );
			}

			const double NormFreq = DstSampleRate * SrcSRDiv / SrcSampleRate;
			const EDSPResamplingMode rsm =
				( UsePower2 && NormFreq == 0.5 ? rsmDownsample2X : rsmNone );

			Convs[ ConvCount ] = new CDSPBlockConvolver(
				CDSPFIRFilterCache :: getLPFilter( NormFreq, ReqTransBand,
				ReqAtten, ReqPhase ), rsm, PrevLatencyFrac );

			MaxOutLen = Convs[ ConvCount ] -> getMaxOutLen( MaxOutLen );
			PrevLatencyFrac = Convs[ ConvCount ] -> getLatencyFrac();
			ConvCount++;

			if( rsm == rsmDownsample2X )
			{
				return; // No interpolator is needed.
			}
		}

		Interp = new CInterpClass( SrcSampleRate * SrcSRMult / SrcSRDiv,
			DstSampleRate, PrevLatencyFrac );

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
	 * Function clears (resets) the state of *this object and returns it to
	 * the state after construction. All input data accumulated in the
	 * internal buffer so far will be discarded.
	 *
	 * This function makes it possible to use *this object for converting
	 * separate streams from the same source sample rate to the same
	 * destination sample rate without reconstructing the object. It is more
	 * efficient to clear the state of the resampler object than to destroy it
	 * and create a new object.
	 */

	void clear()
	{
		int i;

		for( i = 0; i < ConvCount; i++ )
		{
			Convs[ i ] -> clear();
		}

		if( Interp != NULL )
		{
			Interp -> clear();
		}
	}

	/**
	 * @return Interpolator object used by *this resampler. This function
	 * returns NULL if no interpolator is in use.
	 */

	CInterpClass* getInterpolator() const
	{
		return( Interp );
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
	 * data. This pointer may point to the address within the "ip0" input
	 * buffer, or to *this object's internal buffer. In real-time applications
	 * it is suggested to pass this pointer to the next output audio block and
	 * consume any data left from the previous output audio block first before
	 * calling the process() function again.
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
	CPtrKeeper< CInterpClass* > Interp; ///< Fractional interpolator object.
		///< Equals NULL if no fractional interpolation is required meaning
		///< the "power of 2" resampling is performed or no resampling is
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
		///< pointer points to either ConvBufs or TmpBuf.
		///<

	CDSPResampler()
	{
	}
};

/**
 * @brief The resampler class for 16-bit resampling.
 *
 * This class defines resampling parameters suitable for 16-bit resampling,
 * using linear-phase low-pass filter. See the r8b::CDSPResampler class for
 * details.
 */

class CDSPResampler16 :
	public CDSPResampler< CDSPFracInterpolator< 14, 67, 9 > >
{
public:
	/**
	 * Constructor initializes the 16-bit resampler. See the
	 * r8b::CDSPResampler class for details.
	 *
	 * @param SrcSampleRate Source signal sample rate.
	 * @param DstSampleRate Destination signal sample rate.
	 * @param MaxInLen The maximal planned length of the input buffer (in
	 * samples) that will be passed to the resampler.
	 * @param ReqTransBand Required transition band, in percent.
	 */

	CDSPResampler16( const double SrcSampleRate, const double DstSampleRate,
		const int MaxInLen, const double ReqTransBand = 2.0 )
		: CDSPResampler< CDSPFracInterpolator< 14, 67, 9 > >( SrcSampleRate,
			DstSampleRate, MaxInLen, ReqTransBand, 96.33, fprLinearPhase,
			true )
	{
	}
};

/**
 * @brief The resampler class for 24-bit resampling.
 *
 * This class defines resampling parameters suitable for 24-bit resampling
 * (including 32-bit floating point resampling), using linear-phase low-pass
 * filter. See the r8b::CDSPResampler class for details.
 */

class CDSPResampler24 :
	public CDSPResampler< CDSPFracInterpolator< 20, 211, 9 > >
{
public:
	/**
	 * Constructor initializes the 24-bit resampler (including 32-bit floating
	 * point). See the r8b::CDSPResampler class for details.
	 *
	 * @param SrcSampleRate Source signal sample rate.
	 * @param DstSampleRate Destination signal sample rate.
	 * @param MaxInLen The maximal planned length of the input buffer (in
	 * samples) that will be passed to the resampler.
	 * @param ReqTransBand Required transition band, in percent.
	 */

	CDSPResampler24( const double SrcSampleRate, const double DstSampleRate,
		const int MaxInLen, const double ReqTransBand = 2.0 )
		: CDSPResampler< CDSPFracInterpolator< 20, 211, 9 > >( SrcSampleRate,
			DstSampleRate, MaxInLen, ReqTransBand, 144.50, fprLinearPhase,
			true )
	{
	}
};

} // namespace r8b

#endif // R8B_CDSPRESAMPLER_INCLUDED