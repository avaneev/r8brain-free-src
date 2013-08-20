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
 * Note that, even if technically possible, resampling should not be performed
 * if the source and destination sample rates are equal since the low-pass
 * pass (reconstruction) filter is always applied to the signal.
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
	 * @param SrcSampleRate Source signal sample rate.
	 * @param DstSampleRate Destination signal sample rate.
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

		int SrcSRMult;
		int SrcSRDiv;
		int MaxOutLen = MaxInLen;

		if( DstSampleRate * 2.0 > SrcSampleRate )
		{
			// Only a single convolver with 2X upsampling is required.

			SrcSRMult = 2;
			SrcSRDiv = 1;

			const double NormFreq = ( DstSampleRate >= SrcSampleRate ? 0.5 :
				0.5 * DstSampleRate / SrcSampleRate );

			ConvCount = 1;
			Convs[ 0 ] = new CDSPBlockConvolver(
				CDSPFIRFilterCache :: getLPFilter( NormFreq, ReqTransBand,
				ReqAtten ), CDSPResamplingMode :: rsUpsample2X );

			ConvBuffer.alloc( MaxInLen * 2 );
			MaxOutLen = Convs[ 0 ] -> getMaxOutLen( MaxOutLen );
		}
		else
		{
			SrcSRMult = 1;
		}

		Interp = new CDSPFracInterpolator( SrcSampleRate * SrcSRMult /
			SrcSRDiv, DstSampleRate );

		MaxOutLen = Interp -> getMaxOutLen( MaxOutLen );

		if( ConvBuffer != NULL && MaxOutLen <= MaxInLen * 2 )
		{
			InterpBuffer = ConvBuffer;
		}
		else
		if( MaxOutLen <= MaxInLen )
		{
			InterpBuffer = NULL;
		}
		else
		{
			TmpBuffer.alloc( MaxOutLen );
			InterpBuffer = TmpBuffer;
		}
	}

	/**
	 * Function performs sample rate conversion.
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

		double* const op1 = ( ConvBuffer == NULL ? ip0 : ConvBuffer );
		double* ip = ip0;
		double* op;
		int i;

		for( i = 0; i < ConvCount; i++ )
		{
			op = op1;
			l = Convs[ i ] -> process( ip, op, l );
			ip = op;
		}

		op = ( InterpBuffer == NULL ? ip0 : InterpBuffer );
		op0 = op;

		return( Interp -> process( ip, op, l ));
	}

private:
	CPtrKeeper< CDSPBlockConvolver* > Convs[ 8 ]; ///< 8 convolvers with the
		///< built-in 2x downsampling is enough for 256x downsampling.
		///<
	int ConvCount; ///< The number of objects defined in the Convs[] array.
		///<
	CPtrKeeper< CDSPFracInterpolator* > Interp; ///< Fractional interpolator
		///< object.
		///<
	CFixedBuffer< double > ConvBuffer; ///< Intermediate convolution buffer to
		///< use, used only when 2x upsampling is performed. If NULL then the
		///< input buffer will be used instead.
		///<
	CFixedBuffer< double > TmpBuffer; ///< Additional output buffer.
		///<
	double* InterpBuffer; ///< Final output interpolation buffer to use.
		///< If NULL then then input buffer will be used instead. Otherwise
		///< this pointer points to either ConvBuffer or TmpBuffer.
		///<

	CDSPResampler()
	{
	}
};

} // namespace r8b

#endif // R8B_CDSPRESAMPLER_INCLUDED
