//$ nocpp

/**
 * \file CDSPResampler.h
 * \brief The master sample rate converter class.
 *
 * This file includes the master sample rate converter class that combines all
 * elements of this library into a single-end resampling engine.
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
 * \brief The master sample rate converter class.
 *
 * This class can be considered the "master" sample rate converter class since
 * it combines all functionality of this library into a single-end class to
 * perform sample rate conversion from/to any sample rate, including
 * non-integer sample rates.
 *
 * Note that, even if technically possible, resampling should not be performed
 * if the source and destination sample rates are equal since the low-pass
 * pass (reconstruction) filter is always applied to the signal.
 *
 * Note that objects of this class can be constructed on stack as it has a
 * small member data size.
 */

class CDSPResampler : public R8B_BASECLASS
{
	R8BNOCTOR( CDSPResampler );

public:
	/**
	 * Constructor initalizes the resampler object.
	 *
	 * @param SrcSampleRate Source signal sample rate.
	 * @param DstSampleRate Destination signal sample rate.
	 * @param ReqTransBand Required transition band, in percent of the
	 * spectral space of the input signal (or the output signal if
	 * downsampling is performed) between filter's -3 dB point and the Nyquist
	 * frequency. The range is from CDSPFIRFilter::getLPMinTransBand() to
	 * CDSPFIRFilter::getLPMaxTransBand(), inclusive.
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
		const double ReqTransBand, const double ReqAtten,
		const int MaxInLen )
	{
		int SrcSRMult;
		int SrcSRDiv;
		int MaxOutLen = MaxInLen;

		if( DstSampleRate * 2.0 > SrcSampleRate )
		{
			SrcSRMult = 2;
			SrcSRDiv = 1;

			// Only a single convolver with 2X upsampling is required.

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
