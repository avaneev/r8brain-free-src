//$ nocpp

/**
 * \file CDSPFIRFilter.h
 * \brief FIR filter generator and filter cache classes.
 *
 * This file includes low-pass FIR filter generator and filter cache.
 *
 * r8brain-free-src Copyright (c) 2013 Aleksey Vaneev
 * See the "License.txt" file for license.
 */

#ifndef R8B_CDSPFIRFILTER_INCLUDED
#define R8B_CDSPFIRFILTER_INCLUDED

#include "CDSPSincFilterGen.h"
#include "CDSPRealFFT.h"

namespace r8b {

/**
 * \brief Calculation and storage class for FIR filters.
 *
 * Class that implements calculation and storing of a FIR filter (currently
 * contains low-pass filter calculation routines designed for sample rate
 * conversion. Objects of this class cannot be created directly, but can be
 * obtained via the CDSPFilterCache object.
 */

class CDSPFIRFilter : public R8B_BASECLASS
{
	friend class CDSPBlockConvolver;
	friend class CDSPFIRFilterCache;

public:
	/**
	 * @return The minimal allowed low-pass filter's transition band, in
	 * percent.
	 */

	static double getLPMinTransBand()
	{
		return( 0.5 );
	}

	/**
	 * @return The maximal allowed low-pass filter's transition band, in
	 * percent.
	 */

	static double getLPMaxTransBand()
	{
		return( 25.0 );
	}

	/**
	 * @return The minimal allowed low-pass filter's stop-band attenuation, in
	 * decibel.
	 */

	static double getLPMinAtten()
	{
		return( 38.0 );
	}

	/**
	 * @return The maximal allowed low-pass filter's stop-band attenuation, in
	 * decibel.
	 */

	static double getLPMaxAtten()
	{
		return( 220.0 );
	}

private:
	double ReqTransBand; ///< Required transition band in percent, as passed
		/// by the user.
	double ReqAtten; ///< Required stop-band attenuation in decibel, as passed
		/// by the user (positive value).
	double ReqNormFreq; ///< Required normalized frequency, 0 to 1 inclusive.
	int Latency; ///< Filter's latency in samples.
	int BlockSizeBits; ///< Block size used to store filter response,
		/// expressed as Nth power of 2. This value is used directly by the
		/// convolver.
	CFixedBuffer< double > KernelBlock; ///< Filter response buffer, size
		/// equals to 1 << ( BlockSizeBits + 1 ). Second part of the buffer
		/// contains zero-padding to allow alias-free convolution.

	/**
	 * Function builds filter kernel based on the ReqTransBand, ReqAtten and
	 * ReqNormFreq parameters.
	 */

	void buildLPFilter()
	{
		const double tb = ReqTransBand * 0.01;
		double atten = -ReqAtten;

		// Apply correction.

		if( atten >= -74.035088 )
		{
			atten -= 0.499762495989797 + 5.14372798398936e-5 *
				sqr( atten ) + 1072.6794503315 * tb / ( 1097.57988452953 -
				sqr( atten )) + 0.055976849615706 * cos( 0.355505067732265 *
				atten ) * cos( 0.228340476535502 * atten ) -
				0.388481865573685  *cos( 0.355505067732265 * atten ) -
				0.034075657453079 * atten * cos( 0.228340476535502 * atten );
		}
		else
		if( atten >= -98.666667 )
		{
			atten -= 0.323720648767925 + 0.573420106154753 * sin(
				-0.416282878764675 * atten ) + 0.220918140226294 * atan( sin(
				0.540772530335019 * atten )) - tb * atan2( 1.18464210151739 -
				sin( -0.346828095957449 * atten ), 4.98198649859924 * tb );
		}
		else
		{
			atten -= 0.71361770778327 * sin( 0.080859249933107 * atten ) *
				atan2( 0.186753898792276, 8.55050852098357 +
				0.080859249933107 * atten ) - 0.40224810474143 -
				0.72626322840629 * sin( -0.12386247667656 * atten );
		}

		const double hl = ( -10.2934576604145 - 0.285144433965609 * atten ) /
			tb + -211.072141982443 / ( tb * atten - 12.3281863268897 * tb -
			0.350469717177663 * tb * atten * atan( sin( 4.22528055809605 -
			0.0953582696744781 * atten )));

		const double hltb = hl * tb;
		const double pwr = 0.1716055279404 + 0.116379153860107 * hltb +
			0.0183229986637533 * sin( 0.128970625216538 * hltb ) -
			0.756472509630929 * gauss( gauss( sin( 6.1588739570383 +
			0.128970625216538 * hltb ) - 0.137337051614019 * hltb ) -
			0.326334088682917 );

		const double fo = atan2( 1.37248462316218 * asinh( 0.676453431881951 +
			pwr ), hl ) + gauss( cos( 2.38417547027863 * pwr ) /
			( 0.532351950498439 * hl + 6.94666497679782 * exp( sqr( pwr )))) -
			0.976759641634279 * tb;

		CDSPSincFilterGen sinc;
		sinc.Len2 = 0.25 * hl / NormFreq;
		sinc.Freq1 = 0.0;
		sinc.Freq2 = M_PI * fo * NormFreq;
		sinc.initBand();

		Latency = sinc.fl2;
		BlockSizeBits = getBitOccupancy( sinc.KernelLen - 1 );
		const int BlockSize = 1 << BlockSizeBits;
		KernelBlock.alloc( BlockSize * 2 );

		CDSPRealFFTKeeper ffto( BlockSizeBits + 1 );

		sinc.generateBandPow( &KernelBlock[ 0 ], pwr );
normalizeFIRFilter( &KernelBlock[ 0 ], sinc.KernelLen, 1.0 );
//		normalizeFIRFilter( &KernelBlock[ 0 ], sinc.KernelLen,
//			ffto -> getInvMulConst() );

		memset( &KernelBlock[ sinc.KernelLen ], 0,
			( BlockSize * 2 - sinc.KernelLen ) * sizeof( double ));

//		ffto -> forward( KernelBlock );
	}
};

/**
 * \brief FIR filters cache class.
 *
 * Class that implements cache for calculated FIR filters. The required FIR
 * filter should be obtained via the getLPFilter() static function.
 */

class CDSPFIRFilterCache : public R8B_BASECLASS
{
public:
	/**
	 * @param ReqTransBand Required transition band, in percent of the full
	 * bandwidth, in the range CDSPFIRFilter.getLPMinTransBand() to
	 * CDSPFIRFilter.getLPMaxTransBand(), inclusive. The transition band
	 * specifies part of the spectrum between the -3 dB and Nyquist points.
	 * The resulting -3 dB point varies in the range from -2.98 and -3.08 dB,
	 * but is generally very close to -3 dB.
	 * @param ReqAtten Required stop-band attenuation in decibel, in the range
	 * CDSPFIRFilter.getLPMinAtten() to CDSPFIRFilter.getLPMaxAtten(),
	 * inclusive. Note that the stop-band attenuation of the resulting filter
	 * may vary in the range +/- 1 decibel in comparison to the required
	 * value.
	 * @param ReqNormFreq Required normalized frequency, in the range 0 to 1,
	 * inclusive.
	 * @return A reference to a new or a previously calculated low-pass FIR
	 * filter object with the required characteristics. A reference count is
	 * incremented in the returned filter object which should be released
	 * after use.
	 */

	static CDSPFIRFilter& getLPFilter( const double ReqTransBand,
		const double ReqAtten, const double ReqNormFreq )
	{
	}
};

} // namespace r8b

#endif // R8B_CDSPFIRFILTER_INCLUDED
