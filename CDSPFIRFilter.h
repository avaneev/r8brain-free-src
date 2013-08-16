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
	R8BNOCTOR( CDSPFIRFilter );

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

	/**
	 * This function should be called when the filter obtained via the
	 * cache is no longer needed.
	 */

	void unref();

//private:
	double ReqTransBand; ///< Required transition band in percent, as passed
		///< by the user.
		///<
	double ReqAtten; ///< Required stop-band attenuation in decibel, as passed
		///< by the user (positive value).
		///<
	double ReqNormFreq; ///< Required normalized frequency, 0 to 1 inclusive.
	CDSPFIRFilter* Next; ///< Next FIR filter in cache's list.
	int RefCount; ///< The number of references made to *this FIR filter.
	int Latency; ///< Filter's latency in samples.
	int BlockSizeBits; ///< Block size used to store filter response,
		///< expressed as Nth power of 2. This value is used directly by the
		///< convolver.
		///<
	CFixedBuffer< double > KernelBlock; ///< Filter response buffer, size
		///< equals to 1 << ( BlockSizeBits + 1 ). Second part of the buffer
		///< contains zero-padding to allow alias-free convolution.
		///<

	CDSPFIRFilter()
		: RefCount( 1 )
	{
	}

	~CDSPFIRFilter()
	{
		R8BASSERT( RefCount == 0 );

		delete Next;
	}

	/**
	 * Function builds filter kernel based on the ReqTransBand, ReqAtten and
	 * ReqNormFreq parameters.
	 */

	void buildLPFilter()
	{
		const double tb = ReqTransBand * 0.01;
		double atten = -ReqAtten;

		// Apply 3-part correction to the "attenuation" for higher precision.

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
			atten -= 0.538238440603965 * sin( -0.1233515698667 * atten ) *
				atan2( 0.0143246704986262 * atten, cos( 0.0935923685623112 *
				atten )) + 0.696034798186251 * sin( 6.00222816844784 +
				0.0790029450006683 * atten ) * atan2( 0.257260396929507,
				8.32685945641861 + 0.0790029450006683 * atten ) -
				0.30920824764881;
		}

		// A set of "magic formulas" that calculate "half-band filter kernel
		// length" (hl), windowing function's power factor (pwr) and
		// "-3 dB frequency offset" (fo).

		const double hl = ( -10.3864830460583 - 0.28561497095004 * atten ) /
			tb + -215.422486627309 / ( tb * atten - 12.2962925122936 * tanh(
			tanh( tb )) - 0.340884649462501 * tb * atten * atan( sin(
			4.22129490874673 - 0.0953516921538414 * atten )));

		const double hltb = hl * tb;
		const double pwr = 0.185341477676462 + 0.116234005923982 * hltb +
			0.0183525684531612 * sin( 0.132513576412821 * hltb + gauss(
			0.132513576412821 * hltb )) - 0.768004731482082 * gauss( gauss(
			sin( 6.10653526213251 + 0.132513576412821 * hltb ) -
			0.137241895475179 * hltb ) - 0.329587582364288 );

		const double fo = atan2( 1.36778230157297 * asinh( 0.680688824721814 +
			pwr ), hl ) + gauss( cos( 2.37736468090842 * pwr ) /
			( 0.527409384163176 * hl + 6.42102056930156 * exp( sqr( pwr )))) -
			0.976478203607702 * tb;

		CDSPSincFilterGen sinc;
		sinc.Len2 = 0.25 * hl / ReqNormFreq;
		sinc.Freq1 = 0.0;
		sinc.Freq2 = M_PI * fo * ReqNormFreq;
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
	friend class CDSPFIRFilter;

public:
	/**
	 * @return The number of filters present in the cache now. This value can
	 * be monitored for debugging "forgotten" filters.
	 */

	static int getFilterCount()
	{
		R8BSYNC( StateSync );

		return( FilterCount );
	}

	/**
	 * @param ReqTransBand Required transition band, in percent of the full
	 * bandwidth, in the range CDSPFIRFilter::getLPMinTransBand() to
	 * CDSPFIRFilter::getLPMaxTransBand(), inclusive. The transition band
	 * specifies part of the spectrum between the -3 dB and Nyquist points.
	 * The resulting -3 dB point varies in the range from -2.98 and -3.08 dB,
	 * but is generally very close to -3 dB.
	 * @param ReqAtten Required stop-band attenuation in decibel, in the range
	 * CDSPFIRFilter::getLPMinAtten() to CDSPFIRFilter::getLPMaxAtten(),
	 * inclusive. Note that the stop-band attenuation of the resulting filter
	 * may vary in the range +/- 1 decibel in comparison to the required
	 * value.
	 * @param ReqNormFreq Required normalized frequency, in the range 0 to 1,
	 * inclusive.
	 * @return A reference to a new or a previously calculated low-pass FIR
	 * filter object with the required characteristics. A reference count is
	 * incremented in the returned filter object which should be released
	 * after use via the CDSPFIRFilter::unref() function.
	 */

	static CDSPFIRFilter& getLPFilter( const double ReqTransBand,
		const double ReqAtten, const double ReqNormFreq )
	{
		R8BSYNC( StateSync );

		CDSPFIRFilter* PrevFilter = NULL;
		CDSPFIRFilter* CurFilter = Filters;

		while( CurFilter != NULL )
		{
			if( CurFilter -> ReqTransBand == ReqTransBand &&
				CurFilter -> ReqAtten == ReqAtten &&
				CurFilter -> ReqNormFreq == ReqNormFreq )
			{
				break;
			}

			if( CurFilter -> Next == NULL &&
				FilterCount >= R8B_FILTER_CACHE_MAX )
			{
				if( CurFilter -> RefCount == 0 )
				{
					// Delete the last filter which is not used.

					PrevFilter -> Next = NULL;
					delete CurFilter;
					FilterCount--;
				}
				else
				{
					// Move the last filter to the top of the list since it
					// seems to be in use for a long time.

					PrevFilter -> Next = NULL;
					CurFilter -> Next = Filters.unkeep();
					Filters = CurFilter;
				}

				CurFilter = NULL;
				break;
			}

			PrevFilter = CurFilter;
			CurFilter = CurFilter -> Next;
		}

		if( CurFilter != NULL )
		{
			CurFilter -> RefCount++;

			if( PrevFilter == NULL )
			{
				return( *CurFilter );
			}

			// Remove the filter from the list temporarily.

			PrevFilter -> Next = CurFilter -> Next;
		}
		else
		{
			// Create a new filter object (with RefCount == 1) and build
			// filter kernel.

			CurFilter = new CDSPFIRFilter();
			CurFilter -> ReqTransBand = ReqTransBand;
			CurFilter -> ReqAtten = ReqAtten;
			CurFilter -> ReqNormFreq = ReqNormFreq;
			CurFilter -> buildLPFilter();
			FilterCount++;
		}

		// Insert the filter at the start of the list.

		CurFilter -> Next = Filters.unkeep();
		Filters = CurFilter;

		return( *CurFilter );
	}

private:
	static CSyncObject StateSync; ///< Cache state synchronizer.
	static CPtrKeeper< CDSPFIRFilter* > Filters; ///< The chain of cached
		///< filters.
		///<
	static int FilterCount; ///< The number of filters currently preset in the
		///< cache.
		///<
};

// ---------------------------------------------------------------------------
// CDSPFIRFilter PUBLIC
// ---------------------------------------------------------------------------

inline void CDSPFIRFilter :: unref()
{
	R8BSYNC( CDSPFIRFilterCache :: StateSync );

	RefCount--;
}

// ---------------------------------------------------------------------------

} // namespace r8b

#endif // R8B_CDSPFIRFILTER_INCLUDED
