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
	~CDSPFIRFilter()
	{
		R8BASSERT( RefCount == 0 );

		delete Next;
	}

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

private:
	double ReqNormFreq; ///< Required normalized frequency, 0 to 1 inclusive.
		///<
	double ReqTransBand; ///< Required transition band in percent, as passed
		///< by the user.
		///<
	double ReqAtten; ///< Required stop-band attenuation in decibel, as passed
		///< by the user (positive value).
		///<
	CDSPFIRFilter* Next; ///< Next FIR filter in cache's list.
		///<
	int RefCount; ///< The number of references made to *this FIR filter.
		///<
	int Latency; ///< Filter's latency in samples.
		///<
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

	/**
	 * Function builds filter kernel based on the ReqTransBand, ReqAtten and
	 * ReqNormFreq parameters.
	 */

	void buildLPFilter()
	{
		const double tb = ReqTransBand * 0.01;
		double atten = -ReqAtten;

		// Apply 4-range correction to the "attenuation" for +/- 1 dB
		// precision, without this correction the error is +/- 3 dB.

		if( atten >= -74.491228 )
		{
			double atten_b = 0.0;

			if( atten >= -50.651267 )
			{
				atten_b = 0.699900982886959 * cos( 0.222669162355806 *
					atten );
			}

			atten -= 0.438075758491934 + 1.72682101344218 * cos(
				0.230227912668262 * atten ) + 0.0102506237901691 * atten *
				sin( 0.334121655800277 * atten );

			atten -= atten_b;
		}
		else
		if( atten >= -100.491228 )
		{
			atten -= 2.20298948096653 * cos( 0.400387943600211 * atten ) +
				0.0162671715330083 * atten * cos( 0.400358395063225 *
				atten ) - 1.95268707933334 - 0.0263263087274513 * atten;
		}
		else
		{
			atten -= 12.4238553852216 * sin( 0.08178837006 * atten ) * atan2(
				2.477633394, 0.08178837006 * atten ) + 0.5689731045 * sin(
				-0.1238783946 * atten ) * atan2( 0.01319093007 * atten + sin(
				-0.1238783946 * atten ), cos( 0.6821178425 + 0.0962288404 *
				atten )) - 0.1170035732 - 35.00537255 * sin( 0.08178837006 *
				atten );
		}

		// A set of "magic formulas" that calculate "half-band filter kernel
		// length" (hl), windowing function's power factor (pwr) and
		// "-3 dB frequency offset" (fo).

		const double hl = ( -10.3864830460583 - 0.28561497095004 * atten ) /
			tb + -215.422486627309 / ( tb * atten - 12.2962925122936 * tanh(
			tanh( tb )) - 0.340884649462501 * tb * atten * atan( sin(
			4.22129490874673 - 0.0953516921538414 * atten )));

		const double hltb = hl * tb;
		const double pwr = 0.201981441582014 + 0.115979123679821 * hltb +
			0.0200435272760573 * sin( 0.234743434639915 + 0.130876326345944 *
			hltb ) - 0.780453063766734 * gauss( gauss( sin( 6.10822300874378 +
			0.130876326345944 * hltb ) - 0.136709763359556 * hltb ) -
			0.339451943355516 );

		const double fo = atan2( 1.36789104485968 * asinh( 0.681921659358349 +
			pwr ), hl ) + gauss(( 0.681921659358349 * cos( 2.2957524045949 *
			pwr ) + pwr * cos( 2.2957524045949 * pwr )) / ( 0.590803811989741*
			hl + 4.50625485711854 * exp( pwr + sqr( pwr )))) -
			0.976490268819982 * tb;

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
		normalizeFIRFilter( &KernelBlock[ 0 ], sinc.KernelLen,
			ffto -> getInvMulConst() );

		memset( &KernelBlock[ sinc.KernelLen ], 0,
			( BlockSize * 2 - sinc.KernelLen ) * sizeof( double ));

		ffto -> forward( KernelBlock );
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
	 * Function calculates or returns reference to a previously calculated
	 * (cached) low-pass FIR filter. Note that the real transition band and
	 * attenuation achieved by the filter varies with the magnitude of the
	 * required attenuation: the real filter specs can be considered precise,
	 * "to spec" if the required attenuation is above 95 dB.
	 *
	 * @param ReqNormFreq Required normalized frequency, in the range 0 to 1,
	 * inclusive. This is the point after which the stop-band spans.
	 * @param ReqTransBand Required transition band, in percent of the
	 * 0 to ReqNormFreq spectral bandwidth, in the range
	 * CDSPFIRFilter::getLPMinTransBand() to
	 * CDSPFIRFilter::getLPMaxTransBand(), inclusive. The transition band
	 * specifies part of the spectrum between the -3 dB and ReqNormFreq
	 * points. The real resulting -3 dB point varies in the range from -2.98
	 * to -3.09 dB, but is generally very close to -3 dB.
	 * @param ReqAtten Required stop-band attenuation in decibel, in the range
	 * CDSPFIRFilter::getLPMinAtten() to CDSPFIRFilter::getLPMaxAtten(),
	 * inclusive. Note that the stop-band attenuation of the resulting filter
	 * varies in the range +/- 1 decibel in comparison to the required value.
	 * @return A reference to a new or a previously calculated low-pass FIR
	 * filter object with the required characteristics. A reference count is
	 * incremented in the returned filter object which should be released
	 * after use via the CDSPFIRFilter::unref() function.
	 */

	static CDSPFIRFilter& getLPFilter( const double ReqNormFreq,
		const double ReqTransBand, const double ReqAtten )
	{
		R8BASSERT( ReqNormFreq >= 0.0 && ReqNormFreq <= 1.0 );
		R8BASSERT( ReqTransBand >= CDSPFIRFilter::getLPMinTransBand() );
		R8BASSERT( ReqTransBand <= CDSPFIRFilter::getLPMaxTransBand() );
		R8BASSERT( ReqAtten >= CDSPFIRFilter::getLPMinAtten() );
		R8BASSERT( ReqAtten <= CDSPFIRFilter::getLPMaxAtten() );

		R8BSYNC( StateSync );

		CDSPFIRFilter* PrevFilter = NULL;
		CDSPFIRFilter* CurFilter = Filters;

		while( CurFilter != NULL )
		{
			if( CurFilter -> ReqNormFreq == ReqNormFreq &&
				CurFilter -> ReqTransBand == ReqTransBand &&
				CurFilter -> ReqAtten == ReqAtten )
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
			CurFilter -> ReqNormFreq = ReqNormFreq;
			CurFilter -> ReqTransBand = ReqTransBand;
			CurFilter -> ReqAtten = ReqAtten;
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
		///<
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
