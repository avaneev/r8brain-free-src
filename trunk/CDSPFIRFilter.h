//$ nocpp

/**
 * @file CDSPFIRFilter.h
 *
 * @brief FIR filter generator and filter cache classes.
 *
 * This file includes low-pass FIR filter generator and filter cache.
 *
 * r8brain-free-src Copyright (c) 2013 Aleksey Vaneev
 * See the "License.txt" file for license.
 */

#ifndef R8B_CDSPFIRFILTER_INCLUDED
#define R8B_CDSPFIRFILTER_INCLUDED

#include "CDSPFIRFilterParams.h"
#include "CDSPRealFFT.h"

namespace r8b {

/**
 * Enumeration of filter's phase responses.
 */

enum EDSPFilterPhaseResponse
{
	fprLinearPhase = 0, ///< Linear-phase response. Features a linear-phase
		///< high-latency response, with the latency expressed as integer
		///< value.
	fprMinPhase ///< Minimum-phase response. Features a minimal latency
		///< response, but the response's phase is non-linear. The latency is
		///< usually expressed as non-integer value, and usually is small, but
		///< is never equal to zero. The minimum-phase filter is transformed
		///< from the linear-phase filter. The transformation has precision
		///< limits which may skew both the -3 dB point and attenuation of the
		///< filter being transformed: as it was measured, the skew happens
		///< purely at random, and in most cases it is within tolerable range.
		///< In a small (1%) random subset of cases the skew is bigger and
		///< cannot be predicted.
};

/**
 * @brief Calculation and storage class for FIR filters.
 *
 * Class that implements calculation and storing of a FIR filter (currently
 * contains low-pass filter calculation routines designed for sample rate
 * conversion). Objects of this class cannot be created directly, but can be
 * obtained via the CDSPFilterCache::getLPFilter() static function.
 */

class CDSPFIRFilter : public R8B_BASECLASS
{
	R8BNOCTOR( CDSPFIRFilter );

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
		return( 45.0 );
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
	 * @return Filter's latency, in samples (integer part).
	 */

	int getLatency() const
	{
		return( Latency );
	}

	/**
	 * @return Filter's latency, in samples (fractional part). Always zero for
	 * linear-phase filters.
	 */

	double getLatencyFrac() const
	{
		return( LatencyFrac );
	}

	/**
	 * @return Filter kernel length, in samples. Not to be confused with the
	 * block length.
	 */

	int getKernelLen() const
	{
		return( KernelLen );
	}

	/**
	 * @return Filter's block length, espressed as Nth power of 2. The actual
	 * length is twice as large due to zero-padding.
	 */

	int getBlockLenBits() const
	{
		return( BlockLenBits );
	}

	/**
	 * @return Filter's kernel block, in complex-numbered form obtained via
	 * the CDSPRealFFT::forward() function call, zero-padded, gain-adjusted
	 * with the CDSPRealFFT::getInvMulConst() constant, immediately suitable
	 * for convolution.
	 */

	const double* getKernelBlock() const
	{
		return( KernelBlock );
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
	EDSPFilterPhaseResponse ReqPhase; ///< Required filter's phase response.
		///<
	CDSPFIRFilter* Next; ///< Next FIR filter in cache's list.
		///<
	int RefCount; ///< The number of references made to *this FIR filter.
		///<
	int Latency; ///< Filter's latency in samples (integer part).
		///<
	double LatencyFrac; ///< Filter's latency in samples (fractional part).
		///<
	int KernelLen; ///< Filter kernel length, in samples.
		///<
	int BlockLenBits; ///< Block length used to store *this FIR filter,
		///< expressed as Nth power of 2. This value is used directly by the
		///< convolver.
		///<
	CFixedBuffer< double > KernelBlock; ///< FIR filter buffer, capacity
		///< equals to 1 << ( BlockLenBits + 1 ). Second part of the buffer
		///< contains zero-padding to allow alias-free convolution.
		///<

	CDSPFIRFilter()
		: RefCount( 1 )
	{
	}

	/**
	 * Function builds filter kernel based on the supplied parameters.
	 *
	 * @param Params "Standardized" filter parameters.
	 */

	void buildLPFilter( const CDSPFIRFilterParams& Params )
	{
		CDSPSincFilterGen sinc;
		sinc.Len2 = 0.25 * Params.hl / ReqNormFreq;
		sinc.Freq1 = 0.0;
		sinc.Freq2 = M_PI * Params.fo * ReqNormFreq;
		sinc.initBand();

		KernelLen = sinc.KernelLen;
		BlockLenBits = getBitOccupancy( KernelLen - 1 );
		const int BlockLen = 1 << BlockLenBits;

		KernelBlock.alloc( BlockLen * 2 );
		sinc.generateBandPow( &KernelBlock[ 0 ], Params.pwr, Params.WinFunc );

		if( ReqPhase == fprLinearPhase )
		{
			Latency = sinc.fl2;
			LatencyFrac = 0.0;
		}
		else
		{
			double DCGroupDelay;

			calcMinPhaseTransform( &KernelBlock[ 0 ], KernelLen, 3, false,
				&DCGroupDelay );

			Latency = (int) DCGroupDelay;
			LatencyFrac = DCGroupDelay - Latency;
		}

		CDSPRealFFTKeeper ffto( BlockLenBits + 1 );

		normalizeFIRFilter( &KernelBlock[ 0 ], KernelLen,
			ffto -> getInvMulConst() );

		memset( &KernelBlock[ KernelLen ], 0,
			( BlockLen * 2 - KernelLen ) * sizeof( double ));

		ffto -> forward( KernelBlock );
	}
};

/**
 * @brief FIR filter cache class.
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
	 * "to spec" if the required attenuation is above 95 dB. Transition bands
	 * between 25% and 45% should be considered "coarse" as they have a higher
	 * error in both the -3 dB point position and attenuation.
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
	 * @param ReqPhase Required filter's phase response.
	 */

	static CDSPFIRFilter& getLPFilter( const double ReqNormFreq,
		const double ReqTransBand, const double ReqAtten,
		const EDSPFilterPhaseResponse ReqPhase )
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
				CurFilter -> ReqAtten == ReqAtten &&
				CurFilter -> ReqPhase == ReqPhase )
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
			// Create a new filter object (with RefCount == 1) and build the
			// filter kernel.

			CurFilter = new CDSPFIRFilter();
			CurFilter -> ReqNormFreq = ReqNormFreq;
			CurFilter -> ReqTransBand = ReqTransBand;
			CurFilter -> ReqAtten = ReqAtten;
			CurFilter -> ReqPhase = ReqPhase;
			FilterCount++;

			CDSPFIRFilterParams Params;
			Params.calcLPFBlackman( ReqTransBand, ReqAtten );
			CurFilter -> buildLPFilter( Params );
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
