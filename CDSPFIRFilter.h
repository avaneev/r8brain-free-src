//$ nobt
//$ nocpp

/**
 * @file CDSPFIRFilter.h
 *
 * @brief FIR filter generator and filter cache classes.
 *
 * This file includes low-pass FIR filter generator and filter cache.
 *
 * r8brain-free-src Copyright (c) 2013-2026 Aleksey Vaneev
 *
 * See the "LICENSE" file for license.
 */

#ifndef R8B_CDSPFIRFILTER_INCLUDED
#define R8B_CDSPFIRFILTER_INCLUDED

#include "CDSPSincFilterGen.h"
#include "CDSPRealFFT.h"

namespace r8b {

/**
 * @brief Enumeration of filter's phase responses.
 */

enum EDSPFilterPhaseResponse
{
	fprLinearPhase = 0, ///< Linear-phase response. Features a linear-phase,
		///< high-latency response, with the latency expressed as an integer
		///< value.
	fprMinPhase ///< Minimum-phase response. Features a minimal-latency
		///< response, but the response's phase is non-linear. The latency is
		///< usually expressed as a non-integer value, and is usually small,
		///< but is never equal to zero. The minimum-phase filter is obtained
		///< from a linear-phase filter. Note that since in the context of
		///< r8brain-free-src other filters (interpolation, half-band) remain
		///< linear-phase, the resulting phase will be "intermediate". The
		///< minimum-phase transformation has precision limits: this may skew
		///< both the -3 dB point and attenuation of the filter being
		///< transformed: as it was measured, the skew happens purely at
		///< random, and in most cases is within tolerable range. In a small
		///< (1%) random subset of cases the skew is bigger and cannot be
		///< predicted. Minimum-phase transform requires 64-bit floating-point
		///< FFT; results with 32-bit float FFT are far from optimal.
};

/**
 * @brief Calculation and storage class for FIR filters.
 *
 * Class that implements calculation and storing of a FIR filter (currently
 * contains low-pass filter calculation routine designed for sample rate
 * conversion). Objects of this class cannot be created directly, but can be
 * obtained via the CDSPFilterCache::getLPFilter() static function.
 */

class CDSPFIRFilter : public R8B_BASECLASS,
	private CSinglyLinkedListItem< CDSPFIRFilter >
{
	R8BNOCTOR( CDSPFIRFilter )

	friend class CSinglyLinkedListItem< CDSPFIRFilter >;
	friend class CDSPFIRFilterCache;

public:
	/**
	 * @brief Returns the minimal allowed low-pass filter's transition band,
	 * in percent.
	 */

	static double getLPMinTransBand()
	{
		return( 0.5 );
	}

	/**
	 * @brief Returns the maximal allowed low-pass filter's transition band,
	 * in percent.
	 */

	static double getLPMaxTransBand()
	{
		return( 45.0 );
	}

	/**
	 * @brief Returns the minimal allowed low-pass filter's stop-band
	 * attenuation, in decibel.
	 */

	static double getLPMinAtten()
	{
		return( 49.0 );
	}

	/**
	 * @brief Returns the maximal allowed low-pass filter's stop-band
	 * attenuation, in decibel.
	 */

	static double getLPMaxAtten()
	{
		return( 218.0 );
	}

	/**
	 * @brief Returns `true` if kernel block of *this* filter has zero-phase
	 * response.
	 */

	bool isZeroPhase() const
	{
		return( IsZeroPhase );
	}

	/**
	 * @brief Returns filter's latency, in samples (integer part).
	 */

	int getLatency() const
	{
		return( Latency );
	}

	/**
	 * @brief Returns filter's latency, in samples (fractional part). Always
	 * zero for linear-phase filters.
	 */

	double getLatencyFrac() const
	{
		return( LatencyFrac );
	}

	/**
	 * @brief Returns filter kernel length, in samples. Not to be confused
	 * with the block length.
	 */

	int getKernelLen() const
	{
		return( KernelLen );
	}

	/**
	 * @brief Returns filter's block length, expressed as Nth power of 2.
	 *
	 * The actual length is twice as large due to zero-padding. Also it does
	 * not account for length reduction due to `ReqDownShift`.
	 */

	int getBlockLenBits() const
	{
		return( BlockLenBits );
	}

	/**
	 * @brief Returns pointer to filter's kernel block, in complex-numbered
	 * form obtained via the CDSPRealFFT::forward() function call.
	 *
	 * Zero-padded, gain-adjusted with the CDSPRealFFT::getInvMulConst()
	 * multiplied by `ReqGain` constant, immediately suitable for convolution.
	 * Kernel block may have "zero-phase" response, depending on the
	 * isZeroPhase() function's result.
	 *
	 * If `ReqDownShift` is positive, the actual length of the kernel block is
	 * smaller by `2^ReqDownShift`.
	 */

	const realfft_t* getKernelBlock() const
	{
		return( KernelBlockFFT );
	}

	/**
	 * @brief Calculates Nyquist bin value for downsampling using power-of-2
	 * down-shifting.
	 *
	 * Multiplies the forward transform's complex component with filter's
	 * response at the edge bin.
	 *
	 * @param re Input value, real part.
	 * @param im Input value, imaginary part.
	 * @return The value of the Nyquist bin.
	 */

	realfft_t calcDownShiftNyquist( const realfft_t re,
		const realfft_t im ) const
	{
		return( dsnRe * re - dsnIm * im );
	}

	/**
	 * @brief Reduces reference count to *this* object.
	 *
	 * This function should be called when the filter obtained via the filter
	 * cache is no longer needed.
	 */

	void unref();

private:
	double ReqNormFreq; ///< Required normalized frequency, 0 to 1 inclusive.
	double ReqTransBand; ///< Required transition band in percent, as passed
		///< by the user.
	double ReqAtten; ///< Required stop-band attenuation in decibel, as passed
		///< by the user (positive value).
	EDSPFilterPhaseResponse ReqPhase; ///< Required filter's phase response.
	double ReqGain; ///< Required overall filter's gain.
	int ReqDownShift; ///< Required downsampling shift, if a positive value.
	int RefCount; ///< The number of references made to *this* FIR filter.
	bool IsZeroPhase; ///< `true` if kernel block of *this* filter has
		///< zero-phase response.
	int Latency; ///< Filter's latency in samples (integer part).
	double LatencyFrac; ///< Filter's latency in samples (fractional part).
	int KernelLen; ///< Filter kernel length, in samples.
	int BlockLenBits; ///< Block length used to store *this* FIR filter,
		///< expressed as Nth power of 2. This value is used directly by the
		///< convolver.
	CFixedBuffer< double > KernelBlock; ///< FIR filter buffer, capacity
		///< equals to `1 << ( BlockLenBits + 1 )`. Second half of the buffer
		///< contains zero-padding to allow alias-free convolution.
		///< Address-aligned. Reused to store the FFT of the filter.
	realfft_t* KernelBlockFFT; ///< Pointer to FFT-transformed kernel block in
		///< native resolution.
	realfft_t dsnRe; ///< Down-shift Nyquist bin's real value.
	realfft_t dsnIm; ///< Down-shift Nyquist bin's imaginary value.

	CDSPFIRFilter()
		: RefCount( 1 )
	{
	}

	/**
	 * @brief Computes fixed-degree polynomial via Horner evaluation
	 * (ascending powers).
	 *
	 * @param c Polynomial coefficient.
	 * @param n Number of coefficient.
	 * @param x Argument [-1; 1].
	 */

	static double hpoly( const double* const c, const int n, const double x )
	{
		double r = c[ n - 1 ];
		int i;

		for( i = n - 2; i >= 0; i-- )
		{
			r = r * x + c[ i ];
		}

		return( r );
	}

	/**
	 * @brief Computes `pwr`, `hl`, `fo1` for CDSPFIRFilter::buildLPFilter().
	 *
	 * ## The formulas
	 *
	 * a = (A - 133.5)/84.5                                   in [-1, 1]
	 * b = P9(a) + 0.0523 + 4.07*T^2                          in [4.49, 24.72]
	 * u = (ln b - 2.3555)/0.8585                             in [-1, 1]
	 * s = (T/0.45)^2                                         in (0, 1]
	 * delta = exp(D0(u) + s^2*D2(u) + s^3*D3(u) + s^4*D4(u)) in [0.16, 0.72]
	 * Y = exp(Y0(u) + s^2*Y2(u) + s^3*Y3(u))                 in [1.59, 3.48]
	 * X = (4/pi)*(sqrt(b^2 + pi^2) - delta)                  in [6.19, 31.5]
	 *
	 * hl  = (X + Y)/T          fo1 = X*T/(X + Y)          pwr = b/125
	 *
	 * `P9` and the `D/Y` pieces are fixed-degree polynomials (Horner), all
	 * coefficients in the header. Beta stays 125.0; `pwr` lands in
	 * [0.036, 0.198] so the window's `pow_a` sees an effective beta
	 * `125*pwr = b`.
	 *
	 * ## Why this shape
	 *
	 * The three outputs decouple into three 1-D physical roles:
	 *
	 * - `X = fo1*hl` positions the transition-band main-lobe descent so it
	 *   lands at Nyquist at exactly -A dB. The landing sits just *before* the
	 *   Kaiser window's first null at OmegaL = sqrt(b^2+pi^2); in lobe units
	 *   X = (4/pi)(sqrt(b^2+pi^2) - delta), where delta in (0, 1] is the gap
	 *   to the null. The `delta(A,T)` polynomial surface is the fit of the
	 *   exact solve.
	 *
	 * - `Y = (T-fo1)*hl` sets the -3 dB point at the transition start;
	 *   its dependence is gentle (~=3 dB/ln-unit), so a small surface
	 *   suffices.
	 *
	 * - `b = 125*pwr` sets every sidelobe level via the measured
	 *   first-sidelobe law of the power-raised Kaiser-125 family;
	 *   the `0.0523 + 4.07*T^2` bump is the smooth margin that keeps all
	 *   sidelobe peaks below the landing at wide T, where tap-count
	 *   quantization (`fl2 = floor(L)`, L as small as ~9) makes the exact
	 *   locus ripple by +/-4 dB -- a non-smooth effect no closed form can
	 *   track, so it's covered by margin instead.
	 *
	 * ## Singularity-freeness
	 *
	 * Every operation is provably safe on the closed domain T in
	 * [0.005, 0.45], A in [49, 218]: `log` of b >= 4.49 > 0; `exp` of bounded
	 * polynomials (finite values, strictly positive results);
	 * `sqrt(b^2+pi^2)` >= 5.48 > delta >= 0.16 so X >= 6.19 > 0; divisions
	 * only by T >= 0.005 and X+Y >= 7.78. All polynomial arguments are
	 * affinely normalized to [-1, 1]. Consequently hl in [21.2, 6581]
	 * (so `Len2 = 0.25*hl/Nf >= 5.3 >= 2` for any Nf <= 1), fo1 in (0, T)
	 * subset (0, 0.45), pwr in [0.036, 0.198]. No table lookups, no branches,
	 * no extrapolation regimes.
	 *
	 * @param T Required transition band, normalized [0.005, 0.45].
	 * @param ReqAtten Required attenuation in decibel [49.0, 218.0].
	 */

	static void calcDesignParams( const double T, const double ReqAtten,
		double& pwr, double& hl, double& fo1 )
	{
		// b0(a): effective Kaiser beta vs attenuation, deg-9 in
		// a = (A-133.5)/84.5
		static const double LPF_CB[ 10 ] = {
			14.409935832973984, 9.6626843693138884, -0.17435120318333536,
			0.47264599329014079, 0.4940426137302133, -1.5662364718692721,
			-0.71006021971136335, 2.0294051576357957, 0.22363464796569277,
			-0.78254591651630068 };

		// log(delta) = D0(u) + s^2*D2(u) + s^3*D3(u) + s^4*D4(u),
		// s = (T/0.45)^2
		static const double LPF_D0[ 7 ] = {
			-1.0753246533566314, -0.75662395688736628, -0.0047894065672392377,
			0.08577566962378394, -0.17001888073403126, -0.013628259555023822,
			0.093950306915379783 };
		static const double LPF_D2[ 3 ] = {
			1.1482209047708865, -0.052152189285352819, -0.020187508756994623 };
		static const double LPF_D3[ 3 ] = {
			-1.4201039645189557, 0.1242900889796041, -0.028854770456694875 };
		static const double LPF_D4[ 3 ] = {
			0.55500996383110457, -0.062635292370927881, 0.023661985245151548 };

		// log(Y) = Y0(u) + s^2*Y2(u) + s^3*Y3(u)
		static const double LPF_Y0[ 5 ] = {
			0.8399122624879678, 0.39839640583388353, 0.015028625919905413,
			-0.005569192458538463, 0.001494097359614131 };
		static const double LPF_Y2[ 3 ] = {
			0.00069542490312080887, -0.00062890127106954123,
			-0.00071460219184551251 };
		static const double LPF_Y3[ 2 ] = {
			-0.00065713366734133005, 0.0011106064299533634 };

		const double a = ( ReqAtten - 133.5 ) / 84.5;
		const double b = hpoly( LPF_CB, 10, a ) +
			0.0523 + 4.07 * T * T;

		const double u = ( log( b ) - 2.3555 ) / 0.8585;

		const double s = ( T / 0.45 ) * ( T / 0.45 );
		const double s2 = s * s;
		const double s3 = s2 * s;
		const double s4 = s2 * s2;

		const double D = exp( hpoly( LPF_D0, 7, u ) +
			s2 * hpoly( LPF_D2, 3, u ) +
			s3 * hpoly( LPF_D3, 3, u ) +
			s4 * hpoly( LPF_D4, 3, u ));

		const double Y = exp( hpoly( LPF_Y0, 5, u ) +
			s2 * hpoly( LPF_Y2, 3, u ) +
			s3 * hpoly( LPF_Y3, 2, u ));

		const double X = ( 4.0 / R8B_PI ) *
			( sqrt( b * b + R8B_PI * R8B_PI ) - D );

		hl = ( X + Y ) / T;
		fo1 = X * T / ( X + Y );
		pwr = b / 125.0;
	}

	/**
	 * @brief Builds filter kernel based on the "Req" parameters.
	 */

	void buildLPFilter()
	{
		double pwr, fo1, hl;

		calcDesignParams( ReqTransBand * 0.01, ReqAtten, pwr, hl, fo1 );

		double WinParams[ 2 ];
		WinParams[ 0 ] = 125.0;
		WinParams[ 1 ] = pwr;

		CDSPSincFilterGen sinc;
		sinc.Len2 = 0.25 * hl / ReqNormFreq;
		sinc.Freq1 = 0.0;
		sinc.Freq2 = R8B_PI * ( 1.0 - fo1 ) * ReqNormFreq;
		sinc.initBand( CDSPSincFilterGen :: wftKaiser, WinParams, true );

		KernelLen = sinc.KernelLen;
		BlockLenBits = getBitOccupancy( KernelLen - 1 ) + R8B_EXTFFT;
		const int BlockLen = 1 << BlockLenBits;

		KernelBlock.alloc( BlockLen * 2 );
		sinc.generateBand( &KernelBlock[ 0 ],
			&CDSPSincFilterGen :: calcWindowKaiser );

		if( ReqPhase == fprLinearPhase )
		{
			IsZeroPhase = true;
			Latency = sinc.fl2;
			LatencyFrac = 0.0;
		}
		else
		{
			IsZeroPhase = false;
			double DCGroupDelay;

			calcMinPhaseTransform( &KernelBlock[ 0 ], KernelLen, 16, false,
				&DCGroupDelay );

			Latency = (int) DCGroupDelay;
			LatencyFrac = DCGroupDelay - Latency;
		}

		CDSPRealFFTKeeper ffto( BlockLenBits + 1 );

		dsnRe = 0;
		dsnIm = 0;

		if( IsZeroPhase )
		{
			// Calculate DC gain.

			double s = 0.0;
			int i;

			for( i = 0; i < KernelLen; i++ )
			{
				s += KernelBlock[ i ];
			}

			s = ffto -> getInvMulConst() * ReqGain / s;

			// Time-shift the filter so that zero-phase response is produced.
			// Simultaneously multiply by "s".

			for( i = 0; i <= sinc.fl2; i++ )
			{
				KernelBlock[ i ] = KernelBlock[ sinc.fl2 + i ] * s;
			}

			for( i = 1; i <= sinc.fl2; i++ )
			{
				KernelBlock[ BlockLen * 2 - i ] = KernelBlock[ i ];
			}

			memset( &KernelBlock[ sinc.fl2 + 1 ], 0,
				(size_t) ( BlockLen * 2 - KernelLen ) *
				sizeof( KernelBlock[ 0 ]));

			KernelBlockFFT = ffto -> forward( KernelBlock );

			if( ReqDownShift > 0 )
			{
				KernelBlockFFT = ffto -> reorderForward( KernelBlockFFT,
					ffto -> getWorkBuf() );

				const int z = ( 2 << BlockLenBits ) >> ReqDownShift;
				dsnRe = KernelBlockFFT[ z ];
				dsnIm = KernelBlockFFT[ z + 1 ];

				CDSPRealFFTKeeper ffto2( BlockLenBits + 1 - ReqDownShift );
				CFixedBuffer< double > NewBlock( ffto2 -> getLen() *
					(int) sizeof( realfft_t ) / (int) sizeof( double ));

				KernelBlockFFT = ffto2 -> reorderInverse( KernelBlockFFT,
					NewBlock );

				ffto2 -> convertToZP( KernelBlockFFT );

				KernelBlock.moveFrom( NewBlock );
			}
			else
			{
				ffto -> convertToZP( KernelBlockFFT );
			}
		}
		else
		{
			normalizeFIRFilter( &KernelBlock[ 0 ], KernelLen,
				ffto -> getInvMulConst() * ReqGain );

			memset( &KernelBlock[ KernelLen ], 0,
				(size_t) ( BlockLen * 2 - KernelLen ) *
				sizeof( KernelBlock[ 0 ]));

			KernelBlockFFT = ffto -> forward( KernelBlock );

			if( ReqDownShift > 0 )
			{
				KernelBlockFFT = ffto -> reorderForward( KernelBlockFFT,
					ffto -> getWorkBuf() );

				const int z = ( 2 << BlockLenBits ) >> ReqDownShift;
				dsnRe = KernelBlockFFT[ z ];
				dsnIm = KernelBlockFFT[ z + 1 ];

				CDSPRealFFTKeeper ffto2( BlockLenBits + 1 - ReqDownShift );
				CFixedBuffer< double > NewBlock( ffto2 -> getLen() *
					(int) sizeof( realfft_t ) / (int) sizeof( double ));

				KernelBlockFFT = ffto2 -> reorderInverse( KernelBlockFFT,
					NewBlock );

				KernelBlock.moveFrom( NewBlock );
			}
		}

		R8BCONSOLE( "CDSPFIRFilter: flt_len=%i latency=%i nfreq=%.4f "
			"tb=%.1f att=%.1f gain=%.3f dsh=%i\n", KernelLen, Latency,
			ReqNormFreq, ReqTransBand, ReqAtten, ReqGain, ReqDownShift );
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
	R8BNOCTOR( CDSPFIRFilterCache )

	friend class CDSPFIRFilter;

public:
	/**
	 * @brief Returns the number of filters present in the cache now. This
	 * value can be monitored for debugging "forgotten" filters.
	 */

	static int getObjCount()
	{
		R8BSYNC( getStateSync() );

		return( getObjCountStatic() );
	}

	/**
	 * @brief Calculates or returns reference to a previously calculated
	 * (cached) low-pass FIR filter.
	 *
	 * Note that the real transition band and attenuation achieved by the
	 * filter varies with the magnitude of the required attenuation, and are
	 * never 100% exact.
	 *
	 * @param ReqNormFreq Required normalized frequency, in the range 0 to 1,
	 * inclusive. This is the point after which the stop-band spans.
	 * @param ReqTransBand Required transition band, in percent of the
	 * 0 to `ReqNormFreq` spectral bandwidth, in the range
	 * CDSPFIRFilter::getLPMinTransBand() to
	 * CDSPFIRFilter::getLPMaxTransBand(), inclusive. The transition band
	 * specifies the part of the spectrum between the -3 dB and `ReqNormFreq`
	 * points. The real resulting -3 dB point varies in the range from -3.00
	 * to -3.05 dB, but is generally very close to -3 dB.
	 * @param ReqAtten Required stop-band attenuation in decibel, in the range
	 * CDSPFIRFilter::getLPMinAtten() to CDSPFIRFilter::getLPMaxAtten(),
	 * inclusive. Note that the actual stop-band attenuation of the resulting
	 * filter may be 0.40-4.46 dB higher.
	 * @param ReqPhase Required filter's phase response.
	 * @param ReqGain Required overall filter's gain (1.0 for unity gain).
	 * @param ReqDownFactor Required downsampling factor; it must be positive.
	 * @param AttenCorrs Attentuation correction table, to pass to the filter
	 * generation function. For internal use.
	 * @see EDSPFilterPhaseResponse
	 * @return A reference to a new or a previously calculated low-pass FIR
	 * filter object with the required characteristics. A reference count is
	 * incremented in the returned filter object which should be released
	 * after use via the CDSPFIRFilter::unref() function.
	 */

	static CDSPFIRFilter& getLPFilter( const double ReqNormFreq,
		const double ReqTransBand, const double ReqAtten,
		const EDSPFilterPhaseResponse ReqPhase, const double ReqGain,
		const int ReqDownFactor )
	{
		R8BASSERT( ReqNormFreq > 0.0 && ReqNormFreq <= 1.0 );
		R8BASSERT( ReqTransBand >= CDSPFIRFilter :: getLPMinTransBand() );
		R8BASSERT( ReqTransBand <= CDSPFIRFilter :: getLPMaxTransBand() );
		R8BASSERT( ReqAtten >= CDSPFIRFilter :: getLPMinAtten() );
		R8BASSERT( ReqAtten <= CDSPFIRFilter :: getLPMaxAtten() );
		R8BASSERT( ReqGain > 0.0 );
		R8BASSERT( ReqDownFactor > 0 );

		R8B_EXITDTOR static CPtrKeeper< CDSPFIRFilter > Objects; // The chain
			// of cached objects.

		int ReqDownShift = getBitOccupancy( ReqDownFactor ) - 1;

		if(( 1 << ReqDownShift ) != ReqDownFactor )
		{
			ReqDownShift = 0;
		}

		R8BSYNC( getStateSync() );

		int& ObjCount = getObjCountStatic();
		CDSPFIRFilter* PrevObj = R8B_NULL;
		CDSPFIRFilter* CurObj = Objects;

		while( CurObj != R8B_NULL )
		{
			if( CurObj -> ReqNormFreq == ReqNormFreq &&
				CurObj -> ReqTransBand == ReqTransBand &&
				CurObj -> ReqGain == ReqGain &&
				CurObj -> ReqDownShift == ReqDownShift &&
				CurObj -> ReqAtten == ReqAtten &&
				CurObj -> ReqPhase == ReqPhase )
			{
				break;
			}

			if( CurObj -> Next == R8B_NULL &&
				ObjCount >= R8B_FILTER_CACHE_MAX )
			{
				if( CurObj -> RefCount == 0 )
				{
					// Delete the last filter which is not used.

					PrevObj -> Next = R8B_NULL;
					delete CurObj;
					ObjCount--;
				}
				else
				{
					// Move the last filter to the top of the list since it
					// seems to be in use for a long time.

					PrevObj -> Next = R8B_NULL;
					CurObj -> Next = Objects.unkeep();
					Objects = CurObj;
				}

				CurObj = R8B_NULL;
				break;
			}

			PrevObj = CurObj;
			CurObj = CurObj -> Next;
		}

		if( CurObj != R8B_NULL )
		{
			CurObj -> RefCount++;

			if( PrevObj == R8B_NULL )
			{
				return( *CurObj );
			}

			// Remove the filter from the list temporarily.

			PrevObj -> Next = CurObj -> Next;
		}
		else
		{
			// Create a new filter object (with RefCount == 1) and build the
			// filter kernel.

			CPtrKeeper< CDSPFIRFilter > f( new CDSPFIRFilter() );

			f -> ReqNormFreq = ReqNormFreq;
			f -> ReqTransBand = ReqTransBand;
			f -> ReqAtten = ReqAtten;
			f -> ReqPhase = ReqPhase;
			f -> ReqGain = ReqGain;
			f -> ReqDownShift = ReqDownShift;

			f -> buildLPFilter();

			CurObj = f.unkeep();
			ObjCount++;
		}

		// Insert the filter at the start of the list.

		CurObj -> Next = Objects.unkeep();
		Objects = CurObj;

		return( *CurObj );
	}

protected:
	/**
	 * @brief Returns reference to filter cache sync object.
	 */

	static CSyncObject& getStateSync()
	{
		R8B_EXITDTOR static CSyncObject StateSync;

		return( StateSync );
	}

	/**
	 * @brief Returns reference to variable containing cache object count.
	 */

	static int& getObjCountStatic()
	{
		R8B_EXITDTOR static int ObjCount = 0; // The number of objects
			// currently preset in the cache.

		return( ObjCount );
	}
};

// ---------------------------------------------------------------------------
// CDSPFIRFilter PUBLIC
// ---------------------------------------------------------------------------

inline void CDSPFIRFilter :: unref()
{
	R8BSYNC( CDSPFIRFilterCache :: getStateSync() );

	RefCount--;
}

// ---------------------------------------------------------------------------

} // namespace r8b

#endif // R8B_CDSPFIRFILTER_INCLUDED
