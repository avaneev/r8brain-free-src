//$ nocpp

/**
 * @file CDSPFracInterpolator.h
 *
 * @brief Fractional delay interpolator and filter bank classes.
 *
 * This file includes fractional delay interpolator class.
 *
 * r8brain-free-src Copyright (c) 2013-2019 Aleksey Vaneev
 * See the "License.txt" file for license.
 */

#ifndef R8B_CDSPFRACINTERPOLATOR_INCLUDED
#define R8B_CDSPFRACINTERPOLATOR_INCLUDED

#include "CDSPSincFilterGen.h"
#include "CDSPProcessor.h"

namespace r8b {

/**
 * @brief Sinc function-based fractional delay filter bank class.
 *
 * Class implements storage and initialization of a bank of sinc-based
 * fractional delay filters, expressed as 0, 1st, 2nd or 3rd order polynomial
 * interpolation coefficients. The filters are windowed by either the "Vaneev"
 * or "Kaiser" power-raised window function. The FilterLen and FilterFracs
 * parameters can be varied freely without breaking the resampler.
 */

class CDSPFracDelayFilterBank : public R8B_BASECLASS
{
	R8BNOCTOR( CDSPFracDelayFilterBank );

	friend class CDSPFracDelayFilterBankCache;

public:
	/**
	 * Constructor.
	 *
	 * @param aFilterLen Specifies the number of samples (taps) each
	 * fractional delay filter should have. This must be an even value, the
	 * minimal value for FilterLen is 6, the maximal value is 30. To achieve a
	 * higher resampling precision, the oversampling should be used in the
	 * first place instead of using a higher FilterLen value. The lower this
	 * value is the lower the signal-to-noise performance of the interpolator
	 * will be. Each FilterLen decrease by 2 decreases SNR by approximately 12
	 * to 14 decibel.
	 * @param aFilterFracs The number of fractional delay positions to sample.
	 * For a high signal-to-noise ratio this has to be a larger value. The
	 * larger the FilterLen is the larger the FilterFracs should be.
	 * Approximate FilterLen to FilterFracs correspondence (for 2nd order
	 * interpolation only): 6:11, 8:17, 10:23, 12:41, 14:67, 16:97, 18:137,
	 * 20:211, 22:353, 24:673, 26:1051, 28:1733, 30:2833. The FilterFracs can
	 * be considerably reduced with 3rd order interpolation in use. In order
	 * to get consistent results when resampling to/from different sample
	 * rates, it is suggested to set this parameter to a suitable prime
	 * number.
	 * @param aElementSize The size of each filter's tap, in "double" values.
	 * This parameter corresponds to the complexity of interpolation. 4 should
	 * be set for 3rd order, 3 for 2nd order, 2 for linear interpolation, 1
	 * for whole-numbered stepping.
	 * @param aInterpPoints The number of points the interpolation is based
	 * on. This value should not be confused with the ElementSize. Set to 2
	 * for linear or no interpolation.
	 */

	CDSPFracDelayFilterBank( const int aFilterLen, const int aFilterFracs,
		const int aElementSize, const int aInterpPoints )
		: FilterLen( aFilterLen )
		, FilterFracs( aFilterFracs )
		, ElementSize( aElementSize )
		, InterpPoints( aInterpPoints )
		, FilterSize( FilterLen * ElementSize )
		, Next( NULL )
		, RefCount( 1 )
	{
		R8BASSERT( FilterLen >= 6 );
		R8BASSERT( FilterLen <= 30 );
		R8BASSERT(( FilterLen & 1 ) == 0 );
		R8BASSERT( FilterFracs > 0 );
		R8BASSERT( ElementSize >= 1 && ElementSize <= 4 );
		R8BASSERT( InterpPoints == 2 || InterpPoints == 8 );

		Table.alloc( FilterSize * ( FilterFracs + InterpPoints ));

		calculate();
	}

	~CDSPFracDelayFilterBank()
	{
		delete Next;
	}

	/**
	 * Function calculates the filter bank.
	 *
	 * @param Params Window function's parameters. If NULL then the built-in
	 * table values for the current FilterLen will be used.
	 */

	void calculate( const double* const Params = NULL )
	{
		R8BCONSOLE( "CDSPFracDelayFilterBank: flt_len=%i fracs=%i order=%i\n",
			FilterLen, FilterFracs, ElementSize - 1 );

		CDSPSincFilterGen sinc;
		sinc.Len2 = FilterLen / 2;

		double* p = Table;
		const int pc2 = InterpPoints / 2;
		int i;

		if( FilterLen <= 20 )
		{
			for( i = -pc2 + 1; i <= FilterFracs + pc2; i++ )
			{
				sinc.FracDelay = (double) ( FilterFracs - i ) / FilterFracs;
				sinc.initFrac( CDSPSincFilterGen :: wftVaneev, Params );
				sinc.generateFrac( p, &CDSPSincFilterGen :: calcWindowVaneev,
					ElementSize );

				normalizeFIRFilter( p, FilterLen, 1.0, ElementSize );
				p += FilterSize;
			}
		}
		else
		{
			for( i = -pc2 + 1; i <= FilterFracs + pc2; i++ )
			{
				sinc.FracDelay = (double) ( FilterFracs - i ) / FilterFracs;
				sinc.initFrac( CDSPSincFilterGen :: wftKaiser, Params, true );
				sinc.generateFrac( p, &CDSPSincFilterGen :: calcWindowKaiser,
					ElementSize );

				normalizeFIRFilter( p, FilterLen, 1.0, ElementSize );
				p += FilterSize;
			}
		}

		const int TablePos2 = FilterSize;
		const int TablePos3 = FilterSize * 2;
		const int TablePos4 = FilterSize * 3;
		const int TablePos5 = FilterSize * 4;
		const int TablePos6 = FilterSize * 5;
		const int TablePos7 = FilterSize * 6;
		const int TablePos8 = FilterSize * 7;
		double* const TableEnd = Table + ( FilterFracs + 1 ) * FilterSize;
		p = Table;

		if( InterpPoints == 8 )
		{
			if( ElementSize == 3 )
			{
				// Calculate 2nd order spline (polynomial) interpolation
				// coefficients using 8 points.

				while( p < TableEnd )
				{
					calcSpline2p8Coeffs( p, p[ 0 ], p[ TablePos2 ],
						p[ TablePos3 ], p[ TablePos4 ], p[ TablePos5 ],
						p[ TablePos6 ], p[ TablePos7 ], p[ TablePos8 ]);

					p += ElementSize;
				}
			}
			else
			if( ElementSize == 4 )
			{
				// Calculate 3rd order spline (polynomial) interpolation
				// coefficients using 8 points.

				while( p < TableEnd )
				{
					calcSpline3p8Coeffs( p, p[ 0 ], p[ TablePos2 ],
						p[ TablePos3 ], p[ TablePos4 ], p[ TablePos5 ],
						p[ TablePos6 ], p[ TablePos7 ], p[ TablePos8 ]);

					p += ElementSize;
				}
			}
		}
		else
		{
			if( ElementSize == 2 )
			{
				// Calculate linear interpolation coefficients.

				while( p < TableEnd )
				{
					p[ 1 ] = p[ TablePos2 ] - p[ 0 ];
					p += ElementSize;
				}
			}
		}
	}

	/**
	 * @param i Filter index, in the range 0 to FilterFracs, inclusive.
	 * @return Reference to the filter.
	 */

	const double& operator []( const int i ) const
	{
		R8BASSERT( i >= 0 && i <= FilterFracs );

		return( Table[ i * FilterSize ]);
	}

	/**
	 * This function should be called when the filter obtained via the
	 * filter bank cache is no longer needed.
	 */

	void unref();

private:
	int FilterLen; ///< Filter length.
		///<
	int FilterFracs; ///< Fractional position count.
		///<
	int ElementSize; ///< Filter element size.
		///<
	int InterpPoints; ///< Interpolation points to use.
		///<
	int FilterSize; ///< This constant specifies the "size" of a single filter
		///< in "double" elements.
		///<
	CFixedBuffer< double > Table; ///< The table of fractional delay filters
		///< for all discrete fractional x = 0..1 sample positions, and
		///< interpolation coefficients.
		///<
	CDSPFracDelayFilterBank* Next; ///< Next filter bank in cache's list.
		///<
	int RefCount; ///< The number of references made to *this filter bank.
		///<
};

/**
 * @brief Fractional delay filter cache class.
 *
 * Class implements cache storage of fractional delay filter banks.
 */

class CDSPFracDelayFilterBankCache : public R8B_BASECLASS
{
	R8BNOCTOR( CDSPFracDelayFilterBankCache );

	friend class CDSPFracDelayFilterBank;

public:
	/**
	 * @return The number of filters present in the cache now. This value can
	 * be monitored for debugging "forgotten" filters.
	 */

	static int getObjCount()
	{
		R8BSYNC( StateSync );

		return( ObjCount );
	}

	/**
	 * Function calculates or returns reference to a previously calculated
	 * (cached) fractional delay filter bank.
	 *
	 * @param aFilterLen Specifies the number of samples (taps) each
	 * fractional delay filter should have.
	 * @param aFilterFracs The number of fractional delay positions to sample.
	 * @param aElementSize The size of each filter's tap, in "double" values.
	 * @param aInterpPoints The number of points the interpolation is based
	 * on.
	 */

	static CDSPFracDelayFilterBank& getFilterBank( const int aFilterLen,
		const int aFilterFracs, const int aElementSize,
		const int aInterpPoints )
	{
		R8BSYNC( StateSync );

		CDSPFracDelayFilterBank* PrevObj = NULL;
		CDSPFracDelayFilterBank* CurObj = Objects;

		while( CurObj != NULL )
		{
			if( CurObj -> FilterLen == aFilterLen &&
				CurObj -> FilterFracs == aFilterFracs &&
				CurObj -> ElementSize == aElementSize &&
				CurObj -> InterpPoints == aInterpPoints )
			{
				break;
			}

			if( CurObj -> Next == NULL && ObjCount >= R8B_FRACBANK_CACHE_MAX )
			{
				if( CurObj -> RefCount == 0 )
				{
					// Delete the last filter which is not used.

					PrevObj -> Next = NULL;
					delete CurObj;
					ObjCount--;
				}
				else
				{
					// Move the last filter to the top of the list since it
					// seems to be in use for a long time.

					PrevObj -> Next = NULL;
					CurObj -> Next = Objects.unkeep();
					Objects = CurObj;
				}

				CurObj = NULL;
				break;
			}

			PrevObj = CurObj;
			CurObj = CurObj -> Next;
		}

		if( CurObj != NULL )
		{
			CurObj -> RefCount++;

			if( PrevObj == NULL )
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

			CurObj = new CDSPFracDelayFilterBank( aFilterLen, aFilterFracs,
				aElementSize, aInterpPoints );

			ObjCount++;
		}

		// Insert the filter at the start of the list.

		CurObj -> Next = Objects.unkeep();
		Objects = CurObj;

		return( *CurObj );
	}

private:
	static CSyncObject StateSync; ///< Cache state synchronizer.
		///<
	static CPtrKeeper< CDSPFracDelayFilterBank* > Objects; ///< The chain of
		///< cached objects.
		///<
	static int ObjCount; ///< The number of objects currently preset in the
		///< cache.
		///<
};

// ---------------------------------------------------------------------------
// CDSPFracDelayFilterBank PUBLIC
// ---------------------------------------------------------------------------

inline void CDSPFracDelayFilterBank :: unref()
{
	R8BSYNC( CDSPFracDelayFilterBankCache :: StateSync );

	RefCount--;
}

/**
 * @param l Number 1.
 * @param s Number 2.
 * @param[out] GCD Resulting GCD.
 * @return "True" if the greatest common denominator of 2 numbers was
 * found.
 */

inline bool findGCD( double l, double s, double& GCD )
{
	int it = 0;

	while( it < 50 )
	{
		if( s <= 0.0 )
		{
			GCD = l;
			return( true );
		}

		const double r = l - s;
		l = s;
		s = ( r < 0.0 ? -r : r );
		it++;
	}

	return( false );
}

/**
 * Function evaluates source and destination sample rate ratio and returns
 * the required input and output stepping. Function returns "false" if
 * *this class cannot be used to perform interpolation using these sample
 * rates.
 *
 * @param SSampleRate Source sample rate.
 * @param DSampleRate Destination sample rate.
 * @param[out] ResInStep Resulting input step.
 * @param[out] ResOutStep Resulting output step.
 * @return "True" if stepping was acquired.
 */

inline bool getWholeStepping( const double SSampleRate,
	const double DSampleRate, int& ResInStep, int& ResOutStep )
{
	double GCD;

	if( !findGCD( SSampleRate, DSampleRate, GCD ) || GCD < 1.0 )
	{
		return( false );
	}

	const double InStep0 = SSampleRate / GCD;
	ResInStep = (int) InStep0;
	const double OutStep0 = DSampleRate / GCD;
	ResOutStep = (int) OutStep0;

	if( InStep0 != ResInStep || OutStep0 != ResOutStep )
	{
		return( false );
	}

	if( ResOutStep > 1500 )
	{
		// Do not allow large output stepping due to low cache
		// performance of large filter banks.

		return( false );
	}

	return( true );
}

/**
 * @brief Fractional delay filter bank-based interpolator class.
 *
 * Class implements the fractional delay interpolator. This implementation at
 * first puts the input signal into a ring buffer and then performs
 * interpolation. The interpolation is performed using sinc-based fractional
 * delay filters. These filters are contained in a bank, and for higher
 * precision they are interpolated between adjacent filters.
 *
 * To increase sample timing precision, this class uses "resettable counter"
 * approach. This gives zero overall sample timing error. With the
 * R8B_FASTTIMING configuration option enabled, the sample timing experiences
 * a very minor drift.
 *
 * VERY IMPORTANT: the interpolation step should not exceed FilterLen / 2 + 1
 * samples or the algorithm in its current form will fail. However, this
 * condition can be easily met if the input signal is suitably downsampled
 * first before the interpolation is performed.
 *
 * @param FilterLen Specifies the number of samples (taps) each fractional
 * delay filter should have. See the r8b::CDSPFracDelayFilterBank class for
 * more details.
 * @param FilterFracs The number of fractional delay positions to sample. See
 * the r8b::CDSPFracDelayFilterBank class for more details.
 */

template< int FilterLen, int FilterFracs >
class CDSPFracInterpolator : public CDSPProcessor
{
public:
	/**
	 * Constructor initalizes the interpolator. It is important to call the
	 * getMaxOutLen() function afterwards to obtain the optimal output buffer
	 * length.
	 *
	 * @param aSrcSampleRate Source sample rate.
	 * @param aDstSampleRate Destination sample rate.
	 * @param PrevLatency Latency, in samples (any value >=0), which was left
	 * in the output signal by a previous process. This latency will be
	 * consumed completely.
	 */

	CDSPFracInterpolator( const double aSrcSampleRate,
		const double aDstSampleRate, const double PrevLatency )
		: SrcSampleRate( aSrcSampleRate )
		, DstSampleRate( aDstSampleRate )
	#if R8B_FASTTIMING
		, FracStep( aSrcSampleRate / aDstSampleRate )
	#endif // R8B_FASTTIMING
	#if R8B_FLTTEST
		, FilterBank( FilterLen, FilterFracs, 4, 8 )
	#endif // R8B_FLTTEST
	{
		R8BASSERT( SrcSampleRate > 0.0 );
		R8BASSERT( DstSampleRate > 0.0 );
		R8BASSERT( PrevLatency >= 0.0 );
		R8BASSERT( BufLenBits >= 5 );
		R8BASSERT(( 1 << BufLenBits ) >= FilterLen * 3 );

		InitFracPos = PrevLatency;
		Latency = (int) InitFracPos;
		InitFracPos -= Latency;

		if( getWholeStepping( SrcSampleRate, DstSampleRate, InStep, OutStep ))
		{
			InitFracPosW = (int) ( InitFracPos * OutStep );
			LatencyFrac = InitFracPos - (double) InitFracPosW / OutStep;
			FilterBankW = &CDSPFracDelayFilterBankCache :: getFilterBank(
				FilterLen, OutStep, 1, 2 );
		}
		else
		{
			InitFracPosW = 0;
			LatencyFrac = 0.0;
			FilterBankW = NULL;
		}

		clear();
	}

	virtual ~CDSPFracInterpolator()
	{
		if( FilterBankW != NULL )
		{
			FilterBankW -> unref();
		}
	}

	virtual int getLatency() const
	{
		return( 0 );
	}

	virtual double getLatencyFrac() const
	{
		return( LatencyFrac );
	}

	virtual int getMaxOutLen( const int MaxInLen ) const
	{
		R8BASSERT( MaxInLen >= 0 );

		return( (int) ceil( MaxInLen * DstSampleRate / SrcSampleRate ) + 1 );
	}

	virtual void clear()
	{
		LatencyLeft = Latency;
		BufLeft = 0;
		WritePos = 0;
		ReadPos = BufLen - fll; // Set "read" position to account for filter's
			// latency at zero fractional delay.

		memset( &Buf[ ReadPos ], 0, fll * sizeof( double ));

		if( FilterBankW != NULL )
		{
			InPosFracW = InitFracPosW;
		}
		else
		{
			InPosFrac = InitFracPos;

			#if !R8B_FASTTIMING
				InCounter = 0;
				InPosInt = 0;
				InPosShift = InitFracPos * DstSampleRate / SrcSampleRate;
			#endif // !R8B_FASTTIMING
		}
	}

	virtual int process( double* ip, int l, double*& op0 )
	{
		R8BASSERT( l >= 0 );
		R8BASSERT( ip != op0 || l == 0 || SrcSampleRate > DstSampleRate );

		if( LatencyLeft > 0 )
		{
			if( LatencyLeft >= l )
			{
				LatencyLeft -= l;
				return( 0 );
			}

			l -= LatencyLeft;
			ip += LatencyLeft;
			LatencyLeft = 0;
		}

		double* op = op0;

		while( l > 0 )
		{
			// Add new input samples to both halves of the ring buffer.

			const int b = min( min( l, BufLen - WritePos ),
				BufLeftMax - BufLeft );

			double* const wp1 = Buf + WritePos;
			memcpy( wp1, ip, b * sizeof( double ));

			if( WritePos < flo )
			{
				const int c = min( b, flo - WritePos );
				memcpy( wp1 + BufLen, wp1, c * sizeof( double ));
			}

			ip += b;
			WritePos = ( WritePos + b ) & BufLenMask;
			l -= b;
			BufLeft += b;
			int i;

			// Produce as many output samples as possible.

			if( FilterBankW != NULL )
			{
				// Whole-number stepping.

				while( BufLeft > fl2 )
				{
					const double* const ftp = &(*FilterBankW)[ InPosFracW ];
					const double* const rp = Buf + ReadPos;
					double s = 0.0;

					for( i = 0; i < FilterLen; i++ )
					{
						s += ftp[ i ] * rp[ i ];
					}

					*op = s;
					op++;

					InPosFracW += InStep;
					const int PosIncr = InPosFracW / OutStep;
					InPosFracW -= PosIncr * OutStep;

					ReadPos = ( ReadPos + PosIncr ) & BufLenMask;
					BufLeft -= PosIncr;
				}
			}
			else
			{
				while( BufLeft > fl2 )
				{
					double x = InPosFrac * FilterFracs;
					const int fti = (int) x; // Function table index.
					x -= fti; // Coefficient for interpolation between
						// adjacent fractional delay filters.
					const double x2 = x * x;
					const double* const ftp = &FilterBank[ fti ];
					const double* const rp = Buf + ReadPos;
					double s = 0.0;
					int ii = 0;

				#if R8B_FLTTEST
					const double x3 = x2 * x;
				#endif // R8B_FLTTEST

					for( i = 0; i < FilterLen; i++ )
					{
					#if !R8B_FLTTEST
						s += ( ftp[ ii ] + ftp[ ii + 1 ] * x +
							ftp[ ii + 2 ] * x2 ) * rp[ i ];
					#else // !R8B_FLTTEST
						s += ( ftp[ ii ] + ftp[ ii + 1 ] * x +
							ftp[ ii + 2 ] * x2 + ftp[ ii + 3 ] * x3 ) *
							rp[ i ];
					#endif // !R8B_FLTTEST

						ii += FilterElementSize;
					}

					*op = s;
					op++;

					#if R8B_FASTTIMING

						InPosFrac += FracStep;
						const int PosIncr = (int) InPosFrac;
						InPosFrac -= PosIncr;

					#else // R8B_FASTTIMING

						InCounter++;
						const double NextInPos = ( InCounter + InPosShift ) *
							SrcSampleRate / DstSampleRate;

						const int NextInPosInt = (int) NextInPos;
						const int PosIncr = NextInPosInt - InPosInt;
						InPosInt = NextInPosInt;
						InPosFrac = NextInPos - NextInPosInt;

					#endif // R8B_FASTTIMING

					ReadPos = ( ReadPos + PosIncr ) & BufLenMask;
					BufLeft -= PosIncr;
				}
			}
		}

		#if !R8B_FASTTIMING

			if( FilterBankW == NULL && InCounter > 1000 )
			{
				// Reset the interpolation position counter to achieve a
				// higher sample timing precision.

				InCounter = 0;
				InPosInt = 0;
				InPosShift = InPosFrac * DstSampleRate / SrcSampleRate;
			}

		#endif // !R8B_FASTTIMING

		return( (int) ( op - op0 ));
	}

private:
#if !R8B_FLTTEST
	static const int FilterElementSize = 3; ///< The number of "doubles" a
		///< single filter tap consists of (includes interpolation
		///< coefficients).
		///<
#else // !R8B_FLTTEST
	static const int FilterElementSize = 4; ///< The number of "doubles" a
		///< single filter tap consists of (includes interpolation
		///< coefficients). During filter testing a higher precision
		///< interpolation is used.
		///<
#endif // !R8B_FLTTEST

	static const int fl2 = FilterLen >> 1; ///< Right-side (half) filter
		///< length.
		///<
	static const int fll = fl2 - 1; ///< Input latency.
		///<
	static const int flo = fll + fl2; ///< Overrun length.
		///<
	static const int BufLenBits = 8; ///< The length of the ring buffer,
		///< expressed as Nth power of 2. This value can be reduced if it is
		///< known that only short input buffers will be passed to the
		///< interpolator. The minimum value of this parameter is 5, and
		///< 1 << BufLenBits should be at least 3 times larger than the
		///< FilterLen.
		///<
	static const int BufLen = 1 << BufLenBits; ///< The length of the ring
		///< buffer. The actual length is twice as long to allow "beyond max
		///< position" positioning.
		///<
	static const int BufLenMask = BufLen - 1; ///< Mask used for quick buffer
		///< position wrapping.
		///<
	static const int BufLeftMax = BufLen - fll; ///< The number of new samples
		///< that the ring buffer can hold at most.
		///<
	double Buf[ BufLen + flo ]; ///< The ring buffer, including overrun
		///< protection.
		///<
	double SrcSampleRate; ///< Source sample rate.
		///<
	double DstSampleRate; ///< Destination sample rate.
		///<
	int InStep; ///< Input whole-number stepping.
		///<
	int OutStep; ///< Output whole-number stepping (corresponds to filter bank
		///< size).
		///<
	double InitFracPos; ///< Initial fractional position, in samples, in the
		///< range [0; 1).
		///<
	int InitFracPosW; ///< Initial fractional position for whole-number
		///< stepping.
		///<
	int Latency; ///< Initial latency that should be removed from the input.
		///<
	double LatencyFrac; ///< Left-over fractional latency.
		///<
	int BufLeft; ///< The number of samples left in the buffer to process.
		///<
	int WritePos; ///< The current buffer write position. Incremented together
		///< with the BufLeft variable.
		///<
	int ReadPos; ///< The current buffer read position.
		///<
	int LatencyLeft; ///< Input latency left to remove.
		///<
	double InPosFrac; ///< Interpolation position (fractional part).
		///<
	int InPosFracW; ///< Interpolation position (fractional part) for
		///< whole-number stepping. Corresponds to the index into the filter
		///< bank.
		///<
	CDSPFracDelayFilterBank* FilterBankW; ///< Whole-number stepping filter
		///< bank. NULL if whole-number stepping is not in use.
		///<
#if R8B_FASTTIMING
	double FracStep; ///< Fractional sample timing step.
#else // R8B_FASTTIMING
	int InCounter; ///< Interpolation step counter.
		///<
	int InPosInt; ///< Interpolation position (integer part).
		///<
	double InPosShift; ///< Interpolation position fractional shift.
		///<
#endif // R8B_FASTTIMING

#if !R8B_FLTTEST
	static const CDSPFracDelayFilterBank FilterBank; ///< Filter bank object,
		///< defined statically if no filter test is taking place.
		///<
#else // !R8B_FLTTEST

public:
	CDSPFracDelayFilterBank FilterBank; ///< Filter bank object, defined as a
		///< member variable to allow for recalculation.
		///<
#endif // !R8B_FLTTEST
};

// ---------------------------------------------------------------------------

#if !R8B_FLTTEST
template< int FilterLen, int FilterFracs >
const CDSPFracDelayFilterBank CDSPFracInterpolator<
	FilterLen, FilterFracs > :: FilterBank( FilterLen, FilterFracs, 3, 8 );
#endif // !R8B_FLTTEST

// ---------------------------------------------------------------------------

} // namespace r8b

#endif // R8B_CDSPFRACINTERPOLATOR_INCLUDED
