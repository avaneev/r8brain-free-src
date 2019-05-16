//$ nocpp

/**
 * @file CDSPHBDownsampler.h
 *
 * @brief Half-band downsampling convolver class.
 *
 * This file includes half-band downsampling convolver class.
 *
 * r8brain-free-src Copyright (c) 2019 Aleksey Vaneev
 * See the "License.txt" file for license.
 */

#ifndef R8B_CDSPHBDOWNSAMPLER_INCLUDED
#define R8B_CDSPHBDOWNSAMPLER_INCLUDED

#include "CDSPProcessor.h"

namespace r8b {

/**
 * @brief Half-band downsampler class.
 *
 * Class implements brute-force half-band 2X downsampling that uses small
 * sparse symmetric FIR filters. The output has 2.0 gain.
 */

class CDSPHBDownsampler : public CDSPProcessor
{
public:
	/**
	 * Constructor initalizes the half-band downsampler.
	 *
	 * @param ReqAtten Required half-band filter attentuation.
	 */

	CDSPHBDownsampler( const double ReqAtten )
	{
		static const int FltCount = 11;
		static const double HBKernel_4[ 4 ] = { // att -64.9241 dB, frac 4.0
			6.1073830069265711e-001, -1.4463982571471876e-001,
			4.1136036923118187e-002, -7.4740105856914872e-003 };
		static const double HBKernel_5[ 5 ] = { // att -87.4775 dB, frac 4.0
			6.1553054142504338e-001, -1.5591352339398118e-001,
			5.2404802661298266e-002, -1.4230574348726146e-002,
			2.2457593805377831e-003 };
		static const double HBKernel_6[ 6 ] = { // att -104.5154 dB, frac 4.0
			6.1883766561934184e-001, -1.6396282655874558e-001,
			6.1104571279325129e-002, -2.0317756445217543e-002,
			5.0264527466018826e-003, -6.9392938429507279e-004 };
		static const double HBKernel_7[ 7 ] = { // att -120.6199 dB, frac 4.0
			6.2125313688727779e-001, -1.6999763849273491e-001,
			6.8014108060738196e-002, -2.5679821316697125e-002,
			7.9798828249699784e-003, -1.7871060154498470e-003,
			2.1836606459564009e-004 };
		static const double HBKernel_8[ 8 ] = { // att -136.5151 dB, frac 4.0
			6.2309299085367287e-001, -1.7468969193368433e-001,
			7.3628746444973150e-002, -3.0378268550055314e-002,
			1.0908085227657214e-002, -3.1287343330312556e-003,
			6.3632014609722092e-004, -6.9597139145649502e-005 };
		static const double HBKernel_9[ 9 ] = { // att -152.3240 dB, frac 4.0
			6.2454069594794803e-001, -1.7844303649890664e-001,
			7.8279410808762842e-002, -3.4501119561829857e-002,
			1.3717889826645487e-002, -4.6090109007760798e-003,
			1.2192752061406873e-003, -2.2647618541786664e-004,
			2.2395554542567748e-005 };
		static const double HBKernel_10[ 10 ] = { // att -168.0859 dB, frac 4.0
			6.2570883988448611e-001, -1.8151274643053061e-001,
			8.2191863294185458e-002, -3.8131779329357615e-002,
			1.6367492549512565e-002, -6.1530178832078578e-003,
			1.9277693942420303e-003, -4.7165916432255402e-004,
			8.0491894752808465e-005, -7.2581515842465856e-006 };
		static const double HBKernel_11[ 11 ] = { // att -183.7962 dB, frac 4.0
			6.2667167706646965e-001, -1.8407153341833782e-001,
			8.5529995600327216e-002, -4.1346831452173063e-002,
			1.8844831683400131e-002, -7.7125170314919214e-003,
			2.7268674834296570e-003, -7.9745028391855826e-004,
			1.8116344571770699e-004, -2.8569149678673122e-005,
			2.3667021922861879e-006 };
		static const double HBKernel_12[ 12 ] = { // att -199.4768 dB, frac 4.0
			6.2747849729182659e-001, -1.8623616781335248e-001,
			8.8409755856508648e-002, -4.4207468780136254e-002,
			2.1149175912217915e-002, -9.2551508154301194e-003,
			3.5871562052326249e-003, -1.1923167600753576e-003,
			3.2627812001613326e-004, -6.9106902008709490e-005,
			1.0122897772888322e-005, -7.7531878091292963e-007 };
		static const double HBKernel_13[ 13 ] = { // att -215.1364 dB, frac 4.0
			6.2816416238782957e-001, -1.8809076918442266e-001,
			9.0918539368474965e-002, -4.6765502172995604e-002,
			2.3287520069933797e-002, -1.0760626940880943e-002,
			4.4853921118213676e-003, -1.6438774496992904e-003,
			5.1441308429384374e-004, -1.3211724349740752e-004,
			2.6191316362108199e-005, -3.5802424384280469e-006,
			2.5491272423372411e-007 };
		static const double HBKernel_14[ 14 ] = { // att -230.7526 dB, frac 4.0
			6.2875473147254901e-001, -1.8969942008858576e-001,
			9.3126095475258408e-002, -4.9067252227455962e-002,
			2.5273009767563311e-002, -1.2218646838258702e-002,
			5.4048946497798353e-003, -2.1409921992386689e-003,
			7.4250304371305991e-004, -2.1924546773651068e-004,
			5.3015823597863675e-005, -9.8743070771832892e-006,
			1.2650397198764347e-006, -8.4146728313072455e-008 };
		static const double FltAttens[ FltCount ] = {
			64.9241, 87.4775, 104.5154, 120.6199, 136.5151, 152.3240,
			168.0859, 183.7962, 199.4768, 215.1364, 230.7526 };
		static const double* const FltPtrs[ FltCount ] = { HBKernel_4,
			HBKernel_5, HBKernel_6, HBKernel_7, HBKernel_8, HBKernel_9,
			HBKernel_10, HBKernel_11, HBKernel_12, HBKernel_13,
			HBKernel_14 };
		static const CConvolveFn FltConvFn[ FltCount ] = {
			&CDSPHBDownsampler :: convolve4, &CDSPHBDownsampler :: convolve5,
			&CDSPHBDownsampler :: convolve6, &CDSPHBDownsampler :: convolve7,
			&CDSPHBDownsampler :: convolve8, &CDSPHBDownsampler :: convolve9,
			&CDSPHBDownsampler :: convolve10, &CDSPHBDownsampler :: convolve11,
			&CDSPHBDownsampler :: convolve12, &CDSPHBDownsampler :: convolve13,
			&CDSPHBDownsampler :: convolve14 };

		int k = 0;

		while( k != FltCount - 1 && FltAttens[ k ] < ReqAtten )
		{
			k++;
		}

		flt = FltPtrs[ k ];
		convfn = FltConvFn[ k ];
		const int fltt = 4 + k;
		fll = fltt * 2 - 1;
		fl2 = fll;
		fl4 = fll + fl2;

		R8BCONSOLE( "CDSPHBDownsampler: taps=%i att=%.1f io=1/2\n", fltt,
			FltAttens[ k ]);

		clear();
	}

	virtual int getLatency() const
	{
		return( 0 );
	}

	virtual double getLatencyFrac() const
	{
		return( 0.0 );
	}

	virtual int getMaxOutLen( const int MaxInLen ) const
	{
		R8BASSERT( MaxInLen >= 0 );

		return(( MaxInLen + 1 ) / 2 );
	}

	/**
	 * Function clears (resets) the state of *this object and returns it to
	 * the state after construction. All input data accumulated in the
	 * internal buffer so far will be discarded.
	 */

	virtual void clear()
	{
		BufLeft = 0;
		WritePos = 0;
		ReadPos = BufLen - fll; // Set "read" position to
			// account for filter's latency.

		memset( &Buf[ ReadPos ], 0, fll * sizeof( double ));
	}

	virtual int process( double* ip, int l, double*& op0 )
	{
		R8BASSERT( l >= 0 );

		double* op = op0;

		while( l > 0 )
		{
			// Add new input samples to both halves of the ring buffer.

			const int b = min( min( l, BufLen - WritePos ),
				BufLen - fll - BufLeft );

			double* const wp1 = Buf + WritePos;
			memcpy( wp1, ip, b * sizeof( double ));

			if( WritePos < fl4 )
			{
				const int c = min( b, fl4 - WritePos );
				memcpy( wp1 + BufLen, wp1, c * sizeof( double ));
			}

			ip += b;
			WritePos = ( WritePos + b ) & BufLenMask;
			l -= b;
			BufLeft += b;

			if( BufLeft > fl2 )
			{
				const int c = ( BufLeft - fl2 + 1 ) >> 1;

				( *this.*convfn )( op, c );

				op += c;
				BufLeft -= c + c;
			}
		}

		return( (int) ( op - op0 ));
	}

private:
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
	double Buf[ BufLen + 54 ]; ///< The ring buffer, including overrun
		///< protection for the largest filter.
		///<
	const double* flt; ///< Half-band filter taps.
		///<
	int fll; ///< Input latency.
		///<
	int fl2; ///< Right-side filter length.
		///<
	int fl4; ///< Overrun length.
		///<
	int BufLeft; ///< The number of samples left in the buffer to process.
		///< When this value is below FilterLenD2Plus1, the interpolation
		///< cycle ends.
		///<
	int WritePos; ///< The current buffer write position. Incremented together
		///< with the BufLeft variable.
		///<
	int ReadPos; ///< The current buffer read position.
		///<
	typedef void( CDSPHBDownsampler :: *CConvolveFn )( double* op,
		int c ); ///< Convolution funtion type.
		///<
	CConvolveFn convfn; ///< Convolution function in use.
		///<

#define R8BHBC1( fn ) \
	void fn( double* op, int c ) \
	{ \
		const double* const rp0 = Buf + fll; \
		int rpos = ReadPos; \
		while( c > 0 ) \
		{ \
			const double* const rp = rp0 + rpos; \
			*op = rp[ 0 ] +

#define R8BHBC2 \
			rpos = ( rpos + 2 ) & BufLenMask; \
			op++; \
			c--; \
		} \
		ReadPos = rpos; \
	}

	R8BHBC1( convolve4 )
				flt[ 0 ] * ( rp[ 1 ] + rp[ -1 ]) +
				flt[ 1 ] * ( rp[ 3 ] + rp[ -3 ]) +
				flt[ 2 ] * ( rp[ 5 ] + rp[ -5 ]) +
				flt[ 3 ] * ( rp[ 7 ] + rp[ -7 ]);
	R8BHBC2

	R8BHBC1( convolve5 )
				flt[ 0 ] * ( rp[ 1 ] + rp[ -1 ]) +
				flt[ 1 ] * ( rp[ 3 ] + rp[ -3 ]) +
				flt[ 2 ] * ( rp[ 5 ] + rp[ -5 ]) +
				flt[ 3 ] * ( rp[ 7 ] + rp[ -7 ]) +
				flt[ 4 ] * ( rp[ 9 ] + rp[ -9 ]);
	R8BHBC2

	R8BHBC1( convolve6 )
				flt[ 0 ] * ( rp[ 1 ] + rp[ -1 ]) +
				flt[ 1 ] * ( rp[ 3 ] + rp[ -3 ]) +
				flt[ 2 ] * ( rp[ 5 ] + rp[ -5 ]) +
				flt[ 3 ] * ( rp[ 7 ] + rp[ -7 ]) +
				flt[ 4 ] * ( rp[ 9 ] + rp[ -9 ]) +
				flt[ 5 ] * ( rp[ 11 ] + rp[ -11 ]);
	R8BHBC2

	R8BHBC1( convolve7 )
				flt[ 0 ] * ( rp[ 1 ] + rp[ -1 ]) +
				flt[ 1 ] * ( rp[ 3 ] + rp[ -3 ]) +
				flt[ 2 ] * ( rp[ 5 ] + rp[ -5 ]) +
				flt[ 3 ] * ( rp[ 7 ] + rp[ -7 ]) +
				flt[ 4 ] * ( rp[ 9 ] + rp[ -9 ]) +
				flt[ 5 ] * ( rp[ 11 ] + rp[ -11 ]) +
				flt[ 6 ] * ( rp[ 13 ] + rp[ -13 ]);
	R8BHBC2

	R8BHBC1( convolve8 )
				flt[ 0 ] * ( rp[ 1 ] + rp[ -1 ]) +
				flt[ 1 ] * ( rp[ 3 ] + rp[ -3 ]) +
				flt[ 2 ] * ( rp[ 5 ] + rp[ -5 ]) +
				flt[ 3 ] * ( rp[ 7 ] + rp[ -7 ]) +
				flt[ 4 ] * ( rp[ 9 ] + rp[ -9 ]) +
				flt[ 5 ] * ( rp[ 11 ] + rp[ -11 ]) +
				flt[ 6 ] * ( rp[ 13 ] + rp[ -13 ]) +
				flt[ 7 ] * ( rp[ 15 ] + rp[ -15 ]);
	R8BHBC2

	R8BHBC1( convolve9 )
				flt[ 0 ] * ( rp[ 1 ] + rp[ -1 ]) +
				flt[ 1 ] * ( rp[ 3 ] + rp[ -3 ]) +
				flt[ 2 ] * ( rp[ 5 ] + rp[ -5 ]) +
				flt[ 3 ] * ( rp[ 7 ] + rp[ -7 ]) +
				flt[ 4 ] * ( rp[ 9 ] + rp[ -9 ]) +
				flt[ 5 ] * ( rp[ 11 ] + rp[ -11 ]) +
				flt[ 6 ] * ( rp[ 13 ] + rp[ -13 ]) +
				flt[ 7 ] * ( rp[ 15 ] + rp[ -15 ]) +
				flt[ 8 ] * ( rp[ 17 ] + rp[ -17 ]);
	R8BHBC2

	R8BHBC1( convolve10 )
				flt[ 0 ] * ( rp[ 1 ] + rp[ -1 ]) +
				flt[ 1 ] * ( rp[ 3 ] + rp[ -3 ]) +
				flt[ 2 ] * ( rp[ 5 ] + rp[ -5 ]) +
				flt[ 3 ] * ( rp[ 7 ] + rp[ -7 ]) +
				flt[ 4 ] * ( rp[ 9 ] + rp[ -9 ]) +
				flt[ 5 ] * ( rp[ 11 ] + rp[ -11 ]) +
				flt[ 6 ] * ( rp[ 13 ] + rp[ -13 ]) +
				flt[ 7 ] * ( rp[ 15 ] + rp[ -15 ]) +
				flt[ 8 ] * ( rp[ 17 ] + rp[ -17 ]) +
				flt[ 9 ] * ( rp[ 19 ] + rp[ -19 ]);
	R8BHBC2

	R8BHBC1( convolve11 )
				flt[ 0 ] * ( rp[ 1 ] + rp[ -1 ]) +
				flt[ 1 ] * ( rp[ 3 ] + rp[ -3 ]) +
				flt[ 2 ] * ( rp[ 5 ] + rp[ -5 ]) +
				flt[ 3 ] * ( rp[ 7 ] + rp[ -7 ]) +
				flt[ 4 ] * ( rp[ 9 ] + rp[ -9 ]) +
				flt[ 5 ] * ( rp[ 11 ] + rp[ -11 ]) +
				flt[ 6 ] * ( rp[ 13 ] + rp[ -13 ]) +
				flt[ 7 ] * ( rp[ 15 ] + rp[ -15 ]) +
				flt[ 8 ] * ( rp[ 17 ] + rp[ -17 ]) +
				flt[ 9 ] * ( rp[ 19 ] + rp[ -19 ]) +
				flt[ 10 ] * ( rp[ 21 ] + rp[ -21 ]);
	R8BHBC2

	R8BHBC1( convolve12 )
				flt[ 0 ] * ( rp[ 1 ] + rp[ -1 ]) +
				flt[ 1 ] * ( rp[ 3 ] + rp[ -3 ]) +
				flt[ 2 ] * ( rp[ 5 ] + rp[ -5 ]) +
				flt[ 3 ] * ( rp[ 7 ] + rp[ -7 ]) +
				flt[ 4 ] * ( rp[ 9 ] + rp[ -9 ]) +
				flt[ 5 ] * ( rp[ 11 ] + rp[ -11 ]) +
				flt[ 6 ] * ( rp[ 13 ] + rp[ -13 ]) +
				flt[ 7 ] * ( rp[ 15 ] + rp[ -15 ]) +
				flt[ 8 ] * ( rp[ 17 ] + rp[ -17 ]) +
				flt[ 9 ] * ( rp[ 19 ] + rp[ -19 ]) +
				flt[ 10 ] * ( rp[ 21 ] + rp[ -21 ]) +
				flt[ 11 ] * ( rp[ 23 ] + rp[ -23 ]);
	R8BHBC2

	R8BHBC1( convolve13 )
				flt[ 0 ] * ( rp[ 1 ] + rp[ -1 ]) +
				flt[ 1 ] * ( rp[ 3 ] + rp[ -3 ]) +
				flt[ 2 ] * ( rp[ 5 ] + rp[ -5 ]) +
				flt[ 3 ] * ( rp[ 7 ] + rp[ -7 ]) +
				flt[ 4 ] * ( rp[ 9 ] + rp[ -9 ]) +
				flt[ 5 ] * ( rp[ 11 ] + rp[ -11 ]) +
				flt[ 6 ] * ( rp[ 13 ] + rp[ -13 ]) +
				flt[ 7 ] * ( rp[ 15 ] + rp[ -15 ]) +
				flt[ 8 ] * ( rp[ 17 ] + rp[ -17 ]) +
				flt[ 9 ] * ( rp[ 19 ] + rp[ -19 ]) +
				flt[ 10 ] * ( rp[ 21 ] + rp[ -21 ]) +
				flt[ 11 ] * ( rp[ 23 ] + rp[ -23 ]) +
				flt[ 12 ] * ( rp[ 25 ] + rp[ -25 ]);
	R8BHBC2

	R8BHBC1( convolve14 )
				flt[ 0 ] * ( rp[ 1 ] + rp[ -1 ]) +
				flt[ 1 ] * ( rp[ 3 ] + rp[ -3 ]) +
				flt[ 2 ] * ( rp[ 5 ] + rp[ -5 ]) +
				flt[ 3 ] * ( rp[ 7 ] + rp[ -7 ]) +
				flt[ 4 ] * ( rp[ 9 ] + rp[ -9 ]) +
				flt[ 5 ] * ( rp[ 11 ] + rp[ -11 ]) +
				flt[ 6 ] * ( rp[ 13 ] + rp[ -13 ]) +
				flt[ 7 ] * ( rp[ 15 ] + rp[ -15 ]) +
				flt[ 8 ] * ( rp[ 17 ] + rp[ -17 ]) +
				flt[ 9 ] * ( rp[ 19 ] + rp[ -19 ]) +
				flt[ 10 ] * ( rp[ 21 ] + rp[ -21 ]) +
				flt[ 11 ] * ( rp[ 23 ] + rp[ -23 ]) +
				flt[ 12 ] * ( rp[ 25 ] + rp[ -25 ]) +
				flt[ 13 ] * ( rp[ 27 ] + rp[ -27 ]);
	R8BHBC2

#undef R8BHBC1
#undef R8BHBC2
};

// ---------------------------------------------------------------------------

} // namespace r8b

#endif // R8B_CDSPHBDOWNSAMPLER_INCLUDED
