//$ nocpp

/**
 * @file CDSPHBConvolver.h
 *
 * @brief Half-band downsampling convolver class.
 *
 * This file includes half-band downsampling convolver class.
 *
 * r8brain-free-src Copyright (c) 2019 Aleksey Vaneev
 * See the "License.txt" file for license.
 */

#ifndef R8B_CDSPHBCONVOLVER_INCLUDED
#define R8B_CDSPHBCONVOLVER_INCLUDED

#include "CDSPProcessor.h"

namespace r8b {

/**
 * @brief Half-band downsampling convolver class.
 *
 * Class implements brute-force half-band 2X downsampling that uses small
 * sparse symmetric FIR filters.
 */

class CDSPHBConvolver : public CDSPProcessor
{
public:
	/**
	 * Constructor initalizes the convolver.
	 *
	 * @param ReqAtten Required half-band filter attentuation.
	 */

	CDSPHBConvolver( const double ReqAtten )
	{
		static const int FltCount = 11;
		static const double HBKernel_4[ 4 ] = { // StopAtten = -64.9241 dB
			3.0536915034632850e-001, -7.2319912857359170e-002,
			2.0568018461558868e-002, -3.7370052928455966e-003 };
		static const double HBKernel_5[ 5 ] = { // StopAtten = -87.4775 dB
			3.0776527071252091e-001, -7.7956761696988730e-002,
			2.6202401330647218e-002, -7.1152871743619395e-003,
			1.1228796902684746e-003 };
		static const double HBKernel_6[ 6 ] = { // StopAtten = -104.5154 dB
			3.0941883280967319e-001, -8.1981413279378135e-002,
			3.0552285639668053e-002, -1.0158878222612300e-002,
			2.5132263733023130e-003, -3.4696469214778869e-004 };
		static const double HBKernel_7[ 7 ] = { // StopAtten = -120.6199 dB
			3.1062656844362180e-001, -8.4998819246324686e-002,
			3.4007054030319894e-002, -1.2839910658310139e-002,
			3.9899414124635602e-003, -8.9355300771687407e-004,
			1.0918303229612292e-004 };
		static const double HBKernel_8[ 8 ] = { // StopAtten = -136.5151 dB
			3.1154649542685675e-001, -8.7344845966894846e-002,
			3.6814373222551655e-002, -1.5189134275084837e-002,
			5.4540426138667275e-003, -1.5643671665347414e-003,
			3.1816007305529581e-004, -3.4798569574112648e-005 };
		static const double HBKernel_9[ 9 ] = { // StopAtten = -152.3240 dB
			3.1227034797423214e-001, -8.9221518250120035e-002,
			3.9139705405200953e-002, -1.7250559781632369e-002,
			6.8589449138030112e-003, -2.3045054506348750e-003,
			6.0963760316434801e-004, -1.1323809273311447e-004,
			1.1197777274540854e-005 };
		static const double HBKernel_10[ 10 ] = { // StopAtten = -168.0859 dB
			3.1285441994054719e-001, -9.0756373210797991e-002,
			4.1095931641371264e-002, -1.9065889659329205e-002,
			8.1837462708107545e-003, -3.0765089392719158e-003,
			9.6388469603070348e-004, -2.3582958177437163e-004,
			4.0245947281469796e-005, -3.6290757795542916e-006 };
		static const double HBKernel_11[ 11 ] = { // StopAtten = -183.7962 dB
			3.1333583852618196e-001, -9.2035766690425125e-002,
			4.2764997775709412e-002, -2.0673415702531672e-002,
			9.4224158235471922e-003, -3.8562585043157016e-003,
			1.3634337358552351e-003, -3.9872513956817901e-004,
			9.0581722118621144e-005, -1.4284574681888706e-005,
			1.1833510784597934e-006 };
		static const double HBKernel_12[ 12 ] = { // StopAtten = -199.4768 dB
			3.1373924858707514e-001, -9.3118083749112346e-002,
			4.4204877719434599e-002, -2.2103734183864354e-002,
			1.0574587791389645e-002, -4.6275752986798535e-003,
			1.7935780427078551e-003, -5.9615835307380571e-004,
			1.6313905035447718e-004, -3.4553448402641074e-005,
			5.0614484139455129e-006, -3.8765934695650373e-007 };
		static const double HBKernel_13[ 13 ] = { // StopAtten = -215.1363 dB
			3.1408208078360489e-001, -9.4045383482463718e-002,
			4.5459268183225543e-002, -2.3382749556300153e-002,
			1.1643758756137546e-002, -5.3803125703677512e-003,
			2.2426955190719688e-003, -8.2193845506672858e-004,
			2.5720642977944408e-004, -6.6058584117012487e-005,
			1.3095648579319885e-005, -1.7901195397129754e-006,
			1.2745620847631756e-007 };
		static const double HBKernel_14[ 14 ] = { // StopAtten = -230.7517 dB
			3.1437736642800962e-001, -9.4849711931384806e-002,
			4.6563050335205025e-002, -2.4533628834576215e-002,
			1.2636507245178175e-002, -6.1093251664958714e-003,
			2.7024484370059554e-003, -1.0704967072912019e-003,
			3.7125180403951852e-004, -1.0962284308879688e-004,
			2.6507945896536711e-005, -4.9371616503979254e-006,
			6.3252117721290801e-007, -4.2073474482151596e-008 };
		static const double FltAttens[ FltCount ] = {
			64.9241, 87.4775, 104.5154, 120.6199, 136.5151, 152.3240,
			168.0859, 183.7962, 199.4768, 215.1363, 230.7517 };
		static const double* const FltPtrs[ FltCount ] = { HBKernel_4,
			HBKernel_5, HBKernel_6, HBKernel_7, HBKernel_8, HBKernel_9,
			HBKernel_10, HBKernel_11, HBKernel_12, HBKernel_13,
			HBKernel_14 };
		static const CConvolveFn FltConvFn[ FltCount ] = {
			&CDSPHBConvolver :: convolve4, &CDSPHBConvolver :: convolve5,
			&CDSPHBConvolver :: convolve6, &CDSPHBConvolver :: convolve7,
			&CDSPHBConvolver :: convolve8, &CDSPHBConvolver :: convolve9,
			&CDSPHBConvolver :: convolve10, &CDSPHBConvolver :: convolve11,
			&CDSPHBConvolver :: convolve12, &CDSPHBConvolver :: convolve13,
			&CDSPHBConvolver :: convolve14 };

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

		R8BCONSOLE( "CDSPHBConvolver: taps=%i att=%.1f io=1/2\n", fltt,
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
				BufLen - fl2 - BufLeft );

			double* const wp1 = Buf + WritePos;

			if( WritePos < fl4 )
			{
				const int c = min( b, fl4 - WritePos );
				double* const wp2 = wp1 + BufLen;
				int i;

				for( i = 0; i < c; i++ )
				{
					wp1[ i ] = ip[ i ];
					wp2[ i ] = ip[ i ];
				}

				memcpy( wp1 + c, ip + c, ( b - c ) * sizeof( double ));
			}
			else
			{
				memcpy( wp1, ip, b * sizeof( double ));
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
	typedef void( CDSPHBConvolver :: *CConvolveFn )( double* op, int c ); ///<
		///< Convolution funtion type.
		///<
	CConvolveFn convfn; ///< Convolution function in use.
		///<

#define R8BHBC1( fn ) \
	void fn( double* op, int c ) \
	{ \
		const double* const rp0 = Buf + fl2; \
		int rpos = ReadPos; \
		while( c > 0 ) \
		{ \
			const double* const rp = rp0 + rpos;
#define R8BHBC2 \
			rpos = ( rpos + 2 ) & BufLenMask; \
			op++; \
			c--; \
		} \
		ReadPos = rpos; \
	}

	R8BHBC1( convolve4 )
			*op = 0.5 * rp[ 0 ] +
				flt[ 0 ] * ( rp[ 1 ] + rp[ -1 ]) +
				flt[ 1 ] * ( rp[ 3 ] + rp[ -3 ]) +
				flt[ 2 ] * ( rp[ 5 ] + rp[ -5 ]) +
				flt[ 3 ] * ( rp[ 7 ] + rp[ -7 ]);
	R8BHBC2

	R8BHBC1( convolve5 )
			*op = 0.5 * rp[ 0 ] +
				flt[ 0 ] * ( rp[ 1 ] + rp[ -1 ]) +
				flt[ 1 ] * ( rp[ 3 ] + rp[ -3 ]) +
				flt[ 2 ] * ( rp[ 5 ] + rp[ -5 ]) +
				flt[ 3 ] * ( rp[ 7 ] + rp[ -7 ]) +
				flt[ 4 ] * ( rp[ 9 ] + rp[ -9 ]);
	R8BHBC2

	R8BHBC1( convolve6 )
			*op = 0.5 * rp[ 0 ] +
				flt[ 0 ] * ( rp[ 1 ] + rp[ -1 ]) +
				flt[ 1 ] * ( rp[ 3 ] + rp[ -3 ]) +
				flt[ 2 ] * ( rp[ 5 ] + rp[ -5 ]) +
				flt[ 3 ] * ( rp[ 7 ] + rp[ -7 ]) +
				flt[ 4 ] * ( rp[ 9 ] + rp[ -9 ]) +
				flt[ 5 ] * ( rp[ 11 ] + rp[ -11 ]);
	R8BHBC2

	R8BHBC1( convolve7 )
			*op = 0.5 * rp[ 0 ] +
				flt[ 0 ] * ( rp[ 1 ] + rp[ -1 ]) +
				flt[ 1 ] * ( rp[ 3 ] + rp[ -3 ]) +
				flt[ 2 ] * ( rp[ 5 ] + rp[ -5 ]) +
				flt[ 3 ] * ( rp[ 7 ] + rp[ -7 ]) +
				flt[ 4 ] * ( rp[ 9 ] + rp[ -9 ]) +
				flt[ 5 ] * ( rp[ 11 ] + rp[ -11 ]) +
				flt[ 6 ] * ( rp[ 13 ] + rp[ -13 ]);
	R8BHBC2

	R8BHBC1( convolve8 )
			*op = 0.5 * rp[ 0 ] +
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
			*op = 0.5 * rp[ 0 ] +
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
			*op = 0.5 * rp[ 0 ] +
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
			*op = 0.5 * rp[ 0 ] +
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
			*op = 0.5 * rp[ 0 ] +
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
			*op = 0.5 * rp[ 0 ] +
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
			*op = 0.5 * rp[ 0 ] +
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

#endif // R8B_CDSPHBCONVOLVER_INCLUDED
