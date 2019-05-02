//$ nocpp

/**
 * @file CDSPHBUpsampler.h
 *
 * @brief Half-band upsampling class.
 *
 * This file includes half-band upsampling class.
 *
 * r8brain-free-src Copyright (c) 2019 Aleksey Vaneev
 * See the "License.txt" file for license.
 */

#ifndef R8B_CDSPHBUPSAMPLER_INCLUDED
#define R8B_CDSPHBUPSAMPLER_INCLUDED

#include "CDSPProcessor.h"

namespace r8b {

/**
 * @brief Half-band upsampling class.
 *
 * Class implements brute-force half-band 2X upsampling that uses small
 * sparse symmetric FIR filters. It is very efficient and should be used at
 * latter upsampling steps after initial steep 2X upsampling.
 */

class CDSPHBUpsampler : public CDSPProcessor
{
public:
	/**
	 * Constructor initalizes the convolver.
	 *
	 * @param ReqAtten Required half-band filter attentuation.
	 * @param SteepIndex Steepness index - 0=steepest. Corresponds to general
	 * upsampling ratio, e.g. at 4x upsampling 0 is used, at 8x upsampling 1
	 * is used, etc.
	 */

	CDSPHBUpsampler( const double ReqAtten, const int SteepIndex )
	{
		static const int FltCount = 11;
		static const double HBKernel_4[ 4 ] = { // StopAtten = -64.9241 dB
			6.1073830069265700e-001, -1.4463982571471834e-001,
			4.1136036923117736e-002, -7.4740105856911931e-003 };
		static const double HBKernel_5[ 5 ] = { // StopAtten = -87.4775 dB
			6.1553054142504182e-001, -1.5591352339397746e-001,
			5.2404802661294436e-002, -1.4230574348723879e-002,
			2.2457593805369491e-003 };
		static const double HBKernel_6[ 6 ] = { // StopAtten = -104.5154 dB
			6.1883766561934372e-001, -1.6396282655874983e-001,
			6.1104571279329445e-002, -2.0317756445220159e-002,
			5.0264527466028497e-003, -6.9392938429513329e-004 };
		static const double HBKernel_7[ 7 ] = { // StopAtten = -120.6199 dB
			6.2125313688726314e-001, -1.6999763849269822e-001,
			6.8014108060695744e-002, -2.5679821316663798e-002,
			7.9798828249511011e-003, -1.7871060154426299e-003,
			2.1836606459402219e-004 };
		static const double HBKernel_8[ 8 ] = { // StopAtten = -136.5151 dB
			6.2309299085364334e-001, -1.7468969193360939e-001,
			7.3628746444883930e-002, -3.0378268549981380e-002,
			1.0908085227611775e-002, -3.1287343330108630e-003,
			6.3632014609105170e-004, -6.9597139144672582e-005 };
		static const double HBKernel_9[ 9 ] = { // StopAtten = -152.3240 dB
			6.2454069594801442e-001, -1.7844303649906834e-001,
			7.8279410808936412e-002, -3.4501119561946458e-002,
			1.3717889826688978e-002, -4.6090109007730362e-003,
			1.2192752061257472e-003, -2.2647618540871939e-004,
			2.2395554540199925e-005 };
		static const double HBKernel_10[ 10 ] = { // StopAtten = -168.0859 dB
			6.2570883988448456e-001, -1.8151274643051463e-001,
			8.2191863294133416e-002, -3.8131779329263038e-002,
			1.6367492549395291e-002, -6.1530178831008531e-003,
			1.9277693941690544e-003, -4.7165916428615340e-004,
			8.0491894740575276e-005, -7.2581515822012221e-006 };
		static const double HBKernel_11[ 11 ] = { // StopAtten = -183.7962 dB
			6.2667167705228000e-001, -1.8407153338070459e-001,
			8.5529995551434368e-002, -4.1346831405384199e-002,
			1.8844831647669924e-002, -7.7125170092675610e-003,
			2.7268674722185082e-003, -7.9745027943944891e-004,
			1.8116344436958087e-004, -2.8569149403079308e-005,
			2.3667021631368357e-006 };
		static const double HBKernel_12[ 12 ] = { // StopAtten = -199.4768 dB
			6.2747849716660786e-001, -1.8623616747750149e-001,
			8.8409755409981194e-002, -4.4207468337011502e-002,
			2.1149175555773780e-002, -9.2551505772804354e-003,
			3.5871560727791518e-003, -1.1923166995018164e-003,
			3.2627809786633932e-004, -6.9106895864035067e-005,
			1.0122896608955045e-005, -7.7531866615743184e-007 };
		static const double HBKernel_13[ 13 ] = { // StopAtten = -215.1362 dB
			6.2816416241409012e-001, -1.8809076924909807e-001,
			9.0918539438652246e-002, -4.6765502217762434e-002,
			2.3287520076995527e-002, -1.0760626918534744e-002,
			4.4853920783696211e-003, -1.6438774206595230e-003,
			5.1441306611588189e-004, -1.3211723499129135e-004,
			2.6191313459378307e-005, -3.5802417774899453e-006,
			2.5491264654675660e-007 };
		static const double HBKernel_14[ 14 ] = { // StopAtten = -230.7513 dB
			6.2875473134883197e-001, -1.8969941971554194e-001,
			9.3126094863805076e-002, -4.9067251432719283e-002,
			2.5273008898692328e-002, -1.2218646028834490e-002,
			5.4048940066335760e-003, -2.1409917655095256e-003,
			7.4250279810872399e-004, -2.1924535312756177e-004,
			5.3015780896448206e-005, -9.8742950571661936e-006,
			1.2650374272560327e-006, -8.4146504875093342e-008 };
		static const double FltAttens[ FltCount ] = {
			64.9241, 87.4775, 104.5154, 120.6199, 136.5151, 152.3240,
			168.0859, 183.7962, 199.4768, 215.1362, 230.7513 };
		static const double* const FltPtrs[ FltCount ] = { HBKernel_4,
			HBKernel_5, HBKernel_6, HBKernel_7, HBKernel_8, HBKernel_9,
			HBKernel_10, HBKernel_11, HBKernel_12, HBKernel_13,
			HBKernel_14 };
		static const CConvolveFn FltConvFn[ FltCount ] = {
			&CDSPHBUpsampler :: convolve4, &CDSPHBUpsampler :: convolve5,
			&CDSPHBUpsampler :: convolve6, &CDSPHBUpsampler :: convolve7,
			&CDSPHBUpsampler :: convolve8, &CDSPHBUpsampler :: convolve9,
			&CDSPHBUpsampler :: convolve10, &CDSPHBUpsampler :: convolve11,
			&CDSPHBUpsampler :: convolve12, &CDSPHBUpsampler :: convolve13,
			&CDSPHBUpsampler :: convolve14 };

		static const int FltCountB = 7;
		static const double HBKernel_2b[ 2 ] = { // StopAtten = -46.2556 dB
			5.6965643437574798e-001, -7.0243561190601822e-002 };
		static const double HBKernel_3b[ 3 ] = { // StopAtten = -93.6536 dB
			5.9040785316467748e-001, -1.0462338733801557e-001,
			1.4234900265846395e-002 };
		static const double HBKernel_4b[ 4 ] = { // StopAtten = -123.4514 dB
			6.0140277542879073e-001, -1.2564483854573050e-001,
			2.7446500598030887e-002, -3.2051079559036744e-003 };
		static const double HBKernel_5b[ 5 ] = { // StopAtten = -152.4403 dB
			6.0818642429044178e-001, -1.3981140187082763e-001,
			3.8489164053787661e-002, -7.6218861795043225e-003,
			7.5772358126258155e-004 };
		static const double HBKernel_6b[ 6 ] = { // StopAtten = -181.2501 dB
			6.1278392271870352e-001, -1.5000053763409249e-001,
			4.7575323519283508e-002, -1.2320702806281281e-002,
			2.1462442604041065e-003, -1.8425092396978648e-004 };
		static const double HBKernel_7b[ 7 ] = { // StopAtten = -209.9472 dB
			6.1610372237019151e-001, -1.5767891821295410e-001,
			5.5089690570484962e-002, -1.6895755290596615e-002,
			3.9416641999499014e-003, -6.0603620400878633e-004,
			4.5632598748568398e-005 };
		static const double HBKernel_8b[ 8 ] = { // StopAtten = -238.5612 dB
			6.1861282849648180e-001, -1.6367179296640288e-001,
			6.1369859727781417e-002, -2.1184465440918565e-002,
			5.9623352367661475e-003, -1.2483096884685629e-003,
			1.7099294398059683e-004, -1.1448310399897466e-005 };
		static const double FltAttensB[ FltCountB ] = {
			46.2556, 93.6536, 123.4514, 152.4403, 181.2501, 209.9472,
			238.5612 };
		static const double* const FltPtrsB[ FltCountB ] = { HBKernel_2b,
			HBKernel_3b, HBKernel_4b, HBKernel_5b, HBKernel_6b, HBKernel_7b,
			HBKernel_8b };
		static const CConvolveFn FltConvFnB[ FltCountB ] = {
			&CDSPHBUpsampler :: convolve2, &CDSPHBUpsampler :: convolve3,
			&CDSPHBUpsampler :: convolve4, &CDSPHBUpsampler :: convolve5,
			&CDSPHBUpsampler :: convolve6, &CDSPHBUpsampler :: convolve7,
			&CDSPHBUpsampler :: convolve8 };

		static const int FltCountC = 5;
		static const double HBKernel_2c[ 2 ] = { // StopAtten = -87.9438 dB
			5.6430278013478086e-001, -6.4338068855764208e-002 };
		static const double HBKernel_3c[ 3 ] = { // StopAtten = -130.8862 dB
			5.8706402915553113e-001, -9.9362380958695873e-002,
			1.2298637065878193e-002 };
		static const double HBKernel_4c[ 4 ] = { // StopAtten = -172.3191 dB
			5.9896586135108265e-001, -1.2111680603660968e-001,
			2.4763118077755664e-002, -2.6121758134936002e-003 };
		static const double HBKernel_5c[ 5 ] = { // StopAtten = -213.4984 dB
			6.0626808278478261e-001, -1.3588224019070938e-001,
			3.5544305138258458e-002, -6.5127022013993230e-003,
			5.8255449020627736e-004 };
		static const double HBKernel_6c[ 6 ] = { // StopAtten = -254.5179 dB
			6.1120171157732273e-001, -1.4654486624691154e-001,
			4.4582957343679119e-002, -1.0840542911916273e-002,
			1.7343703931622656e-003, -1.3363015552414481e-004 };
		static const double FltAttensC[ FltCountC ] = {
			87.9438, 130.8862, 172.3191, 213.4984, 254.5179 };
		static const double* const FltPtrsC[ FltCountC ] = { HBKernel_2c,
			HBKernel_3c, HBKernel_4c, HBKernel_5c, HBKernel_6c };
		static const CConvolveFn FltConvFnC[ FltCountC ] = {
			&CDSPHBUpsampler :: convolve2, &CDSPHBUpsampler :: convolve3,
			&CDSPHBUpsampler :: convolve4, &CDSPHBUpsampler :: convolve5,
			&CDSPHBUpsampler :: convolve6 };

		int k = 0;
		int fltt;
		double att;

		if( SteepIndex == 0 )
		{
			while( k != FltCount - 1 && FltAttens[ k ] < ReqAtten )
			{
				k++;
			}

			flt = FltPtrs[ k ];
			convfn = FltConvFn[ k ];
			fltt = 4 + k;
			att = FltAttens[ k ];
		}
		else
		if( SteepIndex == 1 )
		{
			while( k != FltCountB - 1 && FltAttensB[ k ] < ReqAtten )
			{
				k++;
			}

			flt = FltPtrsB[ k ];
			convfn = FltConvFnB[ k ];
			fltt = 2 + k;
			att = FltAttensB[ k ];
		}
		else
		{
			while( k != FltCountC - 1 && FltAttensC[ k ] < ReqAtten )
			{
				k++;
			}

			flt = FltPtrsC[ k ];
			convfn = FltConvFnC[ k ];
			fltt = 2 + k;
			att = FltAttensC[ k ];
		}

		fll = fltt - 1;
		fl2 = fltt;
		fl4 = fll + fl2;

		R8BCONSOLE( "CDSPHBUpsampler: sti=%i taps=%i att=%.1f io=2/1\n",
			SteepIndex, fltt, att );

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

		return( MaxInLen * 2 );
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
				const int c = BufLeft - fl2;

				( *this.*convfn )( op, c );

				op += c + c;
				BufLeft -= c;
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
	double Buf[ BufLen + 27 ]; ///< The ring buffer, including overrun
		///< protection for the largest filter.
		///<
	const double* flt; ///< Half-band filter taps.
		///<
	int fll; ///< Input latency.
		///<
	int fl2; ///< Right-side filter length.
		///<
	int fl4; ///< Overrrun length.
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
	typedef void( CDSPHBUpsampler :: *CConvolveFn )( double* op, int c ); ///<
		///< Convolution funtion type.
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
			op[ 0 ] = rp[ 0 ];
#define R8BHBC2 \
			rpos = ( rpos + 1 ) & BufLenMask; \
			op += 2; \
			c--; \
		} \
		ReadPos = rpos; \
	}

	R8BHBC1( convolve2 )
			op[ 1 ] =
				flt[ 0 ] * ( rp[ 1 ] + rp[ 0 ]) +
				flt[ 1 ] * ( rp[ 2 ] + rp[ -1 ]);
	R8BHBC2

	R8BHBC1( convolve3 )
			op[ 1 ] =
				flt[ 0 ] * ( rp[ 1 ] + rp[ 0 ]) +
				flt[ 1 ] * ( rp[ 2 ] + rp[ -1 ]) +
				flt[ 2 ] * ( rp[ 3 ] + rp[ -2 ]);
	R8BHBC2

	R8BHBC1( convolve4 )
			op[ 1 ] =
				flt[ 0 ] * ( rp[ 1 ] + rp[ 0 ]) +
				flt[ 1 ] * ( rp[ 2 ] + rp[ -1 ]) +
				flt[ 2 ] * ( rp[ 3 ] + rp[ -2 ]) +
				flt[ 3 ] * ( rp[ 4 ] + rp[ -3 ]);
	R8BHBC2

	R8BHBC1( convolve5 )
			op[ 1 ] =
				flt[ 0 ] * ( rp[ 1 ] + rp[ 0 ]) +
				flt[ 1 ] * ( rp[ 2 ] + rp[ -1 ]) +
				flt[ 2 ] * ( rp[ 3 ] + rp[ -2 ]) +
				flt[ 3 ] * ( rp[ 4 ] + rp[ -3 ]) +
				flt[ 4 ] * ( rp[ 5 ] + rp[ -4 ]);
	R8BHBC2

	R8BHBC1( convolve6 )
			op[ 1 ] =
				flt[ 0 ] * ( rp[ 1 ] + rp[ 0 ]) +
				flt[ 1 ] * ( rp[ 2 ] + rp[ -1 ]) +
				flt[ 2 ] * ( rp[ 3 ] + rp[ -2 ]) +
				flt[ 3 ] * ( rp[ 4 ] + rp[ -3 ]) +
				flt[ 4 ] * ( rp[ 5 ] + rp[ -4 ]) +
				flt[ 5 ] * ( rp[ 6 ] + rp[ -5 ]);
	R8BHBC2

	R8BHBC1( convolve7 )
			op[ 1 ] =
				flt[ 0 ] * ( rp[ 1 ] + rp[ 0 ]) +
				flt[ 1 ] * ( rp[ 2 ] + rp[ -1 ]) +
				flt[ 2 ] * ( rp[ 3 ] + rp[ -2 ]) +
				flt[ 3 ] * ( rp[ 4 ] + rp[ -3 ]) +
				flt[ 4 ] * ( rp[ 5 ] + rp[ -4 ]) +
				flt[ 5 ] * ( rp[ 6 ] + rp[ -5 ]) +
				flt[ 6 ] * ( rp[ 7 ] + rp[ -6 ]);
	R8BHBC2

	R8BHBC1( convolve8 )
			op[ 1 ] =
				flt[ 0 ] * ( rp[ 1 ] + rp[ 0 ]) +
				flt[ 1 ] * ( rp[ 2 ] + rp[ -1 ]) +
				flt[ 2 ] * ( rp[ 3 ] + rp[ -2 ]) +
				flt[ 3 ] * ( rp[ 4 ] + rp[ -3 ]) +
				flt[ 4 ] * ( rp[ 5 ] + rp[ -4 ]) +
				flt[ 5 ] * ( rp[ 6 ] + rp[ -5 ]) +
				flt[ 6 ] * ( rp[ 7 ] + rp[ -6 ]) +
				flt[ 7 ] * ( rp[ 8 ] + rp[ -7 ]);
	R8BHBC2

	R8BHBC1( convolve9 )
			op[ 1 ] =
				flt[ 0 ] * ( rp[ 1 ] + rp[ 0 ]) +
				flt[ 1 ] * ( rp[ 2 ] + rp[ -1 ]) +
				flt[ 2 ] * ( rp[ 3 ] + rp[ -2 ]) +
				flt[ 3 ] * ( rp[ 4 ] + rp[ -3 ]) +
				flt[ 4 ] * ( rp[ 5 ] + rp[ -4 ]) +
				flt[ 5 ] * ( rp[ 6 ] + rp[ -5 ]) +
				flt[ 6 ] * ( rp[ 7 ] + rp[ -6 ]) +
				flt[ 7 ] * ( rp[ 8 ] + rp[ -7 ]) +
				flt[ 8 ] * ( rp[ 9 ] + rp[ -8 ]);
	R8BHBC2

	R8BHBC1( convolve10 )
			op[ 1 ] =
				flt[ 0 ] * ( rp[ 1 ] + rp[ 0 ]) +
				flt[ 1 ] * ( rp[ 2 ] + rp[ -1 ]) +
				flt[ 2 ] * ( rp[ 3 ] + rp[ -2 ]) +
				flt[ 3 ] * ( rp[ 4 ] + rp[ -3 ]) +
				flt[ 4 ] * ( rp[ 5 ] + rp[ -4 ]) +
				flt[ 5 ] * ( rp[ 6 ] + rp[ -5 ]) +
				flt[ 6 ] * ( rp[ 7 ] + rp[ -6 ]) +
				flt[ 7 ] * ( rp[ 8 ] + rp[ -7 ]) +
				flt[ 8 ] * ( rp[ 9 ] + rp[ -8 ]) +
				flt[ 9 ] * ( rp[ 10 ] + rp[ -9 ]);
	R8BHBC2

	R8BHBC1( convolve11 )
			op[ 1 ] =
				flt[ 0 ] * ( rp[ 1 ] + rp[ 0 ]) +
				flt[ 1 ] * ( rp[ 2 ] + rp[ -1 ]) +
				flt[ 2 ] * ( rp[ 3 ] + rp[ -2 ]) +
				flt[ 3 ] * ( rp[ 4 ] + rp[ -3 ]) +
				flt[ 4 ] * ( rp[ 5 ] + rp[ -4 ]) +
				flt[ 5 ] * ( rp[ 6 ] + rp[ -5 ]) +
				flt[ 6 ] * ( rp[ 7 ] + rp[ -6 ]) +
				flt[ 7 ] * ( rp[ 8 ] + rp[ -7 ]) +
				flt[ 8 ] * ( rp[ 9 ] + rp[ -8 ]) +
				flt[ 9 ] * ( rp[ 10 ] + rp[ -9 ]) +
				flt[ 10 ] * ( rp[ 11 ] + rp[ -10 ]);
	R8BHBC2

	R8BHBC1( convolve12 )
			op[ 1 ] =
				flt[ 0 ] * ( rp[ 1 ] + rp[ 0 ]) +
				flt[ 1 ] * ( rp[ 2 ] + rp[ -1 ]) +
				flt[ 2 ] * ( rp[ 3 ] + rp[ -2 ]) +
				flt[ 3 ] * ( rp[ 4 ] + rp[ -3 ]) +
				flt[ 4 ] * ( rp[ 5 ] + rp[ -4 ]) +
				flt[ 5 ] * ( rp[ 6 ] + rp[ -5 ]) +
				flt[ 6 ] * ( rp[ 7 ] + rp[ -6 ]) +
				flt[ 7 ] * ( rp[ 8 ] + rp[ -7 ]) +
				flt[ 8 ] * ( rp[ 9 ] + rp[ -8 ]) +
				flt[ 9 ] * ( rp[ 10 ] + rp[ -9 ]) +
				flt[ 10 ] * ( rp[ 11 ] + rp[ -10 ]) +
				flt[ 11 ] * ( rp[ 12 ] + rp[ -11 ]);
	R8BHBC2

	R8BHBC1( convolve13 )
			op[ 1 ] =
				flt[ 0 ] * ( rp[ 1 ] + rp[ 0 ]) +
				flt[ 1 ] * ( rp[ 2 ] + rp[ -1 ]) +
				flt[ 2 ] * ( rp[ 3 ] + rp[ -2 ]) +
				flt[ 3 ] * ( rp[ 4 ] + rp[ -3 ]) +
				flt[ 4 ] * ( rp[ 5 ] + rp[ -4 ]) +
				flt[ 5 ] * ( rp[ 6 ] + rp[ -5 ]) +
				flt[ 6 ] * ( rp[ 7 ] + rp[ -6 ]) +
				flt[ 7 ] * ( rp[ 8 ] + rp[ -7 ]) +
				flt[ 8 ] * ( rp[ 9 ] + rp[ -8 ]) +
				flt[ 9 ] * ( rp[ 10 ] + rp[ -9 ]) +
				flt[ 10 ] * ( rp[ 11 ] + rp[ -10 ]) +
				flt[ 11 ] * ( rp[ 12 ] + rp[ -11 ]) +
				flt[ 12 ] * ( rp[ 13 ] + rp[ -12 ]);
	R8BHBC2

	R8BHBC1( convolve14 )
			op[ 1 ] =
				flt[ 0 ] * ( rp[ 1 ] + rp[ 0 ]) +
				flt[ 1 ] * ( rp[ 2 ] + rp[ -1 ]) +
				flt[ 2 ] * ( rp[ 3 ] + rp[ -2 ]) +
				flt[ 3 ] * ( rp[ 4 ] + rp[ -3 ]) +
				flt[ 4 ] * ( rp[ 5 ] + rp[ -4 ]) +
				flt[ 5 ] * ( rp[ 6 ] + rp[ -5 ]) +
				flt[ 6 ] * ( rp[ 7 ] + rp[ -6 ]) +
				flt[ 7 ] * ( rp[ 8 ] + rp[ -7 ]) +
				flt[ 8 ] * ( rp[ 9 ] + rp[ -8 ]) +
				flt[ 9 ] * ( rp[ 10 ] + rp[ -9 ]) +
				flt[ 10 ] * ( rp[ 11 ] + rp[ -10 ]) +
				flt[ 11 ] * ( rp[ 12 ] + rp[ -11 ]) +
				flt[ 12 ] * ( rp[ 13 ] + rp[ -12 ]) +
				flt[ 13 ] * ( rp[ 14 ] + rp[ -13 ]);
	R8BHBC2

#undef R8BHBC1
#undef R8BHBC2
};

// ---------------------------------------------------------------------------

} // namespace r8b

#endif // R8B_CDSPHBUPSAMPLER_INCLUDED
