//$ nocpp

/**
 * \file CDSPBlockConvolver.h
 *
 * \brief Single-block overlap-save convolution processor class.
 *
 * This file includes single-block overlap-save convolution processor class.
 *
 * r8brain-free-src Copyright (c) 2013 Aleksey Vaneev
 * See the "License.txt" file for license.
 */

#ifndef R8B_CDSPBLOCKCONVOLVER_INCLUDED
#define R8B_CDSPBLOCKCONVOLVER_INCLUDED

#include "CDSPFIRFilter.h"

namespace r8b {

/**
 * Enumeration of the built-in resampling modes.
 */

enum CDSPResamplingMode
{
	rsNone = 0, ///< No resampling.
	rsUpsample2X, ///< 2x upsampling.
	rsDownsample2X ///< 2x downsampling.
};

/**
 * \brief Single-block overlap-save convolution processing class.
 *
 * Class that implements single-block overlap-save convolution processing. The
 * size of the single FFT block used depends on the length of the filter
 * kernel.
 *
 * The rationale behind "single-block" processing is that increasing the FFT
 * block size by 2 is more efficient than performing convolution at the same
 * FFT block size using two blocks.
 *
 * This class also implements a built-in resampling (2X up or 2X down) which
 * simplifies the overall resampling objects topology.
 */

class CDSPBlockConvolver : public R8B_BASECLASS
{
	R8BNOCTOR( CDSPBlockConvolver );

public:
	/**
	 * Constructor initializes internal variables and constants of *this
	 * object.
	 *
	 * @param aFilter Pre-calculated filter data. Reference to this object is
	 * inhertied by *this object, and the object will be released when *this
	 * object is destroyed.
	 * @param aResamplingMode Resampling mode to use: see
	 * CDSPBlockConvolver::CDSPResamplingMode. Note that the filter should
	 * meet alias-free resampling requirements if resampling is used: the
	 * normalized low-pass frequency of the filter should be equal to 0.5 or
	 * lower.
	 */

	CDSPBlockConvolver( CDSPFIRFilter& aFilter,
		const CDSPResamplingMode aResamplingMode =
		CDSPResamplingMode :: rsNone )
		: Filter( &aFilter )
		, ffto( Filter -> getBlockSizeBits() + 1 )
		, ResamplingMode( aResamplingMode )
		, UpShift( ResamplingMode ==
			CDSPResamplingMode :: rsUpsample2X ? 1 : 0 )
		, BlockSize( 1 << Filter -> getBlockSizeBits() )
		, Latency( BlockSize + Filter -> getLatency() )
	{
		const int bs = BlockSize * (int) sizeof( double );

		PrevInput.alloc( bs );
		memset( &PrevInput[ 0 ], 0, bs );

		WorkBlocks[ 0 ].alloc( bs * 2 );
		CurInput = WorkBlocks[ 0 ];

		WorkBlocks[ 1 ].alloc( bs * 2 );
		CurOutput = WorkBlocks[ 1 ];

		InDataLeft = BlockSize;
		LatencyLeft = Latency;
		DownSkip = 0;
	}

	~CDSPBlockConvolver()
	{
		Filter -> unref();
	}

	/**
	 * @return The number of samples that should be passed to *this object
	 * before the actual output starts. This value includes latencies induced
	 * by both the convolver and filter.
	 */

	int getInLenBeforeOutStart() const
	{
		return( UpShift == 0 ? Latency : ( Latency + 1 ) >> 1 );
	}

	/**
	 * @param MaxInLen The number of samples planned to process at once, at
	 * most.
	 * @return The minimal length of the output buffer required when
	 * processing the "MaxInLen" number of input samples.
	 */

	int getMaxOutLen( const int MaxInLen ) const
	{
		R8BASSERT( MaxInLen >= 0 );

		if( ResamplingMode == CDSPResamplingMode :: rsUpsample2X )
		{
			return( MaxInLen * 2 );
		}

		if( ResamplingMode == CDSPResamplingMode :: rsNone )
		{
			return( MaxInLen );
		}

		return(( MaxInLen + 1 ) >> 1 );
	}

	/**
	 * Function performs convolution processing.
	 *
	 * @param ip Input data pointer.
	 * @param[out] op0 Output data pointer (can be equal to "ip" only if
	 * rsNone or rsDownsample2X resampling mode is used). On function's exit
	 * this variable will point to the actual output data within the output
	 * buffer. The length of the output buffer if resampling mode rsUpsample2X
	 * is used should be twice as large as the input buffer. If rsDownsample2X
	 * is used the size should be equal to "l0 / 2 + 1".
	 * @param l0 How many samples to process.
	 * @return The number of output samples available after processing. This
	 * value can be smaller in comparison to the original "l0" value due to
	 * processing and filter's latency compensation that took place, and due
	 * to resampling if it is used.
	 */

	int process( const double* ip, double*& op0, int l0 )
	{
		R8BASSERT( l0 >= 0 );
		R8BASSERT( UpShift == 0 || ip != op0 || l0 == 0 );

		double* op = op0;
		l0 <<= UpShift;
		int l = l0;

		while( l > 0 )
		{
			const int Offs = BlockSize - InDataLeft;

			if( l < InDataLeft )
			{
				InDataLeft -= l;
				acquireInput( ip, Offs, l );
				memcpy( op, &CurOutput[ Offs ], l * sizeof( double ));

				break;
			}

			const int b = InDataLeft;
			InDataLeft = BlockSize;

			acquireInput( ip, Offs, b );
			memcpy( op, &CurOutput[ Offs ], b * sizeof( double ));

			ffto -> forward( CurInput );
			ffto -> multiplyBlocks( Filter -> getKernelBlock(), CurInput );
			ffto -> inverse( CurInput );

			double* const tmp = CurInput;
			CurInput = CurOutput;
			CurOutput = tmp;

			ip += b >> UpShift;
			op += b;
			l -= b;
		}

		if( LatencyLeft > 0 )
		{
			if( LatencyLeft >= l0 )
			{
				LatencyLeft -= l0;
				return( 0 );
			}

			op0 += LatencyLeft;
			l0 -= LatencyLeft;
			LatencyLeft = 0;
		}

		if( ResamplingMode == CDSPResamplingMode :: rsDownsample2X && l0 > 0 )
		{
			// Perform quick 2x downsampling.

			const double* ip = op0 + DownSkip;
			l = ( l0 + 1 - DownSkip ) >> 1;
			DownSkip += ( l << 1 ) - l0;
			l0 = l;
			op = op0;

			while( l > 0 )
			{
				*op = *ip;
				op++;
				ip += 2;
				l--;
			}
		}

		return( l0 );
	}

private:
	CDSPFIRFilter* Filter; ///< Filter in use.
		///<
	CDSPRealFFTKeeper ffto; ///< FFT object, size Filter.BlockSizeBits + 1.
		///<
	CDSPResamplingMode ResamplingMode; ///< Built-in resampling mode to use.
		///<
	int UpShift; ///< Upsampling shift. Equals 1 if 2x upsampling is used,
		///< 0 otherwise.
		///<
	int BlockSize; ///< Block size in samples.
		///<
	int Latency; ///< Processing and kernel's latency.
		///<
	CFixedBuffer< double > PrevInput; ///< Previous input data,
		///< size = BlockSize.
		///<
	CFixedBuffer< double > WorkBlocks[ 2 ]; ///< Input and output blocks,
		///< size = BlockSize * 2 each. Used in flip-flop manner.
		///<
	double* CurInput; ///< Input data buffer.
		///<
	double* CurOutput; ///< Output data buffer.
		///<
	int InDataLeft; ///< Samples left before processing input and output FFT
		///< blocks.
		///<
	int LatencyLeft; ///< Latency in samples left to skip.
		///<
	int DownSkip; ///< The current 2x downsampling sample skip (0 or 1).
		///<

	CDSPBlockConvolver()
	{
	}

	/**
	 * Function acquires input data and distributes it over the work buffers.
	 *
	 * @param ip Input buffer.
	 * @param Offs Current input/output offset, 0..BlockSize-1.
	 * @param l Sample count.
	 */

	void acquireInput( const double* const ip, const int Offs, const int l )
	{
		double* const prev = &PrevInput[ Offs ];
		double* const op = &CurInput[ Offs ];
		double* const opoffs = op + BlockSize;
		int i;

		if( UpShift == 0 )
		{
			for( i = 0; i < l; i++ )
			{
				opoffs[ i ] = prev[ i ];
				prev[ i ] = ip[ i ];
				op[ i ] = ip[ i ];
			}
		}
		else
		{
			// Insert zero after each input sample to perform quick 2x
			// upsampling. "l" is always an even value in such resampling
			// mode.

			for( i = 0; i < l; i += 2 )
			{
				opoffs[ i ] = prev[ i ];
				opoffs[ i + 1 ] = 0.0;
				const double v = ip[ i >> 1 ] * 2.0;
				prev[ i ] = v;
				op[ i ] = v;
				op[ i + 1 ] = 0.0;
			}
		}
	}
};

} // namespace r8b

#endif // R8B_CDSPBLOCKCONVOLVER_INCLUDED
