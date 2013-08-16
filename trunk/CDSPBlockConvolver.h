//$ nocpp

/**
 * \file CDSPBlockConvolver.h
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
 * \brief Single-block overlap-save convolution processing class.
 *
 * Class that implements single-block overlap-save convolution processing. The
 * size of the single FFT block used depends on the length of the filter
 * kernel.
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
	 */

	CDSPBlockConvolver( CDSPFIRFilter& aFilter )
		: Filter( &aFilter )
		, ffto( Filter -> BlockSizeBits + 1 )
		, BlockSize( 1 << Filter -> BlockSizeBits )
		, Latency( BlockSize + Filter -> Latency )
	{
		const int bs = BlockSize * (int) sizeof( double );
		PrevInput.alloc( bs );
		WorkBlocks[ 0 ].alloc( bs * 2 );
		WorkBlocks[ 1 ].alloc( bs * 2 );
		CurInput = WorkBlocks[ 0 ];
		CurOutput = WorkBlocks[ 1 ];

		memset( &PrevInput[ 0 ], 0, bs );
		InDataLeft = BlockSize;
		LatencyLeft = Latency;
	}

	~CDSPBlockConvolver()
	{
		Filter -> unref();
	}

	/**
	 * Function performs convolution processing.
	 *
	 * @param ip Input data pointer.
	 * @param[out] op0 Output data pointer (can be equal to "ip"). On
	 * function's exit this variable will point to the actual output data
	 * within the output buffer.
	 * @param l0 How many samples to process.
	 * @return The number of output samples available after processing. This
	 * value can be smaller in comparison to the original "l0" value due to
	 * processing and filter's latency compensation that took place.
	 */

	int process( const double* ip, double*& op0, const int l0 )
	{
		R8BASSERT( l0 >= 0 );
		R8BASSERT( ip != NULL || l0 == 0 );
		R8BASSERT( op0 != NULL || l0 == 0 );

		double* op = op0;
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
			ffto -> multiplyBlocks( Filter -> KernelBlock, CurInput );
			ffto -> inverse( CurInput );

			double* const tmp = CurInput;
			CurInput = CurOutput;
			CurOutput = tmp;

			ip += b;
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
			const int res = l0 - LatencyLeft;
			LatencyLeft = 0;
			return( res );
		}

		return( l0 );
	}

private:
	CDSPFIRFilter* Filter; ///< Filter in use.
	CDSPRealFFTKeeper ffto; ///< FFT object, size Filter.BlockSizeBits + 1.
	int BlockSize; ///< Block size in samples.
	int Latency; ///< Processing and kernel's latency.
	CFixedBuffer< double > PrevInput; ///< Previous input data,
		///< size = BlockSize.
		///<
	CFixedBuffer< double > WorkBlocks[ 2 ]; ///< Input and output blocks,
		///< size = BlockSize * 2 each. Used in flip-flop manner.
		///<
	double* CurInput; ///< Input data buffer.
	double* CurOutput; ///< Output data buffer.
	int InDataLeft; ///< Samples left before processing input and output FFT
		///< blocks.
		///<
	int LatencyLeft; ///< Latency in samples left to skip.

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

		for( i = 0; i < l; i++ )
		{
			opoffs[ i ] = prev[ i ];
			prev[ i ] = ip[ i ];
			op[ i ] = ip[ i ];
		}
	}
};

} // namespace r8b

#endif // R8B_CDSPBLOCKCONVOLVER_INCLUDED
