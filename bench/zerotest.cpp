#include "/libvox/Sources/Audio/AudioMath.h"
#include "/libvox/Sources/Core/AppMain.h"
#include "/libvox/Sources/Other/CWaveFile.h"
#include "../CDSPResampler.h"

/**
 * @file masstest.cpp
 *
 * @brief Mass zeroing test - checks for sample rate precision by upsampling
 * and downsampling a reference signal.
 *
 * r8brain-free-src Copyright (c) 2019 Aleksey Vaneev
 * See the "License.txt" file for license.
 */

double calcRMS( const double* const p1, const double* const p2, const int l )
{
	double s = 0.0;
	int i;

	for( i = 0; i < l; i++ )
	{
		const double d = p1[ i ] - p2[ i ];
		s += d * d;
	}

	return( sqrt( s / l ));
}

VOXMAIN
{
	Args.setProgramDescription( "This utility program perform mass zeroing "
		"test of sample rate converter." );

	Args.addReqArg( argtFileName, "in-file", "Input WAV filename." );

	VOXCHECK( Args.checkArgs() );

	CWaveFile inf;
	VOXCHECK( inf.loadFile( *Args.getArgValue( "in-file" ).ValueStr ));

	const int InBufSize = (int) min( 50000, inf.SampleCount );

	CInitArray< CFixedBuffer< double > > InBufs( inf.ChannelCount );
	int i;

	for( i = 0; i < inf.ChannelCount; i++ )
	{
		InBufs[ i ].alloc( InBufSize );
	}

	int ReadCount;
	VOXCHECK( inf.readData( InBufs, InBufSize, ReadCount ));

	const int MaxInLen = 521;
	const double tb = 2.0;

	// Create reference signal which has 9/10 bandwidth of the input signal.

	CFixedBuffer< double > Ref( InBufSize );

	if( true )
	{
		const int Ref0Size = (int) ( InBufSize * 9.0 / 10.0 );
		CFixedBuffer< double > Ref0( Ref0Size );

		CPtrKeeper< r8b :: CDSPResampler24* > Resamp1;
		Resamp1 = new r8b :: CDSPResampler24( 10.0, 9.0, MaxInLen, tb );

		CPtrKeeper< r8b :: CDSPResampler24* > Resamp2;
		Resamp2 = new r8b :: CDSPResampler24( 9.0, 10.0, MaxInLen, tb );

		Resamp1 -> oneshot( MaxInLen, &InBufs[ 0 ][ 0 ], InBufSize,
			&Ref0[ 0 ], Ref0Size );

		Resamp2 -> oneshot( MaxInLen, &Ref0[ 0 ], Ref0Size,
			&Ref[ 0 ], InBufSize );
	}

	const double SrcSampleRate = 10.0;
	int k;

	for( k = 11; k < 200; k++ )
	{
		const double DstSampleRate = k;
		const int ol1 = (int) ( InBufSize * DstSampleRate / SrcSampleRate );
		CFixedBuffer< double > OutBuf1( ol1 );
		CFixedBuffer< double > OutBuf2( InBufSize );

		CPtrKeeper< r8b :: CDSPResampler24* > Resamp1;
		Resamp1 = new r8b :: CDSPResampler24( SrcSampleRate,
			DstSampleRate, MaxInLen, tb );

		CPtrKeeper< r8b :: CDSPResampler24* > Resamp2;
		Resamp2 = new r8b :: CDSPResampler24( DstSampleRate,
			SrcSampleRate, MaxInLen, tb );

		Resamp1 -> oneshot( MaxInLen, &Ref[ 0 ], InBufSize,
			&OutBuf1[ 0 ], ol1 );

		Resamp2 -> oneshot( MaxInLen, &OutBuf1[ 0 ], ol1,
			&OutBuf2[ 0 ], InBufSize );

		const double r = calcRMS( &Ref[ 5000 ], &OutBuf2[ 5000 ],
			InBufSize - 10000 );

		printf( "%3.0f/%3.0f: ", DstSampleRate, SrcSampleRate );

		if( r == 0 )
		{
			printf( "-inf\n" );
		}
		else
		{
			printf( "%.2f\n", 20.0 * log( r ) / log( 10.0 ));
		}
	}

	VOXRET;
}
