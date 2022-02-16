#include "../../../libvox/Sources/Audio/AudioMath.h"
#include "../../../libvox/Sources/Core/AppMain.h"
#include "../../../libvox/Sources/Other/CWaveFile.h"
#include "../CDSPResampler.h"

typedef r8b :: CDSPResampler24 CResamp;

/**
 * @file zerotest.cpp
 *
 * @brief Mass zeroing test - checks for sample rate conversion precision by
 * upsampling and downsampling a reference signal.
 *
 * r8brain-free-src Copyright (c) 2013-2021 Aleksey Vaneev
 * See the "LICENSE" file for license.
 */

double calcRMS( const double* const p1, const double* const p2, const int l,
	double& peakd )
{
	double s = 0.0;
	int i;

	for( i = 0; i < l; i++ )
	{
		const double d = fabs( p1[ i ] - p2[ i ]);
		peakd = max( peakd, d );
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

	const int InBufSize = (int) min( (int64_t) 50000, inf.SampleCount );

	CInitArray< CFixedBuffer< double > > InBufs( inf.ChannelCount );
	int i;

	for( i = 0; i < inf.ChannelCount; i++ )
	{
		InBufs[ i ].alloc( InBufSize );
	}

	int ReadCount;
	VOXCHECK( inf.readData( InBufs, InBufSize, ReadCount ));

	CRnd rnd;
	rnd.init( 0 );

	// Create reference signal which has 92% bandwidth of the input signal.

	CFixedBuffer< double > Ref( InBufSize );

	if( true )
	{
		const int MaxInLen = 521;
		const double tb = 2.0;
		const double bw = 9.3;
		const int Ref0Size = (int) ( InBufSize * bw / 10.0 );
		CFixedBuffer< double > Ref0( Ref0Size );

		CPtrKeeper< r8b :: CDSPResampler24* > Resamp1;
		Resamp1 = new r8b :: CDSPResampler24( 10.0, bw, MaxInLen, tb );

		CPtrKeeper< r8b :: CDSPResampler24* > Resamp2;
		Resamp2 = new r8b :: CDSPResampler24( bw, 10.0, MaxInLen, tb );

		Resamp1 -> oneshot( &InBufs[ 0 ][ 0 ], InBufSize, &Ref0[ 0 ],
			Ref0Size );

		Resamp2 -> oneshot( &Ref0[ 0 ], Ref0Size, &Ref[ 0 ], InBufSize );
	}

	const double SrcSampleRate = 20.0;
	double avgr = 0.0;
	double peakd = 0.0;
	int avgc = 0;
	double avgperf = 0.0;
	double avglatency = 0.0;
	int k;

	for( k = 21; k <= 640; k++ )
	{
		const int MaxInLen = (int) ( 50 + rnd.getUniform() * 1500 );
		const double tb = 0.5 + rnd.getUniform() * 4.5;

		const double DstSampleRate = k;
		const int ol1 = (int) ( InBufSize * DstSampleRate / SrcSampleRate );
		CFixedBuffer< double > OutBuf1( ol1 );
		CFixedBuffer< double > OutBuf2( InBufSize );

		CPtrKeeper< CResamp* > Resamp1;
		Resamp1 = new CResamp( SrcSampleRate, DstSampleRate, MaxInLen, tb );
		avglatency += Resamp1 -> getInLenBeforeOutStart();

		CPtrKeeper< CResamp* > Resamp2;
		Resamp2 = new CResamp( DstSampleRate, SrcSampleRate, MaxInLen, tb );

		const TClock t1( CSystem :: getClock() );
		Resamp1 -> oneshot( &Ref[ 0 ], InBufSize, &OutBuf1[ 0 ], ol1 );
		double perf = 1e-6 * InBufSize /
			CSystem :: getClockDiffSec( t1 );

		const TClock t2( CSystem :: getClock() );
		Resamp2 -> oneshot( &OutBuf1[ 0 ], ol1, &OutBuf2[ 0 ], InBufSize );
		perf = ( perf + 1e-6 * InBufSize /
			CSystem :: getClockDiffSec( t2 )) * 0.5;

		const double r = calcRMS( &Ref[ 5000 ], &OutBuf2[ 5000 ],
			InBufSize - 10000, peakd );

		avgr += r * r;
		avgperf += perf * perf;
		avgc++;

		printf( "%3.0f/%3.0f\t", DstSampleRate, SrcSampleRate );

		if( r == 0 )
		{
			printf( "-inf");
		}
		else
		{
			printf( "%7.2f", 20.0 * log( r ) / log( 10.0 ));
		}

		printf( "\t%.2f\tMrops\n", perf );
	}

	printf( "Average rms %.2f\n", 10.0 * log( avgr / avgc ) / log( 10.0 ));
	printf( "Peak diff %.2f\n", 20.0 * log( peakd ) / log( 10.0 ));
	printf( "Average perf %.2f Mrops\n", sqrt( avgperf / avgc ));
	printf( "Average latency %.0f\n", avglatency / avgc );

	VOXRET;
}
