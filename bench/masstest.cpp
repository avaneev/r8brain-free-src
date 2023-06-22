#include "../../../libvox/Sources/Audio/AudioMath.h"
#include "../../../libvox/Sources/Core/AppMain.h"
#include "../../../libvox/Sources/Other/CWaveFile.h"
#include "../CDSPResampler.h"

/**
 * @file masstest.cpp
 *
 * @brief Mass randomized/stochastic test of various combinations of sample
 * rate conversions, designed for Dr.Memory debugger.
 *
 * r8brain-free-src Copyright (c) 2013-2023 Aleksey Vaneev
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

void addSine( double* const p, const int l, const double Freq,
	const double SampleRate )
{
	CSineGen sg;
	sg.init( R8B_2PI * Freq / SampleRate, 0.5 );
	int i;

	for( i = 0; i < l; i++ )
	{
		p[ i ] += sg.gen();
	}
}

VOXMAIN
{
	Args.setProgramDescription( "This utility program perform mass "
		"randomized test of sample rate converter." );

	Args.addReqArg( argtFileName, "in-file", "Input WAV filename." );

	VOXCHECK( Args.checkArgs() );

	CWaveFile inf;
	VOXCHECK( inf.loadFile( *Args.getArgValue( "in-file" ).ValueStr ));

	const int InBufSize = (int) min( (int64_t) 50000, inf.SampleCount );

	CStructArray< CFixedBuffer< double > > InBufs( inf.ChannelCount );
	int i;

	for( i = 0; i < inf.ChannelCount; i++ )
	{
		InBufs[ i ].alloc( InBufSize );
	}

	int ReadCount;
	VOXCHECK( inf.readData( InBufs, InBufSize, ReadCount ));

	// Create reference signal which has 93% bandwidth of the input signal.

	CFixedBuffer< double > Ref( InBufSize );

	if( true )
	{
		const int MaxInLen = 521;
		const double tb = 2.0;
		const double bw = 9.3;
		const int Ref0Size = (int) ( InBufSize * bw / 10.0 );
		CFixedBuffer< double > Ref0( Ref0Size );

		CPtrKeeper< r8b :: CDSPResampler > Resamp1;
		Resamp1 = new r8b :: CDSPResampler( 10.0, bw, MaxInLen, tb );

		CPtrKeeper< r8b :: CDSPResampler > Resamp2;
		Resamp2 = new r8b :: CDSPResampler( bw, 10.0, MaxInLen, tb );

		Resamp1 -> oneshot( &InBufs[ 0 ][ 0 ], InBufSize, &Ref0[ 0 ],
			Ref0Size );

		Resamp2 -> oneshot( &Ref0[ 0 ], Ref0Size, &Ref[ 0 ], InBufSize );
	}

	CRnd rnd;
	rnd.init( RndInitValue.get() );

	const int TestCount = 1000;
	CFixedBuffer< double > OutBuf2( InBufSize );
	double peakd = 0.0;
	double maxr = 0.0;
	double avgr = 0.0;
	double avgperf = 0.0;
	double avglatency = 0.0;
	int minld1 = 10000000;
	int maxld1 = -10000000;
	int minld2 = 10000000;
	int maxld2 = -10000000;
	int k;

	for( k = 0; k < TestCount; k++ )
	{
		const double SrcSampleRate = 1.0;
		const double DstSampleRate = 1.0 + 44.0 * rnd.getUniform();
		const double ReqAtten = 180.15;// - 130 * rnd.getUniform();

		const double tb = 0.5 + 4.5 * rnd.getUniform();
		const int MaxInLen = 50 + (int) ( 2000 * rnd.getUniform() );
		const int rlen1 = 1 + (int) ( 1024 * rnd.getUniform() );
		const int rlen2 = 1 + (int) ( 1024 * rnd.getUniform() );

		printf( "Iteration %3i dst=%8.4f tb=%7.4f inlen=%4i ",
			k, DstSampleRate, tb, MaxInLen );

		const int ol = (int) ( InBufSize * DstSampleRate / SrcSampleRate );
		CFixedBuffer< double > OutBuf( ol );

		CPtrKeeper< r8b :: CDSPResampler > Resamp;
		Resamp = new r8b :: CDSPResampler( SrcSampleRate, DstSampleRate,
			MaxInLen, tb, ReqAtten );

		const int inlen1a = Resamp -> getInLenBeforeOutStart( rlen1-1 )+1;
		const int inlen1b = Resamp -> getInputRequiredForOutput( rlen1 );
		avglatency += inlen1a;
		minld1 = min( minld1, inlen1b - inlen1a );
		maxld1 = max( maxld1, inlen1b - inlen1a );

		const TClock t1( CSystem :: getClock() );
		Resamp -> oneshot( &Ref[ 0 ], InBufSize, &OutBuf[ 0 ], ol );
		double perf = 1e-6 * InBufSize /
			CSystem :: getClockDiffSec( t1 );

//		addSine( OutBuf, ol, ( SrcSampleRate + DstSampleRate ) * 0.25,
//			DstSampleRate );

		Resamp = new r8b :: CDSPResampler( DstSampleRate, SrcSampleRate,
			MaxInLen, tb, ReqAtten );

		const int inlen2a = Resamp -> getInLenBeforeOutStart( rlen2-1 )+1;
		const int inlen2b = Resamp -> getInputRequiredForOutput( rlen2 );
		minld2 = min( minld2, inlen2b - inlen2a );
		maxld2 = max( maxld2, inlen2b - inlen2a );

		const TClock t2( CSystem :: getClock() );
		Resamp -> oneshot( &OutBuf[ 0 ], ol, &OutBuf2[ 0 ], InBufSize );
		perf = ( perf + 1e-6 * InBufSize /
			CSystem :: getClockDiffSec( t2 )) * 0.5;

		const double r = calcRMS( &Ref[ 5000 ], &OutBuf2[ 5000 ],
			InBufSize - 10000, peakd );

		if( r > maxr )
		{
			maxr = r;
		}

		avgr += r * r;
		avgperf += perf * perf;

		printf( "z=%7.2f perf=%6.2f %i %i\n", 20.0 * log( r ) / log( 10.0 ),
			perf, inlen1b - inlen1a, inlen2b - inlen2a );
	}

	printf( "avg rms %.2f\n", 10.0 * log( avgr / TestCount ) / log( 10.0 ));
	printf( "max rms %.2f\n", 20.0 * log( maxr ) / log( 10.0 ));
	printf( "peak diff %.2f\n", 20.0 * log( peakd ) / log( 10.0 ));
	printf( "avg perf %.2f Mrops\n", sqrt( avgperf / TestCount ));
	printf( "avg latency %.0f\n", avglatency / TestCount );
	printf( "inlen diffs: %i %i | %i %i\n", minld1, maxld1, minld2, maxld2 );

	VOXRET;
}
