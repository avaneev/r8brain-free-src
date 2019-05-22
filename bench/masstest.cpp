#include "/libvox/Sources/Audio/AudioMath.h"
#include "/libvox/Sources/Core/AppMain.h"
#include "/libvox/Sources/Other/CWaveFile.h"
#include "../CDSPResampler.h"

/**
 * @file masstest.cpp
 *
 * @brief Mass randomized test of various combinations of sample rate
 * conversions, designed for Dr.Memory debugger.
 *
 * r8brain-free-src Copyright (c) 2019 Aleksey Vaneev
 * See the "License.txt" file for license.
 */

VOXMAIN
{
	Args.setProgramDescription( "This utility program perform mass "
		"randomized test of sample rate converter." );

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

	CRnd rnd;
	rnd.init( getRndInitValue() );

	const int TestCount = 1000;
	int k;

	for( k = 0; k < TestCount; k++ )
	{
		const double SrcSampleRate = (int) ( 10 + 200 * rnd.getUniform() );
		const double DstSampleRate = (int) ( 10 + 200 * rnd.getUniform() );
		const double tb = 1.0 + 9.0 * rnd.getUniform();
		const int MaxInLen = 50 + (int) ( 2000 * rnd.getUniform() );

		printf( "Iteration %i src=%f dst=%f tb=%f MaxInLen=%i\n",
			k + 1, SrcSampleRate, DstSampleRate, tb, MaxInLen );

		const int ol = (int) ( InBufSize * SrcSampleRate / DstSampleRate );
		CFixedBuffer< double > OutBuf( ol );
		CPtrKeeper< r8b :: CDSPResampler24* > Resamp;
		Resamp = new r8b :: CDSPResampler24( SrcSampleRate,
			DstSampleRate, MaxInLen, tb );

		Resamp -> oneshot( MaxInLen, &InBufs[ 0 ][ 0 ], InBufSize,
			&OutBuf[ 0 ], ol );
	}

	VOXRET;
}
