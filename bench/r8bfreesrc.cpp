///$ bin "Win64/r8bfreesrc"

#include "../../../libvox/Sources/Core/AppMain.h"
#include "../../../libvox/Sources/Other/CWaveFile.h"
#include "../CDSPResampler.h"

/**
 * @file r8bfreesrc.cpp
 *
 * @brief Test/benchmarking utility.
 *
 * This is an example file which you won't be able to compile as it includes
 * some undisclosed program code. Please consider this example as a
 * pseudo-code demonstrating the use of the library. Here you can find an
 * example implementation of the simplest sample rate converter utility.
 *
 * r8brain-free-src Copyright (c) 2013-2023 Aleksey Vaneev
 * See the "LICENSE" file for license.
 */

VOXMAIN
{
	Args.setProgramDescription( "This utility program perform sample rate "
		"conversion of the source WAV/Wave64/AIF sound file." );

	Args.addReqArg( argtFileName, "in-file",
		"The input WAV/Wave64/AIF filename." );

	Args.addReqArg( argtFileName, "out-file", "The output WAV filename." );
	Args.addReqArg( argtDouble, "out-sample-rate", "The output sample rate "
		"(e.g. 96000)." );

	Args.addOptArg( argtDouble, "trans-band", "t", "2.0",
		"Transition band in percent." );

	VOXCHECK( Args.checkArgs() );

	CWaveFile inf;
	VOXCHECK( inf.loadFile( *Args.getArgValue( "in-file" ).ValueStr ));

	double OutSampleRate = Args.getArgValue( "out-sample-rate" ).Value.Double;

	if( OutSampleRate <= 0.0 )
	{
		VOXERROR( "@voxstr_r8bfreesrc_InvalidOutSampleRate "
			"The specified output sample rate is invalid." );
	}

	const double InSampleRate = inf.SampleRate;
	const double d = ( OutSampleRate > InSampleRate ?
		OutSampleRate / InSampleRate : InSampleRate / OutSampleRate );

	if( d >= ( 1 << 29 ))
	{
		VOXERROR( "@voxstr_r8bfreesrc_InvalidOutSampleRate2 "
			"The specified output sample rate differs from the source sample "
			"rate too much." );
	}

	CWaveFile outf;
	outf.inheritFormat( inf );
	outf.SampleRate = OutSampleRate;
	outf.inheritCuePoints( inf );
	outf.inheritInfo( inf );

	VOXCHECK( outf.saveFile( *Args.getArgValue( "out-file" ).ValueStr ));

	const int InBufCapacity = 2048;
	CStructArray< CFixedBuffer< double > > InBufs( inf.ChannelCount );
	CInitArray< CPtrKeeper< r8b :: CDSPResampler > > Resamps(
		inf.ChannelCount );

	const double tb = Args.getArgValue( "trans-band" ).Value.Double;

	if( tb < r8b :: CDSPFIRFilter :: getLPMinTransBand() ||
		tb > r8b :: CDSPFIRFilter :: getLPMaxTransBand() )
	{
		VOXERROR( "@voxstr_r8bfreesrc_InvalidTransBand "
			"The specified transition band is invalid." );
	}

	int i;

	for( i = 0; i < inf.ChannelCount; i++ )
	{
		InBufs[ i ].alloc( InBufCapacity );

		Resamps[ i ] = new r8b :: CDSPResampler24( InSampleRate,
			OutSampleRate, InBufCapacity, tb );
	}

	int64_t ol = (int64_t) ( inf.SampleCount * OutSampleRate /
		InSampleRate );

	const int64_t ol0 = inf.SampleCount * inf.ChannelCount;
	int64_t ool = 0;
	double srct = 0.0;
	CArray< double* > opp( inf.ChannelCount );

	while( ol > 0 )
	{
		int ReadCount;
		VOXCHECK( inf.readData( InBufs, InBufCapacity, ReadCount ));

		if( ReadCount == -1 )
		{
			ReadCount = InBufCapacity;

			for( i = 0; i < inf.ChannelCount; i++ )
			{
				memset( &InBufs[ i ][ 0 ], 0, ReadCount * sizeof( double ));
			}
		}

		int WriteCount; // At initial steps this variable can be equal to 0
			// after resampler. Same number for all channels.

		const TClock t1 = CSystem :: getClock();

		for( i = 0; i < inf.ChannelCount; i++ )
		{
			WriteCount = Resamps[ i ] -> process(
				InBufs[ i ], ReadCount, opp[ i ]);
		}

		srct += CSystem :: getClockDiffSec( t1 );
		ool += WriteCount * inf.ChannelCount;

		if( WriteCount > ol )
		{
			WriteCount = (int) ol;
		}

		VOXCHECK( outf.writeData( &opp[ 0 ], WriteCount ));
		ol -= WriteCount;
	}

	VOXCHECK( outf.finalize() );

	printf( "Resampled in %.4f s, %.3f Mrops (excluding IO operations)\n",
		srct, 1e-6 * ol0 / srct );

	VOXRET;
}
