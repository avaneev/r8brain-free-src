//$ bin "Win64/rmscompare"

#include "/libvox/Sources/Core/AppMain.h"
#include "/libvox/Sources/Other/CWaveFile.h"
#include <stdint.h>

/**
 * @file rmscompare.cpp
 *
 * @brief RMS comparison utility.
 *
 * r8brain-free-src Copyright (c) 2018 Aleksey Vaneev
 * See the "License.txt" file for license.
 */

VOXMAIN
{
	Args.setProgramDescription( "This utility program perform RMS comparison "
		"of two WAV files for the purposes of sample rate conversion "
		"implementation quality assurance." );

	Args.addReqArg( argtFileName, "in-file1", "Input WAV filename 1." );
	Args.addReqArg( argtFileName, "in-file2", "Input WAV filename 2." );

	VOXCHECK( Args.checkArgs() );

	CWaveFile inf1;
	VOXCHECK( inf1.loadFile( *Args.getArgValue( "in-file1" ).ValueStr ));

	CWaveFile inf2;
	VOXCHECK( inf2.loadFile( *Args.getArgValue( "in-file2" ).ValueStr ));

	if( inf1.ChannelCount != inf2.ChannelCount )
	{
		VOXERROR( "@voxstr_rmscompare_ChannelCountMismatch "
			"Input files' channel count mismatch." );
	}

	if( inf1.SampleRate != inf2.SampleRate )
	{
		VOXERROR( "@voxstr_rmscompare_SampleRateMismatch "
			"Input files' sample rate mismatch." );
	}

	const int InBufCapacity = 4096;
	CInitArray< CFixedBuffer< double > > InBufs1( inf1.ChannelCount );
	CInitArray< CFixedBuffer< double > > InBufs2( inf1.ChannelCount );
	CArray< double > rms( inf1.ChannelCount );
	int i;

	for( i = 0; i < inf1.ChannelCount; i++ )
	{
		InBufs1[ i ].alloc( InBufCapacity );
		InBufs2[ i ].alloc( InBufCapacity );
		rms[ i ] = 0.0;
	}

	int64_t ll;

	if( inf1.SampleCount != inf2.SampleCount )
	{
		printf( "NOTE: length of files mismatch, using shortest file\n" );

		ll = ( inf1.SampleCount < inf2.SampleCount ?
			inf1.SampleCount : inf2.SampleCount );
	}
	else
	{
		ll = inf1.SampleCount;
	}

	const int64_t ll0 = ll;
	printf( "Total samples: %lli (%.4f s)\n", ll,
		(double) ll / inf1.SampleRate );

	while( ll > 0 )
	{
		int ReadCount1;
		VOXCHECK( inf1.readData( InBufs1, InBufCapacity, ReadCount1 ));

		int ReadCount2;
		VOXCHECK( inf2.readData( InBufs2, InBufCapacity, ReadCount2 ));

		if( ReadCount1 == -1 || ReadCount2 == -1 )
		{
			VOXERROR( "@voxstr_rmscompare_PrematureEOF "
				"Premature EOF." );
		}

		if( ReadCount1 != ReadCount2 )
		{
			VOXERROR( "@voxstr_rmscompare_ReadCountMismatch "
				"Read count mismatch." );
		}

		for( i = 0; i < inf1.ChannelCount; i++ )
		{
			const double* ip1 = InBufs1[ i ];
			const double* ip2 = InBufs2[ i ];
			int k;

			for( k = 0; k < ReadCount1; k++ )
			{
				const double d = ip1[ k ] - ip2[ k ];
				rms[ i ] += d * d;
			}
		}

		ll -= ReadCount1;
	}

	for( i = 0; i < inf1.ChannelCount; i++ )
	{
		printf( "Channel %i RMS of difference: %.2f dB\n", i + 1,
			10.0 * log( rms[ i ] / ll0 ) / log( 10.0 ));
	}

	VOXRET;
}
