/**
 * \file example.cpp
 *
 * \brief An example C++ file that demonstrates resampler's usage.
 *
 * This is an example file which you won't be able to compile as it includes
 * some undisclosed program code. Please consider this example as a
 * pseudo-code demonstrating the use of the library.
 *
 * r8brain-free-src Copyright (c) 2013 Aleksey Vaneev
 * See the "License.txt" file for license.
 */

int main()
{
	const double OutSampleRate = 96000.0;

	CWaveFile inf;
	inf.loadFile( "Audio.wav" );

	CWaveFile outf;
	outf.inheritFormat( inf );
	outf.SampleRate = OutSampleRate;
	outf.saveFile( "Audio2.wav" );

	const int InBufCapacity = 2048;
	CFixedBuffer< double > InBufs[ inf.ChannelCount ];
	CPtrKeeper< CDSPResampler* > Resamps[ inf.ChannelCount ];

	int i;

	for( i = 0; i < inf.ChannelCount; i++ )
	{
		InBufs[ i ].alloc( InBufCapacity );

		Resamps[ i ] = new CDSPResampler( inf.SampleRate, OutSampleRate,
			2, 96, InBufCapacity );
	}

	int64_t ol = inf.SampleCount * OutSampleRate / inf.SampleRate;

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

		const int b = (int) ( ol < ReadCount ? ol : ReadCount );
		double* opp[ inf.ChannelCount ];
		int wc;

		for( i = 0; i < inf.ChannelCount; i++ )
		{
			wc = Resamps[ i ] -> process( InBufs[ i ], b, opp[ i ]);
		}

		outf.writeData( opp, wc );
		ol -= wc;
	}

	outf.finalize();

	return( 0 );
}
