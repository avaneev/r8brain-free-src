///$ bin "Win64/r8bfreesrc"

#include "../../../libvox/Sources/Core/AppMain.h"
#include "../../../libvox/Sources/Other/CWaveFile.h"
#include "../../../libvox/Sources/Audio/AudioMath.h"
#include "../../../libvox/Sources/Audio/DCDAW/DCDAWBiquadRBJ.h"
#include "../CDSPResampler.h"

/**
 * @file sacd.cpp
 *
 * @brief 
 *
 * 
 *
 * r8brain-free-src Copyright (c) 2013-2019 Aleksey Vaneev
 * See the "License.txt" file for license.
 */

class CDither1Bit
{
public:
	CDither1Bit( const double SampleRate )
	{
		const double Freq = 48000.0;
		const double Gain = 6.02 * 14;
		const double BW = 1.7;
		CDCDAWBiquadRBJ flt;
		flt.setParams( SampleRate, Freq, pow( 10.0, Gain * 0.05 ), BW,
			false, true );

		flt.calcHSH();

		b1a = flt.b1 / flt.b0;
		b2a = flt.b2 / flt.b0;
		a1a = flt.a1 / flt.a0;
		a2a = flt.a2 / flt.a0;

		y1a = y2a = 0.0;

		flt.setParams( SampleRate, 200.0, pow( 10.0, -78.0 * 0.05 ), 7.0,
			false, true );

		flt.calcLSH();

		b1b = flt.b1 / flt.b0;
		b2b = flt.b2 / flt.b0;
		a1b = flt.a1 / flt.a0;
		a2b = flt.a2 / flt.a0;

		y1b = y2b = 0.0;

		BitMult = 1.0;
		BitMultI = 1.0 / BitMult;

		rnd.init( RndInitValue.get() );

		prev = 0.0;
		prev2 = 0.0;
	}

	void process( double* const p, const int l )
	{
		int i;

		for( i = 0; i < l; i++ )
		{
			double v = ( p[ i ] - prev ) * 64;
			prev = p[ i ];
			const double fltout = y1a + y1b;
			double res = v + fltout + rnd.get1BitCentered();

			if( res >= 0.0 )
			{
				res = 1.0;
			}
			else
			{
				res = -1.0;
			}

			const double y0b = res - v;
			const double y0a = y0b - y1b;
			const double x0a = y0a - y1a;

			y1a = b1a * x0a - a1a * y0a + y2a;
			y2a = b2a * x0a - a2a * y0a;
			y1b = b1b * y0a - a1b * y0b + y2b;
			y2b = b2b * y0a - a2b * y0b;

			prev2 += res/64;
			p[ i ] = prev2/2;
		}
	}

protected:
	double BitMult; ///< Bit depth multiplier.
	double BitMultI; ///< Bit depth divisor.
	CRnd rnd; ///< Random number generator.
	double b1a, b2a, a1a, a2a;
	double y1a, y2a;
	double b1b, b2b, a1b, a2b;
	double y1b, y2b;
	double prev;
	double prev2;
};

VOXMAIN
{
	Args.setProgramDescription( "This utility program perform sample rate "
		"conversion of the source WAV/Wave64/AIF sound file." );

	Args.addReqArg( argtFileName, "in-file",
		"The input WAV/Wave64/AIF filename." );

	Args.addReqArg( argtFileName, "out-file", "The output WAV filename." );
	Args.addOptArg( argtDouble, "trans-band", "t", "2.0",
		"Transition band in percent." );

	VOXCHECK( Args.checkArgs() );

	CWaveFile inf;
	VOXCHECK( inf.loadFile( *Args.getArgValue( "in-file" ).ValueStr ));

	const double OutSampleRate = 2822400.0;
	const double InSampleRate = inf.SampleRate;
	const double OutSampleRate2 = 88200.0*2;

	CWaveFile outf;
	outf.inheritFormat( inf );
	outf.setSampleDepth( 32, true );
	outf.SampleRate = InSampleRate;//OutSampleRate;
	outf.inheritCuePoints( inf );
	outf.inheritInfo( inf );

	VOXCHECK( outf.saveFile( *Args.getArgValue( "out-file" ).ValueStr ));

	const int InBufCapacity = 1024;
	CInitArray< CFixedBuffer< double > > InBufs( inf.ChannelCount );
	CInitArray< CPtrKeeper< r8b :: CDSPResampler24* > > Resamps(
		inf.ChannelCount );

	CInitArray< CPtrKeeper< CDither1Bit* > > Dithers( inf.ChannelCount );

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

		Dithers[ i ] = new CDither1Bit( OutSampleRate );
	}

	int64_t ol = (int64_t) ( inf.SampleCount * OutSampleRate /
		InSampleRate );

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

		for( i = 0; i < inf.ChannelCount; i++ )
		{
			WriteCount = Resamps[ i ] -> process(
				InBufs[ i ], ReadCount, opp[ i ]);

			Dithers[ i ] -> process( opp[ i ], WriteCount );
		}

		if( WriteCount > ol )
		{
			WriteCount = (int) ol;
		}

		VOXCHECK( outf.writeData( &opp[ 0 ], WriteCount ));
		ol -= WriteCount;
	}

	VOXCHECK( inf.closeInFile() );
	VOXCHECK( outf.finalize() );

	VOXCHECK( inf.loadFile( *Args.getArgValue( "out-file" ).ValueStr ));

	outf.inheritFormat( inf );
	outf.SampleRate = OutSampleRate2;
	outf.inheritCuePoints( inf );
	outf.inheritInfo( inf );

	VOXCHECK( outf.saveFile( "DrumSrc2Res.wav" ));

	for( i = 0; i < inf.ChannelCount; i++ )
	{
		Resamps[ i ] = new r8b :: CDSPResampler24( OutSampleRate,
			OutSampleRate2, InBufCapacity, tb );
	}

	ol = (int64_t) ( inf.SampleCount * OutSampleRate2 / OutSampleRate );

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

		for( i = 0; i < inf.ChannelCount; i++ )
		{
			int z;

			for( z = 0; z < ReadCount; z++ )
			{
				InBufs[ i ][ z ] *= 2.0;
			}

			WriteCount = Resamps[ i ] -> process(
				InBufs[ i ], ReadCount, opp[ i ]);
		}

		if( WriteCount > ol )
		{
			WriteCount = (int) ol;
		}

		VOXCHECK( outf.writeData( &opp[ 0 ], WriteCount ));
		ol -= WriteCount;
	}

	VOXCHECK( inf.closeInFile() );
	VOXCHECK( outf.finalize() );

	VOXRET;
}
