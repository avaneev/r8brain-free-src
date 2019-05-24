// Sinewave resampling test - checks for sample rate conversion precision by
// upsampling a reference sinewave signal. Requires R8B_FLTTEST. Used to find
// optimal FilterFracs values for different filter lengths. When the
// calculated "Avg diff" plateaus it means the optimal FilterFracs value was
// reached.

#include <stdio.h>
#include "../CDSPResampler.h"

typedef r8b :: CDSPResampler24 CResamp;

const double InSampleRate = 44100.0;
const double SineFreq = 10000.0;
const int InBufSize = 70000;
r8b :: CFixedBuffer< double > InBuf( InBufSize );

double calcSineRMS( const double* const p1, const int l, const int o,
	const double SrcSampleRate, const double DstSampleRate )
{
	double s = 0.0;
	int i;

	for( i = 0; i < l; i++ )
	{
		const double v = sin( i * 2.0 * M_PI * SineFreq * SrcSampleRate /
			( InSampleRate * DstSampleRate ));

		if( i >= o )
		{
			const double d = p1[ i ] - v;
			s += d * d;
		}
	}

	return( sqrt( s / ( l - o )));
}

int main()
{
	// Create reference sinewave signal.

	if( true )
	{
		int i;

		for( i = 0; i < InBufSize; i++ )
		{
			InBuf[ i ] = sin( i * 2.0 * M_PI * SineFreq / InSampleRate );
		}
	}

	for( r8b :: InterpFilterFracs = 21; r8b :: InterpFilterFracs<= 800;
		r8b :: InterpFilterFracs += 9 )
	{
		const double SrcSampleRate = 20.0;
		double avgd = 0.0;
		double maxd = 0.0;
		int avgc = 0;
		int k;

		for( k = 21; k < 600; k += 7 )
		{
			const int MaxInLen = 521;
			const double tb = 5.0;

			const double DstSampleRate = k;
			const int ol1 = (int) ( InBufSize * DstSampleRate /
				SrcSampleRate );

			r8b :: CFixedBuffer< double > OutBuf1( ol1 );

			r8b :: CPtrKeeper< CResamp* > Resamp1;
			Resamp1 = new CResamp( SrcSampleRate, DstSampleRate, MaxInLen,
				tb );

			Resamp1 -> oneshot( MaxInLen, &InBuf[ 0 ], InBufSize,
				&OutBuf1[ 0 ], ol1 );

			const int o = (int) ( ol1 * 0.2 );
			const double r = calcSineRMS( OutBuf1, ol1 - o, o,
				SrcSampleRate, DstSampleRate );

			avgd += r * r;
			maxd = max( maxd, r );
			avgc++;
		}

		printf( "InterpFilterLen %i ", r8b :: InterpFilterFracs );
		printf( "Avg diff %.2f ", 10.0 * log( avgd / avgc ) / log( 10.0 ));
		printf( "Max diff %.2f\n", 20.0 * log( maxd ) / log( 10.0 ));
	}
}
