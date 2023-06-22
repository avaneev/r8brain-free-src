#include "../../../libvox/Sources/Audio/AudioMath.h"
#include "../../../libvox/Sources/Core/AppMain.h"
#include "../CDSPResampler.h"

/**
 * @file snrtest.cpp
 *
 * @brief Mass SNR precision test - checks for sample rate conversion
 * precision at the ReqAtten levels. At some ReqAtten settings resampler
 * delivers a lower SNR, thus setting ReqAtten+9 is advisable.
 *
 * r8brain-free-src Copyright (c) 2013-2023 Aleksey Vaneev
 * See the "LICENSE" file for license.
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

const int InBufSize = 100000;
const int BufSkip = 20000;
CFixedBuffer< double > InBuf( InBufSize );

VOXMAIN
{
	CRnd rnd;
	rnd.init( 1000000 );

	int i;

	for( i = 0; i < InBufSize; i++ )
	{
		InBuf[ i ] = ( rnd.getUniform() - 0.5 ) * 2.0;
	}

	CFixedBuffer< double > Ref( InBufSize );

	if( true )
	{
		const int MaxInLen = 521;
		const double tb = 2.0;
		const int Ref0Size = (int) ( InBufSize * 9.0 / 10.0 );
		CFixedBuffer< double > Ref0( Ref0Size );

		CPtrKeeper< r8b :: CDSPResampler24 > Resamp1;
		Resamp1 = new r8b :: CDSPResampler24( 10.0, 9.0, MaxInLen, tb );

		CPtrKeeper< r8b :: CDSPResampler24 > Resamp2;
		Resamp2 = new r8b :: CDSPResampler24( 9.0, 10.0, MaxInLen, tb );

		Resamp1 -> oneshot( &InBuf[ 0 ], InBufSize, &Ref0[ 0 ], Ref0Size );
		Resamp2 -> oneshot( &Ref0[ 0 ], Ref0Size, &Ref[ 0 ], InBufSize );
	}

	const double SrcSampleRate = 20.0;
	double ReqAtten;

	for( ReqAtten = 49; ReqAtten <= 218; ReqAtten += 6 )
	{
		double avgd = 0.0;
		double maxd = 0.0;
		int avgc = 0;
		int k;

		for( k = 21; k < 600; k += 7 )
		{
			const int MaxInLen = (int) ( 50 + rnd.getUniform() * 1500 );
			const double tb = 0.5 + rnd.getUniform() * 4.0;

			const double DstSampleRate = k;
			const int ol1 = (int) ( InBufSize * DstSampleRate / SrcSampleRate );
			CFixedBuffer< double > OutBuf1( ol1 );
			CFixedBuffer< double > OutBuf2( InBufSize );

			CPtrKeeper< r8b :: CDSPResampler > Resamp1;
			Resamp1 = new r8b :: CDSPResampler( SrcSampleRate, DstSampleRate,
				MaxInLen, tb, ReqAtten );

			CPtrKeeper< r8b :: CDSPResampler > Resamp2;
			Resamp2 = new r8b :: CDSPResampler( DstSampleRate, SrcSampleRate,
				MaxInLen, tb, ReqAtten );

			Resamp1 -> oneshot( &Ref[ 0 ], InBufSize, &OutBuf1[ 0 ], ol1 );
			Resamp2 -> oneshot( &OutBuf1[ 0 ], ol1, &OutBuf2[ 0 ],
				InBufSize );

			const double r = calcRMS( &Ref[ BufSkip ], &OutBuf2[ BufSkip ],
				InBufSize - BufSkip * 2 );

			avgd += r * r;
			maxd = max( maxd, r );
			avgc++;
		}

		printf( "ReqAtten=%.2f avg %.2f max %.2f\n", ReqAtten,
			10.0 * log( avgd / avgc ) / log( 10.0 ),
			20.0 * log( maxd ) / log( 10.0 ));
	}

	VOXRET;
}
