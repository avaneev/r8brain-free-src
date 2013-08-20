/**
 * \file calcErrorTable.cpp
 *
 * \brief Utility program that calculates parameter error table.
 *
 * This utility program calculates tabbed parameter error table for the whole
 * range of transition bands and stop-band attenuations. The error is
 * introduced by approximation functions used in the CDSPFIRFilter class's
 * buildLPFilter() function. If the error is small the "atten_err" should be
 * close to 0 and "gtb" should be close to -3.01. Results will be sent to
 * "stdout".
 *
 * r8brain-free-src Copyright (c) 2013 Aleksey Vaneev
 * See the "License.txt" file for license.
 */

#include <stdio.h>
#include "CDSPFIRFilter.h"
using namespace r8b;

int main()
{
	const double MinAtten = CDSPFIRFilter :: getLPMinAtten();
	const double MaxAtten = CDSPFIRFilter :: getLPMaxAtten();
	const int AttenSteps = 300;
	int i;

	printf( "tb\tatten\tg05\tgtb\tatten_err\n" );

	for( i = 0; i < AttenSteps; i++ )
	{
		const double atten = MinAtten + ( MaxAtten - MinAtten ) * i /
			( AttenSteps - 1 );

		double tb;

		for( tb = CDSPFIRFilter :: getLPMinTransBand();
			tb <= CDSPFIRFilter :: getLPMaxTransBand(); tb *= 1.03 )
		{
			CDSPFIRFilter& f =
				CDSPFIRFilterCache :: getLPFilter( 0.5, tb, atten );

			// Perform inverse FFT to obtain time-domain filter response.

			CDSPRealFFTKeeper ffto( f.getBlockSizeBits() + 1 );
			CFixedBuffer< double > Kernel( 2 << f.getBlockSizeBits() );
			memcpy( &Kernel[ 0 ], f.getKernelBlock(),
				( 2 * sizeof( double )) << f.getBlockSizeBits() );

			ffto -> inverse( Kernel );

			// Calculate spectrum power at points of interest.

			const int l = f.getLatency() * 2 + 1;
			const double g = 10.0 / log( 10.0 );
			double re;
			double im;

			calcFIRFilterResponse( &Kernel[ 0 ], l, M_PI * 0.5, re, im );
			const double g05 = g * log( re * re + im * im );

			calcFIRFilterResponse( &Kernel[ 0 ], l,
				M_PI * ( 1.0 - tb * 0.01 ) * 0.5, re, im );

			// Output the error data.

			printf( "%f\t%f\t%f\t%f\t%f\n", tb, -atten, g05,
				g * log( re * re + im * im ), g05 + atten );

			f.unref();
		}
	}

	return( 0 );
}
