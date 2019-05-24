// Window function parameters optimizer - finds parameters for the specified
// filter length for Vaneev and Kaiser window functions using the BiteOptDeep
// derivative-free optimization method.

#include <stdio.h>
#include <math.h>
#include "/projects/biteopt/biteopt.h"
#include "../CDSPSincFilterGen.h"

const bool IsKaiser = true; // Is Kaiser window function?
const double CornerFreq = 1.0; // Corner frequency.
const double LinFraction = 0.5; // Linear part of the spectrum.
const double StopFraction = 1.5; // Stop-band part of the spectrum.
const double StopFractionEnd = 4.0; // Stop-band part of the spectrum end.
int FilterLen = 24; // Filter length.
const int Oversample = 20; // Spectrum oversampling ratio.
const int OptDepth = 9; // BiteOptDeep optimizer depth.

class CWinOpt : public CBiteOptDeep
{
public:
	CWinOpt()
	{
		updateDims(( IsKaiser ? 2 : 5 ), OptDepth );
	}

	virtual void getMinValues( double* const p ) const
	{
		if( IsKaiser )
		{
			p[ 0 ] = 2.0;
			p[ 1 ] = 0.5;
		}
		else
		{
			p[ 0 ] = 0.0;
			p[ 1 ] = 0.0;
			p[ 2 ] = 0.0;
			p[ 3 ] = 0.0;
			p[ 4 ] = 0.0;
		}
	}

	virtual void getMaxValues( double* const p ) const
	{
		if( IsKaiser )
		{
			p[ 0 ] = 30.0;
			p[ 1 ] = 3.0;
		}
		else
		{
			p[ 0 ] = 4.0;
			p[ 1 ] = 4.0;
			p[ 2 ] = 4.0;
			p[ 3 ] = 4.0;
			p[ 4 ] = 4.0;
		}
	}

	virtual double optcost( const double* const p )
	{
		r8b :: CDSPSincFilterGen gen;
		gen.Freq1 = 0.0;
		gen.Freq2 = M_PI * CornerFreq / Oversample;
		gen.Len2 = FilterLen * 0.5 * Oversample;
		gen.initBand(( IsKaiser ?
			r8b :: CDSPSincFilterGen :: wftKaiser :
			r8b :: CDSPSincFilterGen :: wftVaneev ), p, IsKaiser );

		const int KernelLen = gen.KernelLen;

		double KernelBlock[ KernelLen ];
		gen.generateBand( &KernelBlock[ 0 ], ( IsKaiser ?
			&r8b :: CDSPSincFilterGen :: calcWindowKaiser :
			&r8b :: CDSPSincFilterGen :: calcWindowVaneev ));

		r8b :: normalizeFIRFilter( &KernelBlock[ 0 ], KernelLen, 1.0 );

		const double _10ln10 = 10.0 / log( 10.0 );
		const int Count1 = 500;
		double cost1 = 0.0;
		double re;
		double im;
		int i;

		for( i = 0; i <= Count1; i++ )
		{
			const double th = M_PI * LinFraction / Oversample * i / Count1;
			r8b :: calcFIRFilterResponse( KernelBlock, KernelLen, th, re, im );
			double p = fabs( _10ln10 * log( re * re + im * im ));

			cost1 = max( p, cost1 );
		}

		const int Count2 = 2000;
		double cost2 = -1000.0;
		const double th1 = M_PI * StopFraction / Oversample;
		const double th2 = M_PI * StopFractionEnd / Oversample;

		for( i = 0; i <= Count2; i++ )
		{
			const double th = th1 + ( th2 - th1 ) * i / Count2;
			r8b :: calcFIRFilterResponse( KernelBlock, KernelLen, th, re, im );
			double p = _10ln10 * log( re * re + im * im );

			cost2 = max( p, cost2 );
		}

		return( cost1 * 36.0 + cost2 );
	}
};

int main()
{
	CBiteRnd rnd;
	rnd.init( 1 );

	CWinOpt opt;

	for( FilterLen = 6; FilterLen <= 30; FilterLen += 2 )
	{
		opt.init( rnd );

		const int sct = 15 * opt.getInitEvals() / OptDepth; // Plateau thresh.
		int i;

		for( i = 0; i < 1000000; i++ )
		{
			if( opt.optimize( rnd ) > sct )
			{
				break;
			}
		}

		if( IsKaiser )
		{
			printf( "\t\t\t{ %.16f, %.16f }, // %i @ %.2f\n",
				opt.getBestParams()[ 0 ], opt.getBestParams()[ 1 ], FilterLen,
				opt.getBestCost() );
		}
		else
		{
			printf( "\t\t\t{ %.8f, %.8f, %.8f, %.8f, %.8f }, // %i @ %.2f\n",
				opt.getBestParams()[ 0 ], opt.getBestParams()[ 1 ],
				opt.getBestParams()[ 2 ], opt.getBestParams()[ 3 ],
				opt.getBestParams()[ 4 ], FilterLen, opt.getBestCost() );
		}
	}
}
