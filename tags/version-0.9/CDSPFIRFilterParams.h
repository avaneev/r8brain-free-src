//$ nocpp

/**
 * @file CDSPFIRFilterParams.h
 *
 * @brief FIR filter parameter storage and calculator class.
 *
 * This file includes low-pass FIR filter parameter storage and calculator
 * class, used by the r8b::CDSPFIRFilter class.
 *
 * r8brain-free-src Copyright (c) 2013 Aleksey Vaneev
 * See the "License.txt" file for license.
 */

#ifndef R8B_CDSPFIRFILTERPARAMS_INCLUDED
#define R8B_CDSPFIRFILTERPARAMS_INCLUDED

#include "CDSPSincFilterGen.h"

namespace r8b {

/**
 * @brief FIR filter parameter storage and calculator class.
 *
 * This class provides storage and calculation functions for a set of
 * "standardized" filter parameters used to build a low-pass FIR filter. These
 * parameters reside in a separate class because they can be calculated in
 * different ways, using different windowing functions, or even via a look-up
 * table (but this is not implemented).
 */

class CDSPFIRFilterParams
{
public:
	CDSPSincFilterGen :: CWindowFunc WinFunc; ///< Windowing function that
		///< should be used to calculate the FIR filter.
		///<
	double hl; ///< Half-band (meaning ReqNormFreq = 0.5) low-pass filter's
		///< kernel length, in samples.
		///<
	double pwr; ///< Windowing function's power factor.
		///<
	double fo; ///< -3.01 dB frequency offset relative to filter's corner
		///< frequency.
		///<

	/**
	 * Function calculates parameters of a low-pass FIR filter windowed by the
	 * Blackman windowing function.
	 *
	 * @param ReqTransBand Required transition band, in percent. See the
	 * CDSPFIRFilterCache::getLPFilter() for details.
	 * @param ReqAtten Required stop-band attenuation, in decibel. See the
	 * CDSPFIRFilterCache::getLPFilter() for details.
	 */

	void calcLPFBlackman( const double ReqTransBand, const double ReqAtten )
	{
		WinFunc = &CDSPSincFilterGen :: calcWindowBlackman;
		const double tb = ReqTransBand * 0.01;
		double atten = -ReqAtten;

		// The following formulas were obtained using the "Eureqa"
		// application. There is no reason to study these formulas as they
		// were auto-generated from the data sets located in the "other"
		// folder.

		if( tb > 0.25 )
		{
			if( atten >= -60.807018 )
			{
				atten -= 3.14263853546384 * cos( 0.116599753023609 * atten ) -
					2.07860492269578;
			}
			else
			if( atten >= -96.842105 )
			{
				atten -= 0.965767779746293 * cos( -0.326099083143211 *
					atten ) - cos( 0.182983338077055 * atten  )* atan2(
					1.1231156701309 + sin( 0.243426743937334 * atten), cosh(
					cos( -0.329140259667326 * atten )) - 0.366224016392942 );
			}
			else
			{
				atten -= 1.46539656220848 * cos( 14.5772431071805 * asinh(
					0.922759254783399 + atten )) - 0.450755513593313 -
					0.00243933927932254 * atten;
			}

			hl = ( 6.58942058096357 + 1.38946776985556 * gauss(
				5.31881128537823 + 0.0531485203972908 * atten ) +
				0.17345968135883 * atten * atan2( 42.1540673302946 + atten,
				63.4823633354858 )) / tb;

			const double hltb = hl * tb;
			pwr = 0.202611118279376 + 0.115970584828984 * hltb +
				0.0229535108403328 * sin( 0.282160865893246 +
				0.130329732713876 * hltb ) - 0.781818138719999 * gauss( gauss(
				sin( 6.1011988083796 + 0.130329732713876 * hltb ) -
				0.136337187201669 * hltb ) - 0.343878017069483 );

			fo = 1.39046581093177 * asinh( 0.63908684575845 + pwr ) / hl +
				cosh( sin( -5400.22053809293 * tb ) / hl ) -
				0.977535105956537 * tb;
		}
		else
		{
			// Apply 4-range correction to the "attenuation" for +/- 1 dB
			// precision, without this correction the error is +/- 3 dB.

			if( atten >= -74.491228 )
			{
				double atten_b = 0.0;

				if( atten >= -50.651267 )
				{
					atten_b = 0.699900982886959 * cos( 0.222669162355806 *
						atten );
				}

				atten -= 0.438075758491934 + 1.72682101344218 * cos(
					0.230227912668262 * atten ) + 0.0102506237901691 * atten *
					sin( 0.334121655800277 * atten );

				atten -= atten_b;
			}
			else
			if( atten >= -100.491228 )
			{
				atten -= 2.20298948096653 * cos( 0.400387943600211 * atten ) +
					0.0162671715330083 * atten * cos( 0.400358395063225 *
					atten ) - 1.95268707933334 - 0.0263263087274513 * atten;
			}
			else
			{
				atten -= 12.4238553852216 * sin( 0.08178837006 * atten ) *
					atan2( 2.477633394, 0.08178837006 * atten ) +
					0.5689731045 * sin( -0.1238783946 * atten ) * atan2(
					0.01319093007 * atten + sin( -0.1238783946 * atten ), cos(
					0.6821178425 + 0.0962288404 * atten )) - 0.1170035732 -
					35.00537255 * sin( 0.08178837006 * atten );
			}

			hl = ( -10.3864830460583 - 0.28561497095004 * atten ) / tb +
				-215.422486627309 / ( tb * atten - 12.2962925122936 * tanh(
				tanh( tb )) - 0.340884649462501 * tb * atten * atan( sin(
				4.22129490874673 - 0.0953516921538414 * atten )));

			const double hltb = hl * tb;
			pwr = 0.201981441582014 + 0.115979123679821 * hltb +
				0.0200435272760573 * sin( 0.234743434639915 +
				0.130876326345944 * hltb ) - 0.780453063766734 * gauss( gauss(
				sin( 6.10822300874378 + 0.130876326345944 * hltb ) -
				0.136709763359556 * hltb ) - 0.339451943355516 );

			fo = atan2( 1.36789104485968 * asinh( 0.681921659358349 + pwr ),
				hl ) + gauss(( 0.681921659358349 * cos( 2.2957524045949 *
				pwr ) + pwr * cos( 2.2957524045949 * pwr )) /
				( 0.590803811989741 * hl + 4.50625485711854 * exp( pwr + sqr(
				pwr )))) - 0.976490268819982 * tb;
		}
	}
};

} // namespace r8b

#endif // R8B_CDSPFIRFILTERPARAMS_INCLUDED
