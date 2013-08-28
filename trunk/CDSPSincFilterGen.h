//$ nocpp

/**
 * \file CDSPSincFilterGen.h
 *
 * \brief Sinc function-based FIR filter generator class.
 *
 * This file includes the CDSPSincFilterGen class implementation that
 * generates FIR filters.
 *
 * r8brain-free-src Copyright (c) 2013 Aleksey Vaneev
 * See the "License.txt" file for license.
 */

#ifndef R8B_CDSPSINCFILTERGEN_INCLUDED
#define R8B_CDSPSINCFILTERGEN_INCLUDED

#include "r8bbase.h"

namespace r8b {

/**
 * \brief Sinc function-based FIR filter generator class.
 *
 * Structure that holds state used to perform generation of sinc functions of
 * various types, windowed by the Blackman window by default (but the
 * windowing function can be changed if necessary).
 */

class CDSPSincFilterGen
{
public:
	double Len2; ///< Required half filter kernel's length in samples (can be
		///< a fractional value). Final physical kernel length will be
		///< provided in the KernelLen variable. Len2 should be >= 2.
		///<
	double Freq1; ///< Required corner circular frequency 1 [0; pi]. Used only
		///< in the generateBand() function.
		///<
	double Freq2; ///< Required corner circular frequency 2 [0; pi]. Used only
		///< in the generateBand() function. The range [Freq1; Freq2] defines
		///< a pass band for the generateBand() function.
		///<
	double FracDelay; ///< Fractional delay in the range [0; 1], used only
		///< by the generateFrac() function. Note that the FracDelay parameter
		///< is actually inversed. At 0.0 value it produces 1 sample delay
		///< (with latency equal to fl2), at 1.0 value it produces 0 sample
		///< delay (with latency equal to fl2 - 1).
		///<
	int KernelLen; ///< Resulting length of the filter kernel, this variable
		///< is set after the call to one of the "init" functions.
		///<
	int fl2; ///< Internal "half kernel length" value. This value can be used
		///< as filter's latency in samples (taps), this variable is set after
		///< the call to one of the "init" functions.
		///<

	typedef double( CDSPSincFilterGen :: *CWindowFunc )(); ///< Windowing
		///< function pointer type.
		///<

	/**
	 * Function initializes *this structure for generation of band-limited
	 * sinc filter kernel. The generateBand() or generateBandPow() functions
	 * should be used to calculate the filter.
	 */

	void initBand()
	{
		R8BASSERT( Len2 >= 2.0 );

		fl2 = (int) floor( Len2 );
		KernelLen = fl2 + fl2 + 1;
		f1.init( Freq1, 2.0, -Freq1 * fl2 );
		f2.init( Freq2, 2.0, -Freq2 * fl2 );

		const double step1 = M_PI / Len2;
		w1.init( step1, 2.0, M_PI * 0.5 - step1 * fl2 );

		const double step2 = M_2PI / Len2;
		w2.init( step2, 2.0, M_PI * 0.5 - step2 * fl2 );

		const double step3 = M_3PI / Len2;
		w3.init( step3, 2.0, M_PI * 0.5 - step3 * fl2 );
	}

	/**
	 * Function initializes *this structure for Hilbert transformation filter
	 * calculation. Freq1 and Freq2 variables are not used.
	 * The generateHilbert() function should be used to calculate the filter.
	 */

	void initHilbert()
	{
		R8BASSERT( Len2 >= 2.0 );

		fl2 = (int) floor( Len2 );
		KernelLen = fl2 + fl2 + 1;

		const double step1 = M_PI / Len2;
		w1.init( step1, 2.0, M_PI * 0.5 - step1 * fl2 );

		const double step2 = M_2PI / Len2;
		w2.init( step2, 2.0, M_PI * 0.5 - step2 * fl2 );

		const double step3 = M_3PI / Len2;
		w3.init( step3, 2.0, M_PI * 0.5 - step3 * fl2 );
	}

	/**
	 * Function initializes *this structure for generation of full-bandwidth
	 * fractional delay sinc filter kernel. Freq1 and Freq2 variables are not
	 * used. The generateFrac() function should be used to calculate the
	 * filter.
	 */

	void initFrac()
	{
		R8BASSERT( Len2 >= 2.0 );

		fl2 = (int) ceil( Len2 );
		KernelLen = fl2 + fl2;

		const double step1 = M_PI / Len2;
		w1.init( step1, 2.0, M_PI * 0.5 - step1 * fl2 + step1 * FracDelay );

		const double step2 = M_2PI / Len2;
		w2.init( step2, 2.0, M_PI * 0.5 - step2 * fl2 + step2 * FracDelay );

		const double step3 = M_3PI / Len2;
		w3.init( step3, 2.0, M_PI * 0.5 - step3 * fl2 + step3 * FracDelay );
	}

	/**
	 * This function is similar to the initFrac() function, but initializes
	 * cosine wave generators for use with the calcWindowVaneev() windowing
	 * function. The generateFrac() function should be used to calculate the
	 * filter, with the calcWindowVaneev() function as parameter.
	 *
	 * @param FilterLen Filter length in use, must be an even value, in the
	 * range 6 to 38.
	 * @param Params Externally-provides set of parameters. NULL - use table
	 * values. Non-NULL value is usually provided by optimization algorithm.
	 */

	void initFracVaneev( const int FilterLen, const double* Params = NULL )
	{
		R8BASSERT( Len2 >= 2.0 );

		fl2 = (int) ceil( Len2 );
		KernelLen = fl2 + fl2;

		// This set of parameters was obtained via probabilistic optimization.
		// The number after @ shows the noise floor level for the given filter
		// length. The signal-to-noise (SNR) ratio approximately equals noise
		// floor plus 85. SNR can be also decreased by using a filter bank
		// with suboptimal number of interpolated fractional delay filters:
		// thus select FilterFracs with care.

		static const double Coeffs[][ 4 ] = {
			{ 0.78181357, 0.25515445, 0.14160659, 0.00377970 }, // 6 @ -138.4
			{ 0.82842159, 0.46796048, 0.03723786, 0.01647440 }, // 8 @ -154.1
			{ 0.87106048, 0.55487482, 0.11407443, 0.06108434 }, // 10 @ -167.6
			{ 0.90233911, 0.62228911, 0.17793224, 0.09128849 }, // 12 @ -181.0
			{ 0.92888330, 0.67854727, 0.18462272, 0.15750025 }, // 14 @ -194.5
			{ 0.94714902, 0.71154551, 0.36654152, 0.01684227 }, // 16 @ -208.1
			{ 0.95678003, 0.75499443, 0.40432469, 0.09098344 }, // 18 @ -223.1
			{ 0.94601085, 0.78526426, 0.47076653, 0.10718288 }, // 20 @ -234.1
			{ 0.95382623, 0.80840395, 0.50560712, 0.14820130 }, // 22 @ -248.4
			{ 0.79587584, 0.92850223, 0.54364506, 0.17823490 }, // 24 @ -261.5
			{ 0.81692363, 0.93690700, 0.56559066, 0.21066054 }, // 26 @ -273.0
			{ 0.97006382, 0.85948381, 0.58057290, 0.24934339 }, // 28 @ -285.3
			{ 0.97022013, 0.86906266, 0.59856588, 0.28601037 }, // 30 @ -297.6
			{ 0.97374957, 0.88090262, 0.62906580, 0.29736763 }, // 32 @ -309.9
			{ 0.81824544, 0.95632121, 0.67943518, 0.35354651 }, // 34 @ -325.4
			{ 0.83224964, 0.96134546, 0.68783822, 0.37241201 }, // 36 @ -335.0
			{ 0.60217397, 0.96734066, 0.76693619, 0.40644160 }, // 38 @ -347.0
		};

		double p[ 4 ];

		if( Params == NULL )
		{
			Params = Coeffs[ FilterLen / 2 - 3 ];
		}
		else
		{
			p[ 0 ] = clampr( Params[ 0 ], -4.0, 4.0 );
			p[ 1 ] = clampr( Params[ 1 ], -4.0, 4.0 );
			p[ 2 ] = clampr( Params[ 2 ], -4.0, 4.0 );
			p[ 3 ] = clampr( Params[ 3 ], -4.0, 4.0 );
			Params = p;
		}

		const double step1 = Params[ 0 ] * M_PI / Len2;
		w1.init( step1, 2.0, M_PI * 0.5 - step1 * fl2 + step1 * FracDelay );

		const double step2 = Params[ 1 ] * M_PI / Len2;
		w2.init( step2, 2.0, M_PI * 0.5 - step2 * fl2 + step2 * FracDelay );

		const double step3 = Params[ 2 ] * M_PI / Len2;
		w3.init( step3, 2.0, M_PI * 0.5 - step3 * fl2 + step3 * FracDelay );

		const double step4 = Params[ 3 ] * M_PI / Len2;
		w4.init( step4, 2.0, M_PI * 0.5 - step4 * fl2 + step4 * FracDelay );
	}

	/**
	 * @return The next "Hann" windowing function coefficient.
	 */

	double calcWindowHann()
	{
		return( 0.5 + 0.5 * w1.gen() );
	}

	/**
	 * @return The next "Hann raised to power 6" windowing function
	 * coefficient. This windowing function is suitable for fractional delay
	 * filters.
	 */

	double calcWindowHann6()
	{
		return( sqr( sqr( sqr( calcWindowHann() ))));
	}

	/**
	 * @return The next "Vaneev" windowing function coefficient, for use with
	 * the fractional delay filters. This windowing function, considering the
	 * optimized parameters, is suitable for fractional delay filters even
	 * better than the "Hann^6" function.
	 */

	double calcWindowVaneev()
	{
		const double v0 = 0.5 + 0.5 * w1.gen();
		const double v1 = 0.5 + 0.5 * w2.gen();
		const double v2 = 0.5 + 0.5 * w3.gen();
		const double v3 = 0.5 + 0.5 * w4.gen();

		return( v0 * sqr( v1 ) * sqr( sqr( v2 )) * sqr( sqr( sqr( v3 ))));
	}

	/**
	 * @return The next "Hamming" windowing function coefficient.
	 */

	double calcWindowHamming()
	{
		return( 0.54 + 0.46 * w1.gen() );
	}

	/**
	 * @return The next "Blackman" windowing function coefficient.
	 */

	double calcWindowBlackman()
	{
		return( 0.42 + 0.5 * w1.gen() + 0.08 * w2.gen() );
	}

	/**
	 * @return The next "Nuttall" windowing function coefficient.
	 */

	double calcWindowNuttall()
	{
		return( 0.355768 + 0.487396 * w1.gen() + 0.144232 * w2.gen() +
			0.012604 * w3.gen() );
	}

	/**
	 * @return The next "Nuttall squared" windowing function coefficient. This
	 * window is good for high dynamic range spectral analysis.
	 */

	double calcWindowNuttall2()
	{
		return( sqr( calcWindowNuttall() ));
	}

	/**
	 * @return The next "Blackman-Nuttall" windowing function coefficient.
	 * This (together with the "Nuttall") windowing function provides quite
	 * short filters for the given transition band and stop-band attenuation,
	 * 10-20% shorter than the Blackman windowing function. However, the
	 * aliased components if the filter is used for resampling are stronger
	 * than what "Blackman" function offers (-55 dB in "Blackman" vs -30 dB
	 * in "Nuttall" vs -5 dB in "Blackman-Nuttall").
	 */

	double calcWindowBlackmanNuttall()
	{
		return( 0.3635819 + 0.4891775 * w1.gen() + 0.1365995 * w2.gen() +
			0.0106411 * w3.gen() );
	}

	/**
	 * Function calculates band-limited windowed sinc function-based filter
	 * kernel.
	 *
	 * @param[out] op Output buffer, length = KernelLen.
	 * @param wfunc Windowing function to use.
	 */

	template< class T >
	void generateBand( T* op,
		CWindowFunc wfunc = &CDSPSincFilterGen :: calcWindowBlackman )
	{
		int t = -fl2;

		while( t < 0 )
		{
			*op = ( f2.gen() - f1.gen() ) * ( *this.*wfunc )() / t / M_PI;
			op++;
			t++;
		}

		f1.gen();
		f2.gen();
		*op = ( Freq2 - Freq1 ) * ( *this.*wfunc )() / M_PI;

		while( t < fl2 )
		{
			t++;
			op++;
			*op = ( f2.gen() - f1.gen() ) * ( *this.*wfunc )() / t / M_PI;
		}
	}

	/**
	 * Function calculates band-limited windowed sinc function-based filter
	 * kernel, with the window raised to the specified power.
	 *
	 * @param[out] op Output buffer, length = KernelLen.
	 * @param p Power factor.
	 * @param wfunc Windowing function to use.
	 */

	template< class T >
	void generateBandPow( T* op, const double p,
		CWindowFunc wfunc = &CDSPSincFilterGen :: calcWindowBlackman )
	{
		R8BASSERT( p > 0.0 );

		int t = -fl2;

		while( t < 0 )
		{
			*op = ( f2.gen() - f1.gen() ) * pows( ( *this.*wfunc )(), p ) /
				t / M_PI;

			op++;
			t++;
		}

		f1.gen();
		f2.gen();
		*op = ( Freq2 - Freq1 ) * pows( ( *this.*wfunc )(), p ) / M_PI;

		while( t < fl2 )
		{
			t++;
			op++;
			*op = ( f2.gen() - f1.gen() ) * pows( ( *this.*wfunc )(), p ) /
				t / M_PI;
		}
	}

	/**
	 * Function calculates windowed Hilbert transformer filter kernel.
	 *
	 * @param[out] op Output buffer, length = KernelLen.
	 * @param wfunc Windowing function to use.
	 */

	template< class T >
	void generateHilbert( T* op,
		CWindowFunc wfunc = &CDSPSincFilterGen :: calcWindowBlackman )
	{
		static const double fvalues[ 2 ] = { 0.0, -2.0 };
		T* op2 = op + fl2 + fl2;
		int t = -fl2;

		while( t < 0 )
		{
			const double v = fvalues[ t & 1 ] * ( *this.*wfunc )() / t / M_PI;
			*op = v;
			op++;
			*op2 = -v;
			op2--;
			t++;
		}

		*op = 0.0;
	}

	/**
	 * Function calculates windowed fractional delay filter kernel.
	 *
	 * @param[out] op Output buffer, length = KernelLen.
	 * @param opinc Output buffer increment, in "op" elements.
	 * @param wfunc Windowing function to use.
	 */

	template< class T >
	void generateFrac( T* op, const int opinc = 1,
		CWindowFunc wfunc = &CDSPSincFilterGen :: calcWindowBlackman )
	{
		R8BASSERT( opinc > 0 );

		double f[ 2 ];
		f[ 0 ] = sin( FracDelay * M_PI );
		f[ 1 ] = -f[ 0 ];

		int t = -fl2;

		if( t + FracDelay < -Len2 )
		{
			( *this.*wfunc )();
			*op = 0.0;
			op += opinc;
			t++;
		}

		while( t < -1 )
		{
			*op = f[ t & 1 ] * ( *this.*wfunc )() / ( t + FracDelay ) / M_PI;
			op += opinc;
			t++;
		}

		const double m = 2e-13 * Len2;
		double ut = t + FracDelay;

		if( fabs( ut ) < m )
		{
			*op = ( *this.*wfunc )();
		}
		else
		{
			*op = f[ t & 1 ] * ( *this.*wfunc )() / ut / M_PI;
		}

		op += opinc;
		t++;
		ut = t + FracDelay;

		if( fabs( ut ) < m )
		{
			*op = ( *this.*wfunc )();
		}
		else
		{
			*op = f[ t & 1 ] * ( *this.*wfunc )() / ut / M_PI;
		}

		const int mt = fl2 - 2;

		while( t < mt )
		{
			op += opinc;
			t++;
			*op = f[ t & 1 ] * ( *this.*wfunc )() / ( t + FracDelay ) / M_PI;
		}

		op += opinc;
		t++;
		ut = t + FracDelay;
		*op = ( ut > Len2 ? 0.0 :
			f[ t & 1 ] * ( *this.*wfunc )() / ut ) / M_PI;
	}

private:
	CSineGen f1; ///< Sine function 1.
		///<
	CSineGen f2; ///< Sine function 2.
		///<
	CSineGen w1; ///< Cosine wave 1 for windowing function.
		///<
	CSineGen w2; ///< Cosine wave 2 for windowing function.
		///<
	CSineGen w3; ///< Cosine wave 3 for windowing function.
		///<
	CSineGen w4; ///< Cosine wave 4 for windowing function. Initialized and
		///< used only by the "Vaneev" windowing function.
		///<
};

} // namespace r8b

#endif // R8B_CDSPSINCFILTERGEN_INCLUDED
