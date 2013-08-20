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
	double Len2; ///< Required half kernel length in samples (can be a
		///< fractional value). Final physical kernel length will be provided
		///< in the KernelLen variable. Len2 should be >= 2.
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
	int fl2; ///< Internal "half kernel size" value. This value can be used
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
	 * Function calculates band-limited windowed sinc function.
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
	 * Function calculates band-limited windowed sinc function, with the
	 * window raised to the specified power.
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
	 * Function calculates windowed Hilbert transformer filter.
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
	 * Function calculates windowed fractional delay filter.
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
		f[ 0 ] = sin( FracDelay * M_PI ) / M_PI;
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
			*op = f[ t & 1 ] * ( *this.*wfunc )() / ( t + FracDelay );
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
			*op = f[ t & 1 ] * ( *this.*wfunc )() / ut;
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
			*op = f[ t & 1 ] * ( *this.*wfunc )() / ut;
		}

		const int mt = fl2 - 2;

		while( t < mt )
		{
			op += opinc;
			t++;
			*op = f[ t & 1 ] * ( *this.*wfunc )() / ( t + FracDelay );
		}

		op += opinc;
		t++;
		ut = t + FracDelay;
		*op = ( ut > Len2 ? 0.0 : f[ t & 1 ] * ( *this.*wfunc )() / ut );
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
};

} // namespace r8b

#endif // R8B_CDSPSINCFILTERGEN_INCLUDED
