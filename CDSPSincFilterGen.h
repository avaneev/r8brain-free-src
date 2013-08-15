//$ nobt
//$ nocpp

/**
 * \file CDSPSincFilterGen.h
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
		/// fractional value). Final physical kernel length will be provided
		/// in the KernelLen variable. Len2 should be >= 2.
	double Freq1; ///< Required corner circular frequency 1 [0; pi]. Used only
		/// in the generateBand() function.
	double Freq2; ///< Required corner circular frequency 2 [0; pi]. Used only
		/// in the generateBand() function. The range [Freq1; Freq2] defines
		/// a pass band for the generateBand() function.
	double FracDelay; ///< Fractional delay in the range [0; 1], used only
		/// by the generateFrac() function. Note that the FracDelay parameter
		/// is actually inversed. At 0.0 value it produces 1 sample delay
		/// (with latency equal to fl2), at 1.0 value it produces 0 sample
		/// delay (with latency equal to fl2 - 1).
	int KernelLen; ///< Resulting length of the filter kernel, this variable
		/// is set after the call to one of the "init" functions.
	int fl2; ///< Internal "half kernel size" value. This value can be used
		/// as filter's latency in samples (taps), this variable is set after
		/// the call to one of the init functions.

	typedef double( CDSPSincFilterGen :: *CWindowFunc )(); ///< Windowing
		/// function pointer type.

	/**
	 * Function initializes *this structure for generation of band-limited
	 * sinc filter kernel. The generateBand() or generateBandPow() functions
	 * should be used to calculate the filter.
	 */

	void initBand()
	{
		fl2 = (int) floor( Len2 );
		KernelLen = fl2 + fl2 + 1;
		f1.init( Freq1, 2.0, -Freq1 * fl2 );
		f2.init( Freq2, 2.0, -Freq2 * fl2 );

		const double wstep = M_PI / Len2;
		w.init( wstep, 2.0, -wstep * fl2 + M_PI * 0.5 );

		const double wstep2 = wstep + wstep;
		w2.init( wstep2, 2.0, -wstep2 * fl2 + M_PI * 0.5 );

		const double wstep3 = wstep2 + wstep;
		w3.init( wstep3, 2.0, -wstep3 * fl2 + M_PI * 0.5 );
	}

	/**
	 * Function initializes *this structure for Hilbert transformation filter
	 * calculation. Freq1 and Freq2 variables are not used.
	 * The generateHilbert() function should be used to calculate the filter.
	 */

	void initHilbert()
	{
		fl2 = (int) floor( Len2 );
		KernelLen = fl2 + fl2 + 1;
		f2.init( M_PI, 2.0, ( (double) -fl2 + 0.5 ) * M_PI );

		const double wstep = M_PI / Len2;
		w.init( wstep, 2.0, -wstep * fl2 + M_PI * 0.5 );

		const double wstep2 = wstep + wstep;
		w2.init( wstep2, 2.0, -wstep2 * fl2 + M_PI * 0.5 );

		const double wstep3 = wstep2 + wstep;
		w3.init( wstep3, 2.0, -wstep3 * fl2 + M_PI * 0.5 );
	}

	/**
	 * Function initializes *this structure for generation of full-bandwidth
	 * fractional delay sinc filter kernel. Freq1 and Freq2 variables are not
	 * used. The generateFrac() function should be used to calculate the
	 * filter.
	 */

	void initFrac()
	{
		fl2 = (int) ceil( Len2 );
		KernelLen = fl2 + fl2;
		const double fl2b = fl2 - FracDelay;
		f2.init( M_PI, 2.0, -M_PI * fl2b );

		const double wstep = M_PI / Len2;
		w.init( wstep, 2.0, -wstep * fl2b + M_PI * 0.5 );

		const double wstep2 = wstep + wstep;
		w2.init( wstep2, 2.0, -wstep2 * fl2b + M_PI * 0.5 );

		const double wstep3 = wstep2 + wstep;
		w3.init( wstep3, 2.0, -wstep3 * fl2b + M_PI * 0.5 );
	}

	/**
	 * @return The next "Hann" windowing function coefficient.
	 */

	double calcWindowHann()
	{
		return( 0.5 + 0.5 * w.gen() );
	}

	/**
	 * @return The next "Hamming" windowing function coefficient.
	 */

	double calcWindowHamming()
	{
		return( 0.54 + 0.46 * w.gen() );
	}

	/**
	 * @return The next "Blackman" windowing function coefficient.
	 */

	double calcWindowBlackman()
	{
		return( 0.42 + 0.5 * w.gen() + 0.08 * w2.gen() );
	}

	/**
	 * @return The next "Nuttall" windowing function coefficient.
	 */

	double calcWindowNuttall()
	{
		return( 0.355768 + 0.487396 * w.gen() + 0.144232 * w2.gen() +
			0.012604 * w3.gen() );
	}

	/**
	 * Function calculates band-limited windowed sinc function.
	 *
	 * @param[out] op Output buffer, length = KernelLen.
	 * @param wfunc Windowing function to use.
	 */

	void generateBand( double* op,
		CWindowFunc wfunc = &CDSPSincFilterGen :: calcWindowBlackman )
	{
		int t = -fl2;

		while( t < 0 )
		{
			*op = ( f2.gen() - f1.gen() ) * ( *this.*wfunc )() / M_PI / t;
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
			*op = ( f2.gen() - f1.gen() ) * ( *this.*wfunc )() / M_PI / t;
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

	void generateBandPow( double* op, const double p,
		CWindowFunc wfunc = &CDSPSincFilterGen :: calcWindowBlackman )
	{
		int t = -fl2;

		while( t < 0 )
		{
			*op = ( f2.gen() - f1.gen() ) * pows( ( *this.*wfunc )(), p ) /
				M_PI / t;

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
				M_PI / t;
		}
	}

	/**
	 * Function calculates windowed Hilbert transformer filter.
	 *
	 * @param[out] op Output buffer, length = KernelLen.
	 * @param wfunc Windowing function to use.
	 */

	void generateHilbert( double* op,
		CWindowFunc wfunc = &CDSPSincFilterGen :: calcWindowBlackman )
	{
		double t = -M_PI * fl2;
		int l = fl2;

		while( l > 0 )
		{
			*op = ( f2.gen() - 1.0 ) * ( *this.*wfunc )() / t;
			op++;
			t += M_PI;
			l--;
		}

		*op = 0.0;
		double* ip = op;
		l = fl2;

		while( l > 0 )
		{
			op++;
			ip--;
			*op = -*ip;
			l--;
		}
	}

	/**
	 * Function calculates windowed fractional delay filter.
	 *
	 * @param[out] op Output buffer, length = KernelLen.
	 * @param opinc Output buffer increment, in "op" elements.
	 * @param wfunc Windowing function to use.
	 */

	void generateFrac( double* op, const int opinc = 1,
		CWindowFunc wfunc = &CDSPSincFilterGen :: calcWindowBlackman )
	{
		int t = -fl2;

		if( t + FracDelay < -Len2 )
		{
			f2.gen();
			( *this.*wfunc )();
			*op = 0.0;
			op += opinc;
			t++;
		}

		while( t < -1 )
		{
			*op = f2.gen() * ( *this.*wfunc )() / M_PI / ( t + FracDelay );
			op += opinc;
			t++;
		}

		const double m = 2e-13 * Len2;
		double ut = t + FracDelay;

		if( fabs( ut ) < m )
		{
			f2.gen();
			*op = ( *this.*wfunc )();
		}
		else
		{
			*op = f2.gen() * ( *this.*wfunc )() / M_PI / ut;
		}

		op += opinc;
		t++;
		ut = t + FracDelay;

		if( fabs( ut ) < m )
		{
			f2.gen();
			*op = ( *this.*wfunc )();
		}
		else
		{
			*op = f2.gen() * ( *this.*wfunc )() / M_PI / ut;
		}

		const int mt = fl2 - 2;

		while( t < mt )
		{
			op += opinc;
			t++;
			*op = f2.gen() * ( *this.*wfunc )() / M_PI / ( t + FracDelay );
		}

		op += opinc;
		t++;
		ut = t + FracDelay;
		*op = ( ut > Len2 ? 0.0 : f2.gen() * ( *this.*wfunc )() / M_PI / ut );
	}

private:
	CSineGen f1; ///< Sine function 1.
	CSineGen f2; ///< Sine function 2.
	CSineGen w; ///< Windowing function.
	CSineGen w2; ///< Windowing function 2.
	CSineGen w3; ///< Windowing function 3.
};

} // namespace r8b

#endif // R8B_CDSPSINCFILTERGEN_INCLUDED
