/**
 * @file r8bsrc.cpp
 *
 * @brief The "r8bsrc.dll" source file.
 *
 * This source file contains implementation of the DLL functions defined in
 * the "r8bsrc.h" file. The compiled DLL and LIB files of the latest version
 * can be downloaded at project's home page:
 * https://github.com/avaneev/r8brain-free-src
 *
 * Define a project-wide empty R8BSRC_DECL macro to create a non-DLL object
 * file.
 *
 * r8brain-free-src Copyright (c) 2013-2022 Aleksey Vaneev
 * See the "LICENSE" file for license.
 */

//$ def "dll"
//$ bin "*win32|Win32/r8bsrc.dll"
//$ bin "*win64|Win64/r8bsrc.dll"

#ifndef R8BSRC_DECL
	#define R8BSRC_DECL __declspec( dllexport )
#endif // R8BSRC_DECL
#include "r8bsrc.h"

#include "../CDSPResampler.h"
using namespace r8b;

#if defined( __INTEL_COMPILER )

// Auto-dispatch init "hack" in Intel C++ Compiler to make AVX2 auto-dispatch
// code work on AMD processors.

extern "C" {
	extern long long __intel_cpu_feature_indicator; // CPU feature bits.
	extern long long __intel_cpu_feature_indicator_x; // CPU feature bits.

	void __intel_cpu_features_init_x(); // checks CPU features without
		// discriminating by CPU brand.
}

#endif // defined( __INTEL_COMPILER )

class CAutoDispatchInit
{
public:
	CAutoDispatchInit()
	{
		#if defined( __INTEL_COMPILER )
		__intel_cpu_feature_indicator = 0;
		__intel_cpu_feature_indicator_x = 0;
		__intel_cpu_features_init_x();
		__intel_cpu_feature_indicator = __intel_cpu_feature_indicator_x;
		#endif // defined( __INTEL_COMPILER )

		#if R8B_IPP
		ippInit();
		#endif // R8B_IPP
	}
};

CAutoDispatchInit AutoDispatchInit;

extern "C" {

R8BSRC_DECL CR8BResampler r8b_create( const double SrcSampleRate,
	const double DstSampleRate, const int MaxInLen, const double ReqTransBand,
	const ER8BResamplerRes Res )
{
	if( Res == r8brr16 )
	{
		return( new CDSPResampler16( SrcSampleRate, DstSampleRate, MaxInLen,
			ReqTransBand ));
	}

	if( Res == r8brr16IR )
	{
		return( new CDSPResampler16IR( SrcSampleRate, DstSampleRate, MaxInLen,
			ReqTransBand ));
	}

	return( new CDSPResampler24( SrcSampleRate, DstSampleRate, MaxInLen,
		ReqTransBand ));
}

R8BSRC_DECL void r8b_delete( CR8BResampler const rs )
{
	delete (CDSPProcessor*) rs;
}

R8BSRC_DECL int r8b_inlen( CR8BResampler const rs, const int ReqOutSamples )
{
	if( ReqOutSamples < 1 )
	{
		return( 0 );
	}

	return(( (CDSPProcessor*) rs ) -> getInLenBeforeOutPos(
		ReqOutSamples - 1 ));
}

R8BSRC_DECL void r8b_clear( CR8BResampler const rs )
{
	( (CDSPProcessor*) rs ) -> clear();
}

R8BSRC_DECL int r8b_process( CR8BResampler const rs, double* const ip0,
	const int l, double*& op0 )
{
	return(( (CDSPProcessor*) rs ) -> process( ip0, l, op0 ));
}

} // extern "C"
