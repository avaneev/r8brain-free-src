/**
 * @file r8bsrc.cpp
 *
 * @brief The "r8bsrc.dll" source file.
 *
 * This source file contains implementation of the DLL functions defined in
 * the "r8bsrc.h" file. This source file was designed for a proprietary DLL
 * compilation process (includes //$ definitions). The compiled DLL and LIB
 * files of the latest version can be downloaded at project's home page:
 * https://github.com/avaneev/r8brain-free-src
 *
 * r8brain-free-src Copyright (c) 2013-2019 Aleksey Vaneev
 * See the "License.txt" file for license.
 */

//$ def "dll"
//$ bin "*win32|Win32/r8bsrc.dll"
//$ bin "*win64|Win64/r8bsrc.dll"
//$ export "r8b_create"
//$ export "r8b_delete"
//$ export "r8b_clear"
//$ export "r8b_process"

#include "r8bsrc.h"
#include "../CDSPResampler.h"
using namespace r8b;

#if defined( __INTEL_COMPILER )

// Auto-dispatch init "hack" in Intel C++ Compiler to make AVX2 auto-dispatch
// code work on AMD processors.

extern "C" {
	extern long long __intel_cpu_feature_indicator; // CPU feature bits
	extern long long __intel_cpu_feature_indicator_x; // CPU feature bits
	void __intel_cpu_features_init_x(); // checks CPU features without discriminating by CPU brand
}

class CAutoDispatchInit
{
public:
	CAutoDispatchInit()
	{
		__intel_cpu_feature_indicator = 0;
		__intel_cpu_feature_indicator_x = 0;
		__intel_cpu_features_init_x();
		__intel_cpu_feature_indicator = __intel_cpu_feature_indicator_x;
	}
};

CAutoDispatchInit AutoDispatchInit;

#endif // defined( __INTEL_COMPILER )

extern "C" {

CR8BResampler _cdecl r8b_create( const double SrcSampleRate,
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

void _cdecl r8b_delete( CR8BResampler const rs )
{
	delete (CDSPProcessor*) rs;
}

void _cdecl r8b_clear( CR8BResampler const rs )
{
	( (CDSPProcessor*) rs ) -> clear();
}

int _cdecl r8b_process( CR8BResampler const rs, double* const ip0, int l,
	double*& op0 )
{
	return(( (CDSPProcessor*) rs ) -> process( ip0, l, op0 ));
}

} // extern "C"
