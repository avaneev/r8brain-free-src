/**
 * @file r8bsrc.cpp
 *
 * @brief The "r8bsrc.dll" source file.
 *
 * This source file contains implementation of the DLL functions defined in
 * the "r8bsrc.h" file. This source file was designed for a proprietary DLL
 * compilation process (includes //$ definitions). The compiled DLL and LIB
 * files of the latest version can be downloaded at project's home page:
 * https://code.google.com/p/r8brain-free-src/
 */

//$ def "dll"
//$ bin "*win32|Win32/r8bsrc.dll"
//$ bin "*win64|Win64/r8bsrc.dll"
//$ export "r8b_create"
//$ export "r8b_delete"
//$ export "r8b_get_latency"
//$ export "r8b_clear"
//$ export "r8b_process"

#include "r8bsrc.h"
#include "../CDSPResampler.h"
using namespace r8b;

/**
 * A "wrapper" virtual base class for resamplers of varying resolutions.
 */

class CR8BResamplerBase
{
	R8BNOCTOR( CR8BResamplerBase );

public:
	virtual ~CR8BResamplerBase()
	{
	}

	/**
	 * @return The number of samples that should be passed to *this object
	 * before the actual output starts.
	 */

	virtual int getInLenBeforeOutStart() const = 0;

	/**
	 * Function clears (resets) the state of *this object and returns it to
	 * the state after construction. All input data accumulated in the
	 * internal buffer so far will be discarded. See r8b::CDSPResampler for
	 * more details.
	 */

	virtual void clear() = 0;

	/**
	 * Function performs sample rate conversion. See r8b::CDSPResampler for
	 * more details.
	 */

	virtual int process( double* const ip0, int l, double*& op0 ) = 0;

protected:
	CR8BResamplerBase()
	{
	}
};

/**
 * A templated class that holds a resampler object of the specified class.
 *
 * @param CResampClass The resampler class.
 */

template< class CResampClass >
class CR8BResamplerTpl : public CR8BResamplerBase
{
public:
	CR8BResamplerTpl( CResampClass* const aResamp )
		: Resamp( aResamp )
	{
	}

	virtual ~CR8BResamplerTpl()
	{
		delete Resamp;
	}

	virtual int getInLenBeforeOutStart() const
	{
		Resamp -> getInLenBeforeOutStart();
	}

	virtual void clear()
	{
		Resamp -> clear();
	}

	virtual int process( double* const ip0, int l, double*& op0 )
	{
		return( Resamp -> process( ip0, l, op0 ));
	}

private:
	CResampClass* Resamp;
};

extern "C" {

CR8BResampler _cdecl r8b_create( const double SrcSampleRate,
	const double DstSampleRate, const int MaxInLen,
	const double ReqTransBand, const ER8BResamplerRes Res )
{
	if( Res == r8brr16 )
	{
		return( new CR8BResamplerTpl< CDSPResampler16 >( new CDSPResampler16(
			SrcSampleRate, DstSampleRate, MaxInLen, ReqTransBand )));
	}
	else
	if( Res == r8brr16IR )
	{
		return( new CR8BResamplerTpl< CDSPResampler16IR >(
			new CDSPResampler16IR( SrcSampleRate, DstSampleRate, MaxInLen,
			ReqTransBand )));
	}

	return( new CR8BResamplerTpl< CDSPResampler24 >( new CDSPResampler24(
		SrcSampleRate, DstSampleRate, MaxInLen, ReqTransBand )));
}

void _cdecl r8b_delete( CR8BResampler const rs )
{
	delete (CR8BResamplerBase*) rs;
}

int _cdecl r8b_get_latency( CR8BResampler const rs )
{
	return(( (CR8BResamplerBase*) rs ) -> getInLenBeforeOutStart() );
}

void _cdecl r8b_clear( CR8BResampler const rs )
{
	( (CR8BResamplerBase*) rs ) -> clear();
}

int _cdecl r8b_process( CR8BResampler const rs, double* const ip0, int l,
	double*& op0 )
{
	return(( (CR8BResamplerBase*) rs ) -> process( ip0, l, op0 ));
}

} // extern "C"
