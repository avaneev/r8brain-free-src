//$ nocpp

/**
 * @file r8bsrc.h
 *
 * @brief Inclusion file for use with the "r8bsrc.dll".
 *
 * This is the inclusion file for the "r8bsrc.dll" dynamic link library
 * (the "r8bsrc.lib" library should be included into the project). This DLL
 * is designed to be used on Windows, on a processor with SSE2 support.
 * On non-Windows systems it is preferrable to use the C++ library directly.
 *
 * Before using this DLL library please read the description of the
 * r8b::CDSPResampler class and its functions.
 *
 * Note that the "int" and "enum" types have 32-bit size on both 32-bit and
 * 64-bit systems. Pointer types, including the CR8BResampler type, have
 * 32-bit size on 32-bit system and 64-bit size on 64-bit system.
 *
 * r8brain-free-src Copyright (c) 2013 Aleksey Vaneev
 * See the "License.txt" file for license.
 */

#ifndef R8BSRC_INCLUDED
#define R8BSRC_INCLUDED

/**
 * Resampler object handle.
 */

typedef void* CR8BResampler;

/**
 * Possible resampler object resolutions.
 */

enum ER8BResamplerRes
{
	r8brr16 = 0, ///< 16-bit precision resampler.
		///<
	r8brr16IR = 1, ///< 16-bit precision resampler for impulse responses.
		///<
	r8brr24 = 2 ///< 24-bit precision resampler (including 32-bit floating
		///< point).
		///<
};

extern "C" {

/**
 * Function creates a new linear-phase resampler object.
 *
 * @param SrcSampleRate Source signal sample rate.
 * @param DstSampleRate Destination signal sample rate.
 * @param MaxInLen The maximal planned length of the input buffer (in samples)
 * that will be passed to the resampler. The resampler relies on this value as
 * it allocates intermediate buffers. Input buffers longer than this value
 * should never be supplied to the resampler. Note that the resampler may use
 * the input buffer itself for intermediate sample data storage.
 * @param Res Resampler's required resolution.
 */

CR8BResampler _cdecl r8b_create( const double SrcSampleRate,
	const double DstSampleRate, const int MaxInLen,
	const double ReqTransBand, const ER8BResamplerRes Res );

/**
 * Function deletes a resampler previously created via the r8b_create()
 * function.
 *
 * @param rs Resampler object to delete.
 */

void _cdecl r8b_delete( CR8BResampler const rs );

/**
 * Function returns the number of samples that should be passed to the
 * resampler object before the actual output starts.
 *
 * @param rs Resampler object.
 */

int _cdecl r8b_get_latency( CR8BResampler const rs );

/**
 * Function clears (resets) the state of the resampler object and returns it
 * to the state after construction. All input data accumulated in the
 * internal buffer of this resampler object so far will be discarded.
 *
 * @param rs Resampler object to clear.
 */

void _cdecl r8b_clear( CR8BResampler const rs );

/**
 * Function performs sample rate conversion.
 *
 * If the source and destination sample rates are equal, the resampler will do
 * nothing and will simply return the input buffer unchanged.
 *
 * @param ip0 Input buffer. This buffer may be used as output buffer by this
 * function.
 * @param l The number of samples available in the input buffer.
 * @param[out] op0 This variable receives the pointer to the resampled data.
 * This pointer may point to the address within the "ip0" input buffer, or to
 * *this object's internal buffer. In real-time applications it is suggested
 * to pass this pointer to the next output audio block and consume any data
 * left from the previous output audio block first before calling the
 * process() function again.
 * @return The number of samples available in the "op0" output buffer.
 */

int _cdecl r8b_process( CR8BResampler const rs, double* const ip0, int l,
	double*& op0 );

} // extern "C"

#endif // R8BSRC_INCLUDED
