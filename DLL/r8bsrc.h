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
 * r8brain-free-src Copyright (c) 2013-2019 Aleksey Vaneev
 * See the "License.txt" file for license.
 */

#ifndef R8BSRC_INCLUDED
#define R8BSRC_INCLUDED

#ifdef R8B_DLL
#if defined (__GNUC__) && ((__GNUC__ >= 4) || ((__GNUC__ == 3) && (__GNUC_MINOR__ >= 1)))
#define R8BRAIN_API __attribute__((visibility ("default")))
#else
#ifdef R8B_DLL_EXPORT
#define R8BRAIN_API __declspec(dllexport)
#else
#define R8BRAIN_API __declspec(dllimport)
#endif
#endif
#else
#define R8BRAIN_API _cdecl
#endif

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
 * @param SrcSampleRate Source signal sample rate. Both sample rates can
 * be specified as a ratio, e.g. SrcSampleRate = 1.0, DstSampleRate = 2.0.
 * @param DstSampleRate Destination signal sample rate.
 * @param MaxInLen The maximal planned length of the input buffer (in samples)
 * that will be passed to the resampler. The resampler relies on this value as
 * it allocates intermediate buffers. Input buffers longer than this value
 * should never be supplied to the resampler. Note that the resampler may use
 * the input buffer itself for intermediate sample data storage.
 * @param Res Resampler's required resolution.
 */

CR8BResampler R8BRAIN_API r8b_create( const double SrcSampleRate,
	const double DstSampleRate, const int MaxInLen,
	const double ReqTransBand, const ER8BResamplerRes Res );

/**
 * Function deletes a resampler previously created via the r8b_create()
 * function.
 *
 * @param rs Resampler object to delete.
 */

void R8BRAIN_API r8b_delete( CR8BResampler const rs );

/**
 * Function clears (resets) the state of the resampler object and returns it
 * to the state after construction. All input data accumulated in the
 * internal buffer of this resampler object so far will be discarded.
 *
 * @param rs Resampler object to clear.
 */

void R8BRAIN_API r8b_clear( CR8BResampler const rs );

/**
 * Function performs sample rate conversion.
 *
 * If the source and destination sample rates are equal, the resampler will do
 * nothing and will simply return the input buffer unchanged.
 *
 * You do not need to allocate an intermediate output buffer for use with this
 * function. If required, the resampler will allocate a suitable intermediate
 * output buffer itself.
 *
 * @param rs Resampler object that performs processing.
 * @param ip0 Input buffer. This buffer may be used as output buffer by this
 * function.
 * @param l The number of samples available in the input buffer.
 * @param[out] op0 This variable receives the pointer to the resampled data.
 * This pointer may point to the address within the "ip0" input buffer, or to
 * *this object's internal buffer. In real-time applications it is suggested
 * to pass this pointer to the next output audio block and consume any data
 * left from the previous output audio block first before calling the
 * r8b_process() function again. The buffer pointed to by the "op0" on return
 * may be owned by the resampler, so it should not be freed by the caller.
 * @return The number of samples available in the "op0" output buffer. If the
 * data from the output buffer "op0" is going to be written to a bigger output
 * buffer, it is suggested to check the returned number of samples so that no
 * overflow of the bigger output buffer happens.
 */

int R8BRAIN_API r8b_process( CR8BResampler const rs, double* const ip0, int l,
	double*& op0 );

} // extern "C"

#endif // R8BSRC_INCLUDED
