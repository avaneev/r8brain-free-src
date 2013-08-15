/**
 * \file r8b.cpp
 * \brief C++ file that should be compiled and included into your application.
 *
 * This is a single library file that should be compiled and included into the
 * project that uses the "r8brain-free-src" sample rate converter. This file
 * defines several global static objects used by the library.
 *
 * r8brain-free-src Copyright (c) 2013 Aleksey Vaneev
 * See the "License.txt" file for license.
 */

#include "CDSPRealFFT.h"

namespace r8b {

CSyncSpinLock CDSPRealFFTKeeper :: StateSync;
CDSPRealFFT :: CObjKeeper CDSPRealFFTKeeper :: FFTObjects[ 31 ];

} // namespace r8b
