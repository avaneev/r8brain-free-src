/* Copyright (c) 2013  Julien Pommier ( pommier@modartt.com )

   Based on original fortran 77 code from FFTPACKv4 from NETLIB
   (http://www.netlib.org/fftpack), authored by Dr Paul Swarztrauber
   of NCAR, in 1985.

   As confirmed by the NCAR fftpack software curators, the following
   FFTPACKv5 license applies to FFTPACKv4 sources. My changes are
   released under the same terms.

   FFTPACK license:

   http://www.cisl.ucar.edu/css/software/fftpack5/ftpk.html

   Copyright (c) 2004 the University Corporation for Atmospheric
   Research ("UCAR"). All rights reserved. Developed by NCAR's
   Computational and Information Systems Laboratory, UCAR,
   www.cisl.ucar.edu.

   Redistribution and use of the Software in source and binary forms,
   with or without modification, is permitted provided that the
   following conditions are met:

   - Neither the names of NCAR's Computational and Information Systems
   Laboratory, the University Corporation for Atmospheric Research,
   nor the names of its sponsors or contributors may be used to
   endorse or promote products derived from this Software without
   specific prior written permission.  

   - Redistributions of source code must retain the above copyright
   notices, this list of conditions, and the disclaimer below.

   - Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions, and the disclaimer below in the
   documentation and/or other materials provided with the
   distribution.

   THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE CONTRIBUTORS OR COPYRIGHT
   HOLDERS BE LIABLE FOR ANY CLAIM, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS WITH THE
   SOFTWARE.


   PFFFT : a Pretty Fast FFT.

   This file is largerly based on the original FFTPACK implementation, modified in
   order to take advantage of SIMD instructions of modern CPUs.
*/

/*
  ChangeLog: 
  - 2011/10/02, version 1: This is the very first release of this file.
*/

/*
   NOTE: This file is adapted from Julien Pommier's original PFFFT,
   which works on 32 bit floating point precision using SSE instructions,
   to work with 64 bit floating point precision using AVX instructions.
   Author: Dario Mambro @ https://github.com/unevens/pffft
*/

/* detect compiler flavour */
#if defined(_MSC_VER)
#  define COMPILER_MSVC
#elif defined(__GNUC__)
#  define COMPILER_GCC
#endif

#ifdef COMPILER_MSVC
#  define _USE_MATH_DEFINES
#include <malloc.h>
#else
#include <alloca.h>
#endif

#include "pffft-double.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#if defined(COMPILER_GCC)
#  define ALWAYS_INLINE(return_type) inline return_type __attribute__ ((always_inline))
#  define NEVER_INLINE(return_type) return_type __attribute__ ((noinline))
#  define RESTRICT __restrict
#  define VLA_ARRAY_ON_STACK(type__, varname__, size__) type__ varname__[size__];
#elif defined(COMPILER_MSVC)
#  define ALWAYS_INLINE(return_type) __forceinline return_type
#  define NEVER_INLINE(return_type) __declspec(noinline) return_type
#  define RESTRICT __restrict
#  define VLA_ARRAY_ON_STACK(type__, varname__, size__) type__ *varname__ = (type__*)_alloca(size__ * sizeof(type__))
#endif


/* 
   vector support macros: the rest of the code is independant of
   AVX -- adding support for other platforms with 4-element
   vectors should be limited to these macros 
*/


// define PFFFTD_SIMD_DISABLE if you want to use scalar code instead of simd code
//#define PFFFTD_SIMD_DISABLE

/*
  AVX support macros
*/
#if !defined(PFFFTD_SIMD_DISABLE) && defined(__AVX__)

#include <immintrin.h>
typedef __m256d v4sd;
#  define SIMD_SZ 4 // 4 doubles by simd vector 
#  define VZERO() _mm256_setzero_pd()
#  define VMUL(a,b) _mm256_mul_pd(a,b)
#  define VADD(a,b) _mm256_add_pd(a,b)
#  define VMADD(a,b,c) _mm256_add_pd(_mm256_mul_pd(a,b), c)
#  define VSUB(a,b) _mm256_sub_pd(a,b)
#  define LD_PS1(p) _mm256_set1_pd(p)
/* INTERLEAVE2 (in1, in2, out1, out2) pseudo code:
out1 = [ in1[0], in2[0], in1[1], in2[1] ]
out2 = [ in1[2], in2[2], in1[3], in2[3] ]
*/
#  define INTERLEAVE2(in1, in2, out1, out2) {							\
	__m128d low1__ = _mm256_castpd256_pd128(in1);						\
	__m128d low2__ = _mm256_castpd256_pd128(in2);						\
	__m128d high1__ = _mm256_extractf128_pd(in1, 1);					\
	__m128d high2__ = _mm256_extractf128_pd(in2, 1);					\
	__m256d tmp__ = _mm256_insertf128_pd(								\
		_mm256_castpd128_pd256(_mm_shuffle_pd(low1__, low2__, 0)),		\
		_mm_shuffle_pd(low1__, low2__, 3),								\
		1);																\
	out2 = _mm256_insertf128_pd(										\
		_mm256_castpd128_pd256(_mm_shuffle_pd(high1__, high2__, 0)),	\
		_mm_shuffle_pd(high1__, high2__, 3),							\
		1);																\
	out1 = tmp__;														\
}
/*UNINTERLEAVE2(in1, in2, out1, out2) pseudo code:
out1 = [ in1[0], in1[2], in2[0], in2[2] ]
out2 = [ in1[1], in1[3], in2[1], in2[3] ]
*/
#  define UNINTERLEAVE2(in1, in2, out1, out2) {							\
	__m128d low1__ = _mm256_castpd256_pd128(in1);						\
	__m128d low2__ = _mm256_castpd256_pd128(in2);						\
	__m128d high1__ = _mm256_extractf128_pd(in1, 1);					\
	__m128d high2__ = _mm256_extractf128_pd(in2, 1); 					\
	__m256d tmp__ = _mm256_insertf128_pd(								\
		_mm256_castpd128_pd256(_mm_shuffle_pd(low1__, high1__, 0)),		\
		_mm_shuffle_pd(low2__, high2__, 0),								\
		1);																\
	out2 = _mm256_insertf128_pd(										\
		_mm256_castpd128_pd256(_mm_shuffle_pd(low1__, high1__, 3)),		\
		_mm_shuffle_pd(low2__, high2__, 3),								\
		1);																\
	out1 = tmp__;														\
}
#  define VTRANSPOSE4(row0, row1, row2, row3) {				\
        __m256d tmp3, tmp2, tmp1, tmp0;                     \
                                                            \
        tmp0 = _mm256_shuffle_pd((row0),(row1), 0x0);       \
        tmp2 = _mm256_shuffle_pd((row0),(row1), 0xF);       \
        tmp1 = _mm256_shuffle_pd((row2),(row3), 0x0);       \
        tmp3 = _mm256_shuffle_pd((row2),(row3), 0xF);       \
                                                            \
        (row0) = _mm256_permute2f128_pd(tmp0, tmp1, 0x20);	\
        (row1) = _mm256_permute2f128_pd(tmp2, tmp3, 0x20);  \
        (row2) = _mm256_permute2f128_pd(tmp0, tmp1, 0x31);  \
        (row3) = _mm256_permute2f128_pd(tmp2, tmp3, 0x31);  \
    }
/*VSWAPHL(a, b) pseudo code:
return [ b[0], b[1], a[2], a[3] ]
*/
#  define VSWAPHL(a,b)	\
   _mm256_insertf128_pd(_mm256_castpd128_pd256(_mm256_castpd256_pd128(b)), _mm256_extractf128_pd(a, 1), 1)
#  define VALIGNED(ptr) ((((long)(ptr)) & 0x1F) == 0)
#else
#define PFFFTD_SIMD_DISABLE
#endif

// fallback mode for situations where AVX is not available, use scalar mode instead
#ifdef PFFFTD_SIMD_DISABLE
typedef double v4sd;
#  define SIMD_SZ 1
#  define VZERO() 0.f
#  define VMUL(a,b) ((a)*(b))
#  define VADD(a,b) ((a)+(b))
#  define VMADD(a,b,c) ((a)*(b)+(c))
#  define VSUB(a,b) ((a)-(b))
#  define LD_PS1(p) (p)
#  define VALIGNED(ptr) ((((long)(ptr)) & 0x3) == 0)
#endif

// shortcuts for complex multiplcations
#define VCPLXMUL(ar,ai,br,bi) { v4sd tmp; tmp=VMUL(ar,bi); ar=VMUL(ar,br); ar=VSUB(ar,VMUL(ai,bi)); ai=VMUL(ai,br); ai=VADD(ai,tmp); }
#define VCPLXMULCONJ(ar,ai,br,bi) { v4sd tmp; tmp=VMUL(ar,bi); ar=VMUL(ar,br); ar=VADD(ar,VMUL(ai,bi)); ai=VMUL(ai,br); ai=VSUB(ai,tmp); }
#ifndef SVMUL
// multiply a scalar with a vector
#define SVMUL(f,v) VMUL(LD_PS1(f),v)
#endif

#if !defined(PFFFTD_SIMD_DISABLE)
typedef union v4sd_union {
  v4sd  v;
  double f[4];
} v4sd_union;

#include <string.h>

#define assertv4(v,f0,f1,f2,f3) assert(v.f[0] == (f0) && v.f[1] == (f1) && v.f[2] == (f2) && v.f[3] == (f3))

/* detect bugs with the vector support macros */
void validate_pffftd_simd() {
  double f[16] = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 };
  v4sd_union a0, a1, a2, a3, t, u; 
  memcpy(a0.f, f, 4*sizeof(double));
  memcpy(a1.f, f+4, 4*sizeof(double));
  memcpy(a2.f, f+8, 4*sizeof(double));
  memcpy(a3.f, f+12, 4*sizeof(double));

  t = a0; u = a1; t.v = VZERO();
  printf("VZERO=[%2g %2g %2g %2g]\n", t.f[0], t.f[1], t.f[2], t.f[3]); assertv4(t, 0, 0, 0, 0);
  t.v = VADD(a1.v, a2.v);
  printf("VADD(4:7,8:11)=[%2g %2g %2g %2g]\n", t.f[0], t.f[1], t.f[2], t.f[3]); assertv4(t, 12, 14, 16, 18);
  t.v = VMUL(a1.v, a2.v);
  printf("VMUL(4:7,8:11)=[%2g %2g %2g %2g]\n", t.f[0], t.f[1], t.f[2], t.f[3]); assertv4(t, 32, 45, 60, 77);
  t.v = VMADD(a1.v, a2.v,a0.v);
  printf("VMADD(4:7,8:11,0:3)=[%2g %2g %2g %2g]\n", t.f[0], t.f[1], t.f[2], t.f[3]); assertv4(t, 32, 46, 62, 80);

  INTERLEAVE2(a1.v,a2.v,t.v,u.v);
  printf("INTERLEAVE2(4:7,8:11)=[%2g %2g %2g %2g] [%2g %2g %2g %2g]\n", t.f[0], t.f[1], t.f[2], t.f[3], u.f[0], u.f[1], u.f[2], u.f[3]);
  assertv4(t, 4, 8, 5, 9); assertv4(u, 6, 10, 7, 11);
  UNINTERLEAVE2(a1.v,a2.v,t.v,u.v);
  printf("UNINTERLEAVE2(4:7,8:11)=[%2g %2g %2g %2g] [%2g %2g %2g %2g]\n", t.f[0], t.f[1], t.f[2], t.f[3], u.f[0], u.f[1], u.f[2], u.f[3]);
  assertv4(t, 4, 6, 8, 10); assertv4(u, 5, 7, 9, 11);

  t.v=LD_PS1(f[15]);
  printf("LD_PS1(15)=[%2g %2g %2g %2g]\n", t.f[0], t.f[1], t.f[2], t.f[3]);
  assertv4(t, 15, 15, 15, 15);
  t.v = VSWAPHL(a1.v, a2.v);
  printf("VSWAPHL(4:7,8:11)=[%2g %2g %2g %2g]\n", t.f[0], t.f[1], t.f[2], t.f[3]);
  assertv4(t, 8, 9, 6, 7);
  VTRANSPOSE4(a0.v, a1.v, a2.v, a3.v);
  printf("VTRANSPOSE4(0:3,4:7,8:11,12:15)=[%2g %2g %2g %2g] [%2g %2g %2g %2g] [%2g %2g %2g %2g] [%2g %2g %2g %2g]\n", 
         a0.f[0], a0.f[1], a0.f[2], a0.f[3], a1.f[0], a1.f[1], a1.f[2], a1.f[3], 
         a2.f[0], a2.f[1], a2.f[2], a2.f[3], a3.f[0], a3.f[1], a3.f[2], a3.f[3]); 
  assertv4(a0, 0, 4, 8, 12); assertv4(a1, 1, 5, 9, 13); assertv4(a2, 2, 6, 10, 14); assertv4(a3, 3, 7, 11, 15);
}
#endif //!PFFFTD_SIMD_DISABLE

/* SSE and co like 16-bytes aligned pointers */
#define MALLOC_V4SF_ALIGNMENT 64 // with a 64-byte alignment, we are even aligned on L2 cache lines...
void *pffftd_aligned_malloc(size_t nb_bytes) {
  void *p, *p0 = malloc(nb_bytes + MALLOC_V4SF_ALIGNMENT);
  if (!p0) return (void *) 0;
  p = (void *) (((size_t) p0 + MALLOC_V4SF_ALIGNMENT) & (~((size_t) (MALLOC_V4SF_ALIGNMENT-1))));
  *((void **) p - 1) = p0;
  return p;
}

void pffftd_aligned_free(void *p) {
  if (p) free(*((void **) p - 1));
}

int pffftd_simd_size() { return SIMD_SZ; }

/*
  passf2 and passb2 has been merged here, fsign = -1 for passf2, +1 for passb2
*/
static NEVER_INLINE(void) passf2_ps(int ido, int l1, const v4sd *cc, v4sd *ch, const double *wa1, double fsign) {
  int k, i;
  int l1ido = l1*ido;
  if (ido <= 2) {
    for (k=0; k < l1ido; k += ido, ch += ido, cc+= 2*ido) {
      ch[0]         = VADD(cc[0], cc[ido+0]);
      ch[l1ido]     = VSUB(cc[0], cc[ido+0]);
      ch[1]         = VADD(cc[1], cc[ido+1]);
      ch[l1ido + 1] = VSUB(cc[1], cc[ido+1]);
    }
  } else {
    for (k=0; k < l1ido; k += ido, ch += ido, cc += 2*ido) {
      for (i=0; i<ido-1; i+=2) {
        v4sd tr2 = VSUB(cc[i+0], cc[i+ido+0]);
        v4sd ti2 = VSUB(cc[i+1], cc[i+ido+1]);
        v4sd wr = LD_PS1(wa1[i]), wi = VMUL(LD_PS1(fsign), LD_PS1(wa1[i+1]));
        ch[i]   = VADD(cc[i+0], cc[i+ido+0]);
        ch[i+1] = VADD(cc[i+1], cc[i+ido+1]);
        VCPLXMUL(tr2, ti2, wr, wi);
        ch[i+l1ido]   = tr2;
        ch[i+l1ido+1] = ti2;
      }
    }
  }
}

/*
  passf3 and passb3 has been merged here, fsign = -1 for passf3, +1 for passb3
*/
static NEVER_INLINE(void) passf3_ps(int ido, int l1, const v4sd *cc, v4sd *ch,
                                    const double *wa1, const double *wa2, double fsign) {
  static const double taur = -0.5f;
  double taui = 0.866025403784439f*fsign;
  int i, k;
  v4sd tr2, ti2, cr2, ci2, cr3, ci3, dr2, di2, dr3, di3;
  int l1ido = l1*ido;
  double wr1, wi1, wr2, wi2;
  assert(ido > 2);
  for (k=0; k< l1ido; k += ido, cc+= 3*ido, ch +=ido) {
    for (i=0; i<ido-1; i+=2) {
      tr2 = VADD(cc[i+ido], cc[i+2*ido]);
      cr2 = VADD(cc[i], SVMUL(taur,tr2));
      ch[i]    = VADD(cc[i], tr2);
      ti2 = VADD(cc[i+ido+1], cc[i+2*ido+1]);
      ci2 = VADD(cc[i    +1], SVMUL(taur,ti2));
      ch[i+1]  = VADD(cc[i+1], ti2);
      cr3 = SVMUL(taui, VSUB(cc[i+ido], cc[i+2*ido]));
      ci3 = SVMUL(taui, VSUB(cc[i+ido+1], cc[i+2*ido+1]));
      dr2 = VSUB(cr2, ci3);
      dr3 = VADD(cr2, ci3);
      di2 = VADD(ci2, cr3);
      di3 = VSUB(ci2, cr3);
      wr1=wa1[i], wi1=fsign*wa1[i+1], wr2=wa2[i], wi2=fsign*wa2[i+1]; 
      VCPLXMUL(dr2, di2, LD_PS1(wr1), LD_PS1(wi1));
      ch[i+l1ido] = dr2; 
      ch[i+l1ido + 1] = di2;
      VCPLXMUL(dr3, di3, LD_PS1(wr2), LD_PS1(wi2));
      ch[i+2*l1ido] = dr3;
      ch[i+2*l1ido+1] = di3;
    }
  }
} /* passf3 */

static NEVER_INLINE(void) passf4_ps(int ido, int l1, const v4sd *cc, v4sd *ch,
                                    const double *wa1, const double *wa2, const double *wa3, double fsign) {
  /* isign == -1 for forward transform and +1 for backward transform */

  int i, k;
  v4sd ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;
  int l1ido = l1*ido;
  if (ido == 2) {
    for (k=0; k < l1ido; k += ido, ch += ido, cc += 4*ido) {
      tr1 = VSUB(cc[0], cc[2*ido + 0]);
      tr2 = VADD(cc[0], cc[2*ido + 0]);
      ti1 = VSUB(cc[1], cc[2*ido + 1]);
      ti2 = VADD(cc[1], cc[2*ido + 1]);
      ti4 = VMUL(VSUB(cc[1*ido + 0], cc[3*ido + 0]), LD_PS1(fsign));
      tr4 = VMUL(VSUB(cc[3*ido + 1], cc[1*ido + 1]), LD_PS1(fsign));
      tr3 = VADD(cc[ido + 0], cc[3*ido + 0]);
      ti3 = VADD(cc[ido + 1], cc[3*ido + 1]);

      ch[0*l1ido + 0] = VADD(tr2, tr3);
      ch[0*l1ido + 1] = VADD(ti2, ti3);
      ch[1*l1ido + 0] = VADD(tr1, tr4);
      ch[1*l1ido + 1] = VADD(ti1, ti4);
      ch[2*l1ido + 0] = VSUB(tr2, tr3);
      ch[2*l1ido + 1] = VSUB(ti2, ti3);        
      ch[3*l1ido + 0] = VSUB(tr1, tr4);
      ch[3*l1ido + 1] = VSUB(ti1, ti4);
    }
  } else {
    for (k=0; k < l1ido; k += ido, ch+=ido, cc += 4*ido) {
      for (i=0; i<ido-1; i+=2) {
        double wr1, wi1, wr2, wi2, wr3, wi3;
        tr1 = VSUB(cc[i + 0], cc[i + 2*ido + 0]);
        tr2 = VADD(cc[i + 0], cc[i + 2*ido + 0]);
        ti1 = VSUB(cc[i + 1], cc[i + 2*ido + 1]);
        ti2 = VADD(cc[i + 1], cc[i + 2*ido + 1]);
        tr4 = VMUL(VSUB(cc[i + 3*ido + 1], cc[i + 1*ido + 1]), LD_PS1(fsign));
        ti4 = VMUL(VSUB(cc[i + 1*ido + 0], cc[i + 3*ido + 0]), LD_PS1(fsign));
        tr3 = VADD(cc[i + ido + 0], cc[i + 3*ido + 0]);
        ti3 = VADD(cc[i + ido + 1], cc[i + 3*ido + 1]);

        ch[i] = VADD(tr2, tr3);
        cr3    = VSUB(tr2, tr3);
        ch[i + 1] = VADD(ti2, ti3);
        ci3 = VSUB(ti2, ti3);

        cr2 = VADD(tr1, tr4);
        cr4 = VSUB(tr1, tr4);
        ci2 = VADD(ti1, ti4);
        ci4 = VSUB(ti1, ti4);
        wr1=wa1[i], wi1=fsign*wa1[i+1];
        VCPLXMUL(cr2, ci2, LD_PS1(wr1), LD_PS1(wi1));
        wr2=wa2[i], wi2=fsign*wa2[i+1]; 
        ch[i + l1ido] = cr2;
        ch[i + l1ido + 1] = ci2;

        VCPLXMUL(cr3, ci3, LD_PS1(wr2), LD_PS1(wi2));
        wr3=wa3[i], wi3=fsign*wa3[i+1]; 
        ch[i + 2*l1ido] = cr3;
        ch[i + 2*l1ido + 1] = ci3;

        VCPLXMUL(cr4, ci4, LD_PS1(wr3), LD_PS1(wi3));
        ch[i + 3*l1ido] = cr4;
        ch[i + 3*l1ido + 1] = ci4;
      }
    }
  }
} /* passf4 */

/*
  passf5 and passb5 has been merged here, fsign = -1 for passf5, +1 for passb5
*/
static NEVER_INLINE(void) passf5_ps(int ido, int l1, const v4sd *cc, v4sd *ch,
                                    const double *wa1, const double *wa2, 
                                    const double *wa3, const double *wa4, double fsign) {  
  static const double tr11 = .309016994374947f;
  const double ti11 = .951056516295154f*fsign;
  static const double tr12 = -.809016994374947f;
  const double ti12 = .587785252292473f*fsign;

  /* Local variables */
  int i, k;
  v4sd ci2, ci3, ci4, ci5, di3, di4, di5, di2, cr2, cr3, cr5, cr4, ti2, ti3,
    ti4, ti5, dr3, dr4, dr5, dr2, tr2, tr3, tr4, tr5;

  double wr1, wi1, wr2, wi2, wr3, wi3, wr4, wi4;

#define cc_ref(a_1,a_2) cc[(a_2-1)*ido + a_1 + 1]
#define ch_ref(a_1,a_3) ch[(a_3-1)*l1*ido + a_1 + 1]

  assert(ido > 2);
  for (k = 0; k < l1; ++k, cc += 5*ido, ch += ido) {
    for (i = 0; i < ido-1; i += 2) {
      ti5 = VSUB(cc_ref(i  , 2), cc_ref(i  , 5));
      ti2 = VADD(cc_ref(i  , 2), cc_ref(i  , 5));
      ti4 = VSUB(cc_ref(i  , 3), cc_ref(i  , 4));
      ti3 = VADD(cc_ref(i  , 3), cc_ref(i  , 4));
      tr5 = VSUB(cc_ref(i-1, 2), cc_ref(i-1, 5));
      tr2 = VADD(cc_ref(i-1, 2), cc_ref(i-1, 5));
      tr4 = VSUB(cc_ref(i-1, 3), cc_ref(i-1, 4));
      tr3 = VADD(cc_ref(i-1, 3), cc_ref(i-1, 4));
      ch_ref(i-1, 1) = VADD(cc_ref(i-1, 1), VADD(tr2, tr3));
      ch_ref(i  , 1) = VADD(cc_ref(i  , 1), VADD(ti2, ti3));
      cr2 = VADD(cc_ref(i-1, 1), VADD(SVMUL(tr11, tr2),SVMUL(tr12, tr3)));
      ci2 = VADD(cc_ref(i  , 1), VADD(SVMUL(tr11, ti2),SVMUL(tr12, ti3)));
      cr3 = VADD(cc_ref(i-1, 1), VADD(SVMUL(tr12, tr2),SVMUL(tr11, tr3)));
      ci3 = VADD(cc_ref(i  , 1), VADD(SVMUL(tr12, ti2),SVMUL(tr11, ti3)));
      cr5 = VADD(SVMUL(ti11, tr5), SVMUL(ti12, tr4));
      ci5 = VADD(SVMUL(ti11, ti5), SVMUL(ti12, ti4));
      cr4 = VSUB(SVMUL(ti12, tr5), SVMUL(ti11, tr4));
      ci4 = VSUB(SVMUL(ti12, ti5), SVMUL(ti11, ti4));
      dr3 = VSUB(cr3, ci4);
      dr4 = VADD(cr3, ci4);
      di3 = VADD(ci3, cr4);
      di4 = VSUB(ci3, cr4);
      dr5 = VADD(cr2, ci5);
      dr2 = VSUB(cr2, ci5);
      di5 = VSUB(ci2, cr5);
      di2 = VADD(ci2, cr5);
      wr1=wa1[i], wi1=fsign*wa1[i+1], wr2=wa2[i], wi2=fsign*wa2[i+1]; 
      wr3=wa3[i], wi3=fsign*wa3[i+1], wr4=wa4[i], wi4=fsign*wa4[i+1]; 
      VCPLXMUL(dr2, di2, LD_PS1(wr1), LD_PS1(wi1));
      ch_ref(i - 1, 2) = dr2;
      ch_ref(i, 2)     = di2;
      VCPLXMUL(dr3, di3, LD_PS1(wr2), LD_PS1(wi2));
      ch_ref(i - 1, 3) = dr3;
      ch_ref(i, 3)     = di3;
      VCPLXMUL(dr4, di4, LD_PS1(wr3), LD_PS1(wi3));
      ch_ref(i - 1, 4) = dr4;
      ch_ref(i, 4)     = di4;
      VCPLXMUL(dr5, di5, LD_PS1(wr4), LD_PS1(wi4));
      ch_ref(i - 1, 5) = dr5;
      ch_ref(i, 5)     = di5;
    }
  }
#undef ch_ref
#undef cc_ref
}

static NEVER_INLINE(void) radf2_ps(int ido, int l1, const v4sd * RESTRICT cc, v4sd * RESTRICT ch, const double *wa1) {
  static const double minus_one = -1.f;
  int i, k, l1ido = l1*ido;
  for (k=0; k < l1ido; k += ido) {
    v4sd a = cc[k], b = cc[k + l1ido];
    ch[2*k] = VADD(a, b);
    ch[2*(k+ido)-1] = VSUB(a, b);
  }
  if (ido < 2) return;
  if (ido != 2) {
    for (k=0; k < l1ido; k += ido) {
      for (i=2; i<ido; i+=2) {
        v4sd tr2 = cc[i - 1 + k + l1ido], ti2 = cc[i + k + l1ido];
        v4sd br = cc[i - 1 + k], bi = cc[i + k];
        VCPLXMULCONJ(tr2, ti2, LD_PS1(wa1[i - 2]), LD_PS1(wa1[i - 1])); 
        ch[i + 2*k] = VADD(bi, ti2);
        ch[2*(k+ido) - i] = VSUB(ti2, bi);
        ch[i - 1 + 2*k] = VADD(br, tr2);
        ch[2*(k+ido) - i -1] = VSUB(br, tr2);
      }
    }
    if (ido % 2 == 1) return;
  }
  for (k=0; k < l1ido; k += ido) {
    ch[2*k + ido] = SVMUL(minus_one, cc[ido-1 + k + l1ido]);
    ch[2*k + ido-1] = cc[k + ido-1];
  }
} /* radf2 */


static NEVER_INLINE(void) radb2_ps(int ido, int l1, const v4sd *cc, v4sd *ch, const double *wa1) {
  static const double minus_two=-2;
  int i, k, l1ido = l1*ido;
  v4sd a,b,c,d, tr2, ti2;
  for (k=0; k < l1ido; k += ido) {
    a = cc[2*k]; b = cc[2*(k+ido) - 1];
    ch[k] = VADD(a, b);
    ch[k + l1ido] =VSUB(a, b);
  }
  if (ido < 2) return;
  if (ido != 2) {
    for (k = 0; k < l1ido; k += ido) {
      for (i = 2; i < ido; i += 2) {
        a = cc[i-1 + 2*k]; b = cc[2*(k + ido) - i - 1];
        c = cc[i+0 + 2*k]; d = cc[2*(k + ido) - i + 0];
        ch[i-1 + k] = VADD(a, b);
        tr2 = VSUB(a, b);
        ch[i+0 + k] = VSUB(c, d);
        ti2 = VADD(c, d);
        VCPLXMUL(tr2, ti2, LD_PS1(wa1[i - 2]), LD_PS1(wa1[i - 1]));
        ch[i-1 + k + l1ido] = tr2;
        ch[i+0 + k + l1ido] = ti2;
      }
    }
    if (ido % 2 == 1) return;
  }
  for (k = 0; k < l1ido; k += ido) {
    a = cc[2*k + ido-1]; b = cc[2*k + ido];
    ch[k + ido-1] = VADD(a,a);
    ch[k + ido-1 + l1ido] = SVMUL(minus_two, b);
  }
} /* radb2 */

static void radf3_ps(int ido, int l1, const v4sd * RESTRICT cc, v4sd * RESTRICT ch,
                     const double *wa1, const double *wa2) {
  static const double taur = -0.5f;
  static const double taui = 0.866025403784439f;
  int i, k, ic;
  v4sd ci2, di2, di3, cr2, dr2, dr3, ti2, ti3, tr2, tr3, wr1, wi1, wr2, wi2;
  for (k=0; k<l1; k++) {
    cr2 = VADD(cc[(k + l1)*ido], cc[(k + 2*l1)*ido]);
    ch[3*k*ido] = VADD(cc[k*ido], cr2);
    ch[(3*k+2)*ido] = SVMUL(taui, VSUB(cc[(k + l1*2)*ido], cc[(k + l1)*ido]));
    ch[ido-1 + (3*k + 1)*ido] = VADD(cc[k*ido], SVMUL(taur, cr2));
  }
  if (ido == 1) return;
  for (k=0; k<l1; k++) {
    for (i=2; i<ido; i+=2) {
      ic = ido - i;
      wr1 = LD_PS1(wa1[i - 2]); wi1 = LD_PS1(wa1[i - 1]);
      dr2 = cc[i - 1 + (k + l1)*ido]; di2 = cc[i + (k + l1)*ido];
      VCPLXMULCONJ(dr2, di2, wr1, wi1);

      wr2 = LD_PS1(wa2[i - 2]); wi2 = LD_PS1(wa2[i - 1]);
      dr3 = cc[i - 1 + (k + l1*2)*ido]; di3 = cc[i + (k + l1*2)*ido];
      VCPLXMULCONJ(dr3, di3, wr2, wi2);
        
      cr2 = VADD(dr2, dr3);
      ci2 = VADD(di2, di3);
      ch[i - 1 + 3*k*ido] = VADD(cc[i - 1 + k*ido], cr2);
      ch[i + 3*k*ido] = VADD(cc[i + k*ido], ci2);
      tr2 = VADD(cc[i - 1 + k*ido], SVMUL(taur, cr2));
      ti2 = VADD(cc[i + k*ido], SVMUL(taur, ci2));
      tr3 = SVMUL(taui, VSUB(di2, di3));
      ti3 = SVMUL(taui, VSUB(dr3, dr2));
      ch[i - 1 + (3*k + 2)*ido] = VADD(tr2, tr3);
      ch[ic - 1 + (3*k + 1)*ido] = VSUB(tr2, tr3);
      ch[i + (3*k + 2)*ido] = VADD(ti2, ti3);
      ch[ic + (3*k + 1)*ido] = VSUB(ti3, ti2);
    }
  }
} /* radf3 */


static void radb3_ps(int ido, int l1, const v4sd *RESTRICT cc, v4sd *RESTRICT ch,
                     const double *wa1, const double *wa2)
{
  static const double taur = -0.5f;
  static const double taui = 0.866025403784439f;
  static const double taui_2 = 0.866025403784439f*2;
  int i, k, ic;
  v4sd ci2, ci3, di2, di3, cr2, cr3, dr2, dr3, ti2, tr2;
  for (k=0; k<l1; k++) {
    tr2 = cc[ido-1 + (3*k + 1)*ido]; tr2 = VADD(tr2,tr2);
    cr2 = VMADD(LD_PS1(taur), tr2, cc[3*k*ido]);
    ch[k*ido] = VADD(cc[3*k*ido], tr2);
    ci3 = SVMUL(taui_2, cc[(3*k + 2)*ido]);
    ch[(k + l1)*ido] = VSUB(cr2, ci3);
    ch[(k + 2*l1)*ido] = VADD(cr2, ci3);
  }
  if (ido == 1) return;
  for (k=0; k<l1; k++) {
    for (i=2; i<ido; i+=2) {
      ic = ido - i;
      tr2 = VADD(cc[i - 1 + (3*k + 2)*ido], cc[ic - 1 + (3*k + 1)*ido]);
      cr2 = VMADD(LD_PS1(taur), tr2, cc[i - 1 + 3*k*ido]);
      ch[i - 1 + k*ido] = VADD(cc[i - 1 + 3*k*ido], tr2);
      ti2 = VSUB(cc[i + (3*k + 2)*ido], cc[ic + (3*k + 1)*ido]);
      ci2 = VMADD(LD_PS1(taur), ti2, cc[i + 3*k*ido]);
      ch[i + k*ido] = VADD(cc[i + 3*k*ido], ti2);
      cr3 = SVMUL(taui, VSUB(cc[i - 1 + (3*k + 2)*ido], cc[ic - 1 + (3*k + 1)*ido]));
      ci3 = SVMUL(taui, VADD(cc[i + (3*k + 2)*ido], cc[ic + (3*k + 1)*ido]));
      dr2 = VSUB(cr2, ci3);
      dr3 = VADD(cr2, ci3);
      di2 = VADD(ci2, cr3);
      di3 = VSUB(ci2, cr3);
      VCPLXMUL(dr2, di2, LD_PS1(wa1[i-2]), LD_PS1(wa1[i-1]));
      ch[i - 1 + (k + l1)*ido] = dr2;
      ch[i + (k + l1)*ido] = di2;
      VCPLXMUL(dr3, di3, LD_PS1(wa2[i-2]), LD_PS1(wa2[i-1]));
      ch[i - 1 + (k + 2*l1)*ido] = dr3;
      ch[i + (k + 2*l1)*ido] = di3;
    }
  }
} /* radb3 */

static NEVER_INLINE(void) radf4_ps(int ido, int l1, const v4sd *RESTRICT cc, v4sd * RESTRICT ch,
                                   const double * RESTRICT wa1, const double * RESTRICT wa2, const double * RESTRICT wa3)
{
  static const double minus_hsqt2 = (double)-0.7071067811865475;
  int i, k, l1ido = l1*ido;
  {
    const v4sd *RESTRICT cc_ = cc, * RESTRICT cc_end = cc + l1ido; 
    v4sd * RESTRICT ch_ = ch;
    while (cc < cc_end) {
      // this loop represents between 25% and 40% of total radf4_ps cost !
      v4sd a0 = cc[0], a1 = cc[l1ido];
      v4sd a2 = cc[2*l1ido], a3 = cc[3*l1ido];
      v4sd tr1 = VADD(a1, a3);
      v4sd tr2 = VADD(a0, a2);
      ch[2*ido-1] = VSUB(a0, a2);
      ch[2*ido  ] = VSUB(a3, a1);
      ch[0      ] = VADD(tr1, tr2);
      ch[4*ido-1] = VSUB(tr2, tr1);
      cc += ido; ch += 4*ido;
    }
    cc = cc_; ch = ch_;
  }
  if (ido < 2) return;
  if (ido != 2) {
    for (k = 0; k < l1ido; k += ido) {
      const v4sd * RESTRICT pc = (v4sd*)(cc + 1 + k);
      for (i=2; i<ido; i += 2, pc += 2) {
        int ic = ido - i;
        v4sd wr, wi, cr2, ci2, cr3, ci3, cr4, ci4;
        v4sd tr1, ti1, tr2, ti2, tr3, ti3, tr4, ti4;

        cr2 = pc[1*l1ido+0];
        ci2 = pc[1*l1ido+1];
        wr=LD_PS1(wa1[i - 2]);
        wi=LD_PS1(wa1[i - 1]);
        VCPLXMULCONJ(cr2,ci2,wr,wi);

        cr3 = pc[2*l1ido+0];
        ci3 = pc[2*l1ido+1];
        wr = LD_PS1(wa2[i-2]); 
        wi = LD_PS1(wa2[i-1]);
        VCPLXMULCONJ(cr3, ci3, wr, wi);

        cr4 = pc[3*l1ido];
        ci4 = pc[3*l1ido+1];
        wr = LD_PS1(wa3[i-2]); 
        wi = LD_PS1(wa3[i-1]);
        VCPLXMULCONJ(cr4, ci4, wr, wi);

        /* at this point, on SSE, five of "cr2 cr3 cr4 ci2 ci3 ci4" should be loaded in registers */

        tr1 = VADD(cr2,cr4);
        tr4 = VSUB(cr4,cr2); 
        tr2 = VADD(pc[0],cr3);
        tr3 = VSUB(pc[0],cr3);
        ch[i - 1 + 4*k] = VADD(tr1,tr2);
        ch[ic - 1 + 4*k + 3*ido] = VSUB(tr2,tr1); // at this point tr1 and tr2 can be disposed
        ti1 = VADD(ci2,ci4);
        ti4 = VSUB(ci2,ci4);
        ch[i - 1 + 4*k + 2*ido] = VADD(ti4,tr3);
        ch[ic - 1 + 4*k + 1*ido] = VSUB(tr3,ti4); // dispose tr3, ti4
        ti2 = VADD(pc[1],ci3);
        ti3 = VSUB(pc[1],ci3);
        ch[i + 4*k] = VADD(ti1, ti2);
        ch[ic + 4*k + 3*ido] = VSUB(ti1, ti2);
        ch[i + 4*k + 2*ido] = VADD(tr4, ti3);
        ch[ic + 4*k + 1*ido] = VSUB(tr4, ti3);
      }
    }
    if (ido % 2 == 1) return;
  }
  for (k=0; k<l1ido; k += ido) {
    v4sd a = cc[ido-1 + k + l1ido], b = cc[ido-1 + k + 3*l1ido];
    v4sd c = cc[ido-1 + k], d = cc[ido-1 + k + 2*l1ido];
    v4sd ti1 = SVMUL(minus_hsqt2, VADD(a, b));
    v4sd tr1 = SVMUL(minus_hsqt2, VSUB(b, a));
    ch[ido-1 + 4*k] = VADD(tr1, c);
    ch[ido-1 + 4*k + 2*ido] = VSUB(c, tr1);
    ch[4*k + 1*ido] = VSUB(ti1, d); 
    ch[4*k + 3*ido] = VADD(ti1, d); 
  }
} /* radf4 */


static NEVER_INLINE(void) radb4_ps(int ido, int l1, const v4sd * RESTRICT cc, v4sd * RESTRICT ch,
                                   const double * RESTRICT wa1, const double * RESTRICT wa2, const double *RESTRICT wa3)
{
  static const double minus_sqrt2 = (double)-1.414213562373095;
  static const double two = 2.f;
  int i, k, l1ido = l1*ido;
  v4sd ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;
  {
    const v4sd *RESTRICT cc_ = cc, * RESTRICT ch_end = ch + l1ido; 
    v4sd *ch_ = ch;
    while (ch < ch_end) {
      v4sd a = cc[0], b = cc[4*ido-1];
      v4sd c = cc[2*ido], d = cc[2*ido-1];
      tr3 = SVMUL(two,d);
      tr2 = VADD(a,b);
      tr1 = VSUB(a,b);
      tr4 = SVMUL(two,c);
      ch[0*l1ido] = VADD(tr2, tr3);
      ch[2*l1ido] = VSUB(tr2, tr3);
      ch[1*l1ido] = VSUB(tr1, tr4);
      ch[3*l1ido] = VADD(tr1, tr4);
      
      cc += 4*ido; ch += ido;
    }
    cc = cc_; ch = ch_;
  }
  if (ido < 2) return;
  if (ido != 2) {
    for (k = 0; k < l1ido; k += ido) {
      const v4sd * RESTRICT pc = (v4sd*)(cc - 1 + 4*k);
      v4sd * RESTRICT ph = (v4sd*)(ch + k + 1);
      for (i = 2; i < ido; i += 2) {

        tr1 = VSUB(pc[i], pc[4*ido - i]);
        tr2 = VADD(pc[i], pc[4*ido - i]);
        ti4 = VSUB(pc[2*ido + i], pc[2*ido - i]);
        tr3 = VADD(pc[2*ido + i], pc[2*ido - i]);
        ph[0] = VADD(tr2, tr3);
        cr3 = VSUB(tr2, tr3);

        ti3 = VSUB(pc[2*ido + i + 1], pc[2*ido - i + 1]);
        tr4 = VADD(pc[2*ido + i + 1], pc[2*ido - i + 1]);
        cr2 = VSUB(tr1, tr4);
        cr4 = VADD(tr1, tr4);

        ti1 = VADD(pc[i + 1], pc[4*ido - i + 1]);
        ti2 = VSUB(pc[i + 1], pc[4*ido - i + 1]);

        ph[1] = VADD(ti2, ti3); ph += l1ido;
        ci3 = VSUB(ti2, ti3);
        ci2 = VADD(ti1, ti4);
        ci4 = VSUB(ti1, ti4);
        VCPLXMUL(cr2, ci2, LD_PS1(wa1[i-2]), LD_PS1(wa1[i-1]));
        ph[0] = cr2;
        ph[1] = ci2; ph += l1ido;
        VCPLXMUL(cr3, ci3, LD_PS1(wa2[i-2]), LD_PS1(wa2[i-1]));
        ph[0] = cr3;
        ph[1] = ci3; ph += l1ido;
        VCPLXMUL(cr4, ci4, LD_PS1(wa3[i-2]), LD_PS1(wa3[i-1]));
        ph[0] = cr4;
        ph[1] = ci4; ph = ph - 3*l1ido + 2;
      }
    }
    if (ido % 2 == 1) return;
  }
  for (k=0; k < l1ido; k+=ido) {
    int i0 = 4*k + ido;
    v4sd c = cc[i0-1], d = cc[i0 + 2*ido-1];
    v4sd a = cc[i0+0], b = cc[i0 + 2*ido+0];
    tr1 = VSUB(c,d);
    tr2 = VADD(c,d);
    ti1 = VADD(b,a);
    ti2 = VSUB(b,a);
    ch[ido-1 + k + 0*l1ido] = VADD(tr2,tr2);
    ch[ido-1 + k + 1*l1ido] = SVMUL(minus_sqrt2, VSUB(ti1, tr1));
    ch[ido-1 + k + 2*l1ido] = VADD(ti2, ti2);
    ch[ido-1 + k + 3*l1ido] = SVMUL(minus_sqrt2, VADD(ti1, tr1));
  }
} /* radb4 */

static void radf5_ps(int ido, int l1, const v4sd * RESTRICT cc, v4sd * RESTRICT ch, 
                     const double *wa1, const double *wa2, const double *wa3, const double *wa4)
{
  static const double tr11 = .309016994374947f;
  static const double ti11 = .951056516295154f;
  static const double tr12 = -.809016994374947f;
  static const double ti12 = .587785252292473f;

  /* System generated locals */
  int cc_offset, ch_offset;

  /* Local variables */
  int i, k, ic;
  v4sd ci2, di2, ci4, ci5, di3, di4, di5, ci3, cr2, cr3, dr2, dr3, dr4, dr5,
    cr5, cr4, ti2, ti3, ti5, ti4, tr2, tr3, tr4, tr5;
  int idp2;


#define cc_ref(a_1,a_2,a_3) cc[((a_3)*l1 + (a_2))*ido + a_1]
#define ch_ref(a_1,a_2,a_3) ch[((a_3)*5 + (a_2))*ido + a_1]

  /* Parameter adjustments */
  ch_offset = 1 + ido * 6;
  ch -= ch_offset;
  cc_offset = 1 + ido * (1 + l1);
  cc -= cc_offset;

  /* Function Body */
  for (k = 1; k <= l1; ++k) {
    cr2 = VADD(cc_ref(1, k, 5), cc_ref(1, k, 2));
    ci5 = VSUB(cc_ref(1, k, 5), cc_ref(1, k, 2));
    cr3 = VADD(cc_ref(1, k, 4), cc_ref(1, k, 3));
    ci4 = VSUB(cc_ref(1, k, 4), cc_ref(1, k, 3));
    ch_ref(1, 1, k) = VADD(cc_ref(1, k, 1), VADD(cr2, cr3));
    ch_ref(ido, 2, k) = VADD(cc_ref(1, k, 1), VADD(SVMUL(tr11, cr2), SVMUL(tr12, cr3)));
    ch_ref(1, 3, k) = VADD(SVMUL(ti11, ci5), SVMUL(ti12, ci4));
    ch_ref(ido, 4, k) = VADD(cc_ref(1, k, 1), VADD(SVMUL(tr12, cr2), SVMUL(tr11, cr3)));
    ch_ref(1, 5, k) = VSUB(SVMUL(ti12, ci5), SVMUL(ti11, ci4));
    //printf("pffft: radf5, k=%d ch_ref=%f, ci4=%f\n", k, ch_ref(1, 5, k), ci4);
  }
  if (ido == 1) {
    return;
  }
  idp2 = ido + 2;
  for (k = 1; k <= l1; ++k) {
    for (i = 3; i <= ido; i += 2) {
      ic = idp2 - i;
      dr2 = LD_PS1(wa1[i-3]); di2 = LD_PS1(wa1[i-2]);
      dr3 = LD_PS1(wa2[i-3]); di3 = LD_PS1(wa2[i-2]);
      dr4 = LD_PS1(wa3[i-3]); di4 = LD_PS1(wa3[i-2]);
      dr5 = LD_PS1(wa4[i-3]); di5 = LD_PS1(wa4[i-2]);
      VCPLXMULCONJ(dr2, di2, cc_ref(i-1, k, 2), cc_ref(i, k, 2));
      VCPLXMULCONJ(dr3, di3, cc_ref(i-1, k, 3), cc_ref(i, k, 3));
      VCPLXMULCONJ(dr4, di4, cc_ref(i-1, k, 4), cc_ref(i, k, 4));
      VCPLXMULCONJ(dr5, di5, cc_ref(i-1, k, 5), cc_ref(i, k, 5));
      cr2 = VADD(dr2, dr5);
      ci5 = VSUB(dr5, dr2);
      cr5 = VSUB(di2, di5);
      ci2 = VADD(di2, di5);
      cr3 = VADD(dr3, dr4);
      ci4 = VSUB(dr4, dr3);
      cr4 = VSUB(di3, di4);
      ci3 = VADD(di3, di4);
      ch_ref(i - 1, 1, k) = VADD(cc_ref(i - 1, k, 1), VADD(cr2, cr3));
      ch_ref(i, 1, k) = VSUB(cc_ref(i, k, 1), VADD(ci2, ci3));//
      tr2 = VADD(cc_ref(i - 1, k, 1), VADD(SVMUL(tr11, cr2), SVMUL(tr12, cr3)));
      ti2 = VSUB(cc_ref(i, k, 1), VADD(SVMUL(tr11, ci2), SVMUL(tr12, ci3)));//
      tr3 = VADD(cc_ref(i - 1, k, 1), VADD(SVMUL(tr12, cr2), SVMUL(tr11, cr3)));
      ti3 = VSUB(cc_ref(i, k, 1), VADD(SVMUL(tr12, ci2), SVMUL(tr11, ci3)));//
      tr5 = VADD(SVMUL(ti11, cr5), SVMUL(ti12, cr4));
      ti5 = VADD(SVMUL(ti11, ci5), SVMUL(ti12, ci4));
      tr4 = VSUB(SVMUL(ti12, cr5), SVMUL(ti11, cr4));
      ti4 = VSUB(SVMUL(ti12, ci5), SVMUL(ti11, ci4));
      ch_ref(i - 1, 3, k) = VSUB(tr2, tr5);
      ch_ref(ic - 1, 2, k) = VADD(tr2, tr5);
      ch_ref(i, 3, k) = VADD(ti2, ti5);
      ch_ref(ic, 2, k) = VSUB(ti5, ti2);
      ch_ref(i - 1, 5, k) = VSUB(tr3, tr4);
      ch_ref(ic - 1, 4, k) = VADD(tr3, tr4);
      ch_ref(i, 5, k) = VADD(ti3, ti4);
      ch_ref(ic, 4, k) = VSUB(ti4, ti3);
    }
  }
#undef cc_ref
#undef ch_ref
} /* radf5 */

static void radb5_ps(int ido, int l1, const v4sd *RESTRICT cc, v4sd *RESTRICT ch, 
                  const double *wa1, const double *wa2, const double *wa3, const double *wa4)
{
  static const double tr11 = .309016994374947f;
  static const double ti11 = .951056516295154f;
  static const double tr12 = -.809016994374947f;
  static const double ti12 = .587785252292473f;

  int cc_offset, ch_offset;

  /* Local variables */
  int i, k, ic;
  v4sd ci2, ci3, ci4, ci5, di3, di4, di5, di2, cr2, cr3, cr5, cr4, ti2, ti3,
    ti4, ti5, dr3, dr4, dr5, dr2, tr2, tr3, tr4, tr5;
  int idp2;

#define cc_ref(a_1,a_2,a_3) cc[((a_3)*5 + (a_2))*ido + a_1]
#define ch_ref(a_1,a_2,a_3) ch[((a_3)*l1 + (a_2))*ido + a_1]

  /* Parameter adjustments */
  ch_offset = 1 + ido * (1 + l1);
  ch -= ch_offset;
  cc_offset = 1 + ido * 6;
  cc -= cc_offset;

  /* Function Body */
  for (k = 1; k <= l1; ++k) {
    ti5 = VADD(cc_ref(1, 3, k), cc_ref(1, 3, k));
    ti4 = VADD(cc_ref(1, 5, k), cc_ref(1, 5, k));
    tr2 = VADD(cc_ref(ido, 2, k), cc_ref(ido, 2, k));
    tr3 = VADD(cc_ref(ido, 4, k), cc_ref(ido, 4, k));
    ch_ref(1, k, 1) = VADD(cc_ref(1, 1, k), VADD(tr2, tr3));
    cr2 = VADD(cc_ref(1, 1, k), VADD(SVMUL(tr11, tr2), SVMUL(tr12, tr3)));
    cr3 = VADD(cc_ref(1, 1, k), VADD(SVMUL(tr12, tr2), SVMUL(tr11, tr3)));
    ci5 = VADD(SVMUL(ti11, ti5), SVMUL(ti12, ti4));
    ci4 = VSUB(SVMUL(ti12, ti5), SVMUL(ti11, ti4));
    ch_ref(1, k, 2) = VSUB(cr2, ci5);
    ch_ref(1, k, 3) = VSUB(cr3, ci4);
    ch_ref(1, k, 4) = VADD(cr3, ci4);
    ch_ref(1, k, 5) = VADD(cr2, ci5);
  }
  if (ido == 1) {
    return;
  }
  idp2 = ido + 2;
  for (k = 1; k <= l1; ++k) {
    for (i = 3; i <= ido; i += 2) {
      ic = idp2 - i;
      ti5 = VADD(cc_ref(i  , 3, k), cc_ref(ic  , 2, k));
      ti2 = VSUB(cc_ref(i  , 3, k), cc_ref(ic  , 2, k));
      ti4 = VADD(cc_ref(i  , 5, k), cc_ref(ic  , 4, k));
      ti3 = VSUB(cc_ref(i  , 5, k), cc_ref(ic  , 4, k));
      tr5 = VSUB(cc_ref(i-1, 3, k), cc_ref(ic-1, 2, k));
      tr2 = VADD(cc_ref(i-1, 3, k), cc_ref(ic-1, 2, k));
      tr4 = VSUB(cc_ref(i-1, 5, k), cc_ref(ic-1, 4, k));
      tr3 = VADD(cc_ref(i-1, 5, k), cc_ref(ic-1, 4, k));
      ch_ref(i - 1, k, 1) = VADD(cc_ref(i-1, 1, k), VADD(tr2, tr3));
      ch_ref(i, k, 1) = VADD(cc_ref(i, 1, k), VADD(ti2, ti3));
      cr2 = VADD(cc_ref(i-1, 1, k), VADD(SVMUL(tr11, tr2), SVMUL(tr12, tr3)));
      ci2 = VADD(cc_ref(i  , 1, k), VADD(SVMUL(tr11, ti2), SVMUL(tr12, ti3)));
      cr3 = VADD(cc_ref(i-1, 1, k), VADD(SVMUL(tr12, tr2), SVMUL(tr11, tr3)));
      ci3 = VADD(cc_ref(i  , 1, k), VADD(SVMUL(tr12, ti2), SVMUL(tr11, ti3)));
      cr5 = VADD(SVMUL(ti11, tr5), SVMUL(ti12, tr4));
      ci5 = VADD(SVMUL(ti11, ti5), SVMUL(ti12, ti4));
      cr4 = VSUB(SVMUL(ti12, tr5), SVMUL(ti11, tr4));
      ci4 = VSUB(SVMUL(ti12, ti5), SVMUL(ti11, ti4));
      dr3 = VSUB(cr3, ci4);
      dr4 = VADD(cr3, ci4);
      di3 = VADD(ci3, cr4);
      di4 = VSUB(ci3, cr4);
      dr5 = VADD(cr2, ci5);
      dr2 = VSUB(cr2, ci5);
      di5 = VSUB(ci2, cr5);
      di2 = VADD(ci2, cr5);
      VCPLXMUL(dr2, di2, LD_PS1(wa1[i-3]), LD_PS1(wa1[i-2]));
      VCPLXMUL(dr3, di3, LD_PS1(wa2[i-3]), LD_PS1(wa2[i-2]));
      VCPLXMUL(dr4, di4, LD_PS1(wa3[i-3]), LD_PS1(wa3[i-2]));
      VCPLXMUL(dr5, di5, LD_PS1(wa4[i-3]), LD_PS1(wa4[i-2]));

      ch_ref(i-1, k, 2) = dr2; ch_ref(i, k, 2) = di2;
      ch_ref(i-1, k, 3) = dr3; ch_ref(i, k, 3) = di3;
      ch_ref(i-1, k, 4) = dr4; ch_ref(i, k, 4) = di4;
      ch_ref(i-1, k, 5) = dr5; ch_ref(i, k, 5) = di5;
    }
  }
#undef cc_ref
#undef ch_ref
} /* radb5 */

static NEVER_INLINE(v4sd *) rfftf1_ps(int n, const v4sd *input_readonly, v4sd *work1, v4sd *work2, 
                                      const double *wa, const int *ifac) {  
  v4sd *in  = (v4sd*)input_readonly;
  v4sd *out = (in == work2 ? work1 : work2);
  int nf = ifac[1], k1;
  int l2 = n;
  int iw = n-1;
  assert(in != out && work1 != work2);
  for (k1 = 1; k1 <= nf; ++k1) {
    int kh = nf - k1;
    int ip = ifac[kh + 2];
    int l1 = l2 / ip;
    int ido = n / l2;
    iw -= (ip - 1)*ido;
    switch (ip) {
      case 5: {
        int ix2 = iw + ido;
        int ix3 = ix2 + ido;
        int ix4 = ix3 + ido;
        radf5_ps(ido, l1, in, out, &wa[iw], &wa[ix2], &wa[ix3], &wa[ix4]);
      } break;
      case 4: {
        int ix2 = iw + ido;
        int ix3 = ix2 + ido;
        radf4_ps(ido, l1, in, out, &wa[iw], &wa[ix2], &wa[ix3]);
      } break;
      case 3: {
        int ix2 = iw + ido;
        radf3_ps(ido, l1, in, out, &wa[iw], &wa[ix2]);
      } break;
      case 2:
        radf2_ps(ido, l1, in, out, &wa[iw]);
        break;
      default:
        assert(0);
        break;
    }
    l2 = l1;
    if (out == work2) {
      out = work1; in = work2;
    } else {
      out = work2; in = work1;
    }
  }
  return in; /* this is in fact the output .. */
} /* rfftf1 */

static NEVER_INLINE(v4sd *) rfftb1_ps(int n, const v4sd *input_readonly, v4sd *work1, v4sd *work2, 
                                      const double *wa, const int *ifac) {  
  v4sd *in  = (v4sd*)input_readonly;
  v4sd *out = (in == work2 ? work1 : work2);
  int nf = ifac[1], k1;
  int l1 = 1;
  int iw = 0;
  assert(in != out);
  for (k1=1; k1<=nf; k1++) {
    int ip = ifac[k1 + 1];
    int l2 = ip*l1;
    int ido = n / l2;
    switch (ip) {
      case 5: {
        int ix2 = iw + ido;
        int ix3 = ix2 + ido;
        int ix4 = ix3 + ido;
        radb5_ps(ido, l1, in, out, &wa[iw], &wa[ix2], &wa[ix3], &wa[ix4]);
      } break;
      case 4: {
        int ix2 = iw + ido;
        int ix3 = ix2 + ido;
        radb4_ps(ido, l1, in, out, &wa[iw], &wa[ix2], &wa[ix3]);
      } break;
      case 3: {
        int ix2 = iw + ido;
        radb3_ps(ido, l1, in, out, &wa[iw], &wa[ix2]);
      } break;
      case 2:
        radb2_ps(ido, l1, in, out, &wa[iw]);
        break;
      default:
        assert(0);
        break;
    }
    l1 = l2;
    iw += (ip - 1)*ido;

    if (out == work2) {
      out = work1; in = work2;
    } else {
      out = work2; in = work1;
    }
  }
  return in; /* this is in fact the output .. */
}

static int decompose(int n, int *ifac, const int *ntryh) {
  int nl = n, nf = 0, i, j = 0;
  for (j=0; ntryh[j]; ++j) {
    int ntry = ntryh[j];
    while (nl != 1) {
      int nq = nl / ntry;
      int nr = nl - ntry * nq;
      if (nr == 0) {
        ifac[2+nf++] = ntry;
        nl = nq;
        if (ntry == 2 && nf != 1) {
          for (i = 2; i <= nf; ++i) {
            int ib = nf - i + 2;
            ifac[ib + 1] = ifac[ib];
          }
          ifac[2] = 2;
        }
      } else break;
    }
  }
  ifac[0] = n;
  ifac[1] = nf;  
  return nf;
}



static void rffti1_ps(int n, double *wa, int *ifac)
{
  static const int ntryh[] = { 4,2,3,5,0 };
  int k1, j, ii;

  int nf = decompose(n,ifac,ntryh);
  double argh = (2*M_PI) / n;
  int is = 0;
  int nfm1 = nf - 1;
  int l1 = 1;
  for (k1 = 1; k1 <= nfm1; k1++) {
    int ip = ifac[k1 + 1];
    int ld = 0;
    int l2 = l1*ip;
    int ido = n / l2;
    int ipm = ip - 1;
    for (j = 1; j <= ipm; ++j) {
      double argld;
      int i = is, fi=0;
      ld += l1;
      argld = ld*argh;
      for (ii = 3; ii <= ido; ii += 2) {
        i += 2;
        fi += 1;
        wa[i - 2] = cos(fi*argld);
        wa[i - 1] = sin(fi*argld);
      }
      is += ido;
    }
    l1 = l2;
  }
} /* rffti1 */

static void cffti1_ps(int n, double *wa, int *ifac)
{
  static const int ntryh[] = { 5,3,4,2,0 };
  int k1, j, ii;

  int nf = decompose(n,ifac,ntryh);
  double argh = (2*M_PI)/(double)n;
  int i = 1;
  int l1 = 1;
  for (k1=1; k1<=nf; k1++) {
    int ip = ifac[k1+1];
    int ld = 0;
    int l2 = l1*ip;
    int ido = n / l2;
    int idot = ido + ido + 2;
    int ipm = ip - 1;
    for (j=1; j<=ipm; j++) {
      double argld;
      int i1 = i, fi = 0;
      wa[i-1] = 1;
      wa[i] = 0;
      ld += l1;
      argld = ld*argh;
      for (ii = 4; ii <= idot; ii += 2) {
        i += 2;
        fi += 1;
        wa[i-1] = cos(fi*argld);
        wa[i] = sin(fi*argld);
      }
      if (ip > 5) {
        wa[i1-1] = wa[i-1];
        wa[i1] = wa[i];
      }
    }
    l1 = l2;
  }
} /* cffti1 */


static v4sd *cfftf1_ps(int n, const v4sd *input_readonly, v4sd *work1, v4sd *work2, const double *wa, const int *ifac, int isign) {
  v4sd *in  = (v4sd*)input_readonly;
  v4sd *out = (in == work2 ? work1 : work2); 
  int nf = ifac[1], k1;
  int l1 = 1;
  int iw = 0;
  assert(in != out && work1 != work2);
  for (k1=2; k1<=nf+1; k1++) {
    int ip = ifac[k1];
    int l2 = ip*l1;
    int ido = n / l2;
    int idot = ido + ido;
    switch (ip) {
      case 5: {
        int ix2 = iw + idot;
        int ix3 = ix2 + idot;
        int ix4 = ix3 + idot;
        passf5_ps(idot, l1, in, out, &wa[iw], &wa[ix2], &wa[ix3], &wa[ix4], isign);
      } break;
      case 4: {
        int ix2 = iw + idot;
        int ix3 = ix2 + idot;
        passf4_ps(idot, l1, in, out, &wa[iw], &wa[ix2], &wa[ix3], isign);
      } break;
      case 2: {
        passf2_ps(idot, l1, in, out, &wa[iw], isign);
      } break;
      case 3: {
        int ix2 = iw + idot;
        passf3_ps(idot, l1, in, out, &wa[iw], &wa[ix2], isign);
      } break;
      default:
        assert(0);
    }
    l1 = l2;
    iw += (ip - 1)*idot;
    if (out == work2) {
      out = work1; in = work2;
    } else {
      out = work2; in = work1;
    }
  }

  return in; /* this is in fact the output .. */
}


struct PFFFTD_Setup {
  int     N;
  int     Ncvec; // nb of complex simd vectors (N/4 if PFFFTD_COMPLEX, N/8 if PFFFTD_REAL)
  int ifac[15];
  pffftd_transform_t transform;
  v4sd *data; // allocated room for twiddle coefs
  double *e;    // points into 'data' , N/4*3 elements
  double *twiddle; // points into 'data', N/4 elements
};

PFFFTD_Setup *pffftd_new_setup(int N, pffftd_transform_t transform) {
  PFFFTD_Setup *s = (PFFFTD_Setup*)malloc(sizeof(PFFFTD_Setup));
  int k, m;
  /* unfortunately, the fft size must be a multiple of 16 for complex FFTs 
     and 32 for real FFTs -- a lot of stuff would need to be rewritten to
     handle other cases (or maybe just switch to a scalar fft, I don't know..) */
  if (transform == PFFFTD_REAL) { assert((N%(2*SIMD_SZ*SIMD_SZ))==0 && N>0); }
  if (transform == PFFFTD_COMPLEX) { assert((N%(SIMD_SZ*SIMD_SZ))==0 && N>0); }
  //assert((N % 32) == 0);
  s->N = N;
  s->transform = transform;  
  /* nb of complex simd vectors */
  s->Ncvec = (transform == PFFFTD_REAL ? N/2 : N)/SIMD_SZ;
  s->data = (v4sd*)pffftd_aligned_malloc(2*s->Ncvec * sizeof(v4sd));
  s->e = (double*)s->data;
  s->twiddle = (double*)(s->data + (2*s->Ncvec*(SIMD_SZ-1))/SIMD_SZ);  

  if (transform == PFFFTD_REAL) {
    for (k=0; k < s->Ncvec; ++k) {
      int i = k/SIMD_SZ;
      int j = k%SIMD_SZ;
      for (m=0; m < SIMD_SZ-1; ++m) {
        double A = -2*M_PI*(m+1)*k / N;
        s->e[(2*(i*3 + m) + 0) * SIMD_SZ + j] = cos(A);
        s->e[(2*(i*3 + m) + 1) * SIMD_SZ + j] = sin(A);
      }
    }
    rffti1_ps(N/SIMD_SZ, s->twiddle, s->ifac);
  } else {
    for (k=0; k < s->Ncvec; ++k) {
      int i = k/SIMD_SZ;
      int j = k%SIMD_SZ;
      for (m=0; m < SIMD_SZ-1; ++m) {
        double A = -2*M_PI*(m+1)*k / N;
        s->e[(2*(i*3 + m) + 0)*SIMD_SZ + j] = cos(A);
        s->e[(2*(i*3 + m) + 1)*SIMD_SZ + j] = sin(A);
      }
    }
    cffti1_ps(N/SIMD_SZ, s->twiddle, s->ifac);
  }

  /* check that N is decomposable with allowed prime factors */
  for (k=0, m=1; k < s->ifac[1]; ++k) { m *= s->ifac[2+k]; }
  if (m != N/SIMD_SZ) {
    pffftd_destroy_setup(s); s = 0;
  }

  return s;
}


void pffftd_destroy_setup(PFFFTD_Setup *s) {
  pffftd_aligned_free(s->data);
  free(s);
}

#if !defined(PFFFTD_SIMD_DISABLE)

/* [0 0 1 2 3 4 5 6 7 8] -> [0 8 7 6 5 4 3 2 1] */
static void reversed_copy(int N, const v4sd *in, int in_stride, v4sd *out) {
  v4sd g0, g1;
  int k;
  INTERLEAVE2(in[0], in[1], g0, g1); in += in_stride;
  
  *--out = VSWAPHL(g0, g1); // [g0l, g0h], [g1l g1h] -> [g1l, g0h]
  for (k=1; k < N; ++k) {
    v4sd h0, h1;
    INTERLEAVE2(in[0], in[1], h0, h1); in += in_stride;
    *--out = VSWAPHL(g1, h0);
    *--out = VSWAPHL(h0, h1);
    g1 = h1;
  }
  *--out = VSWAPHL(g1, g0);
}

static void unreversed_copy(int N, const v4sd *in, v4sd *out, int out_stride) {
  v4sd g0, g1, h0, h1;
  int k;
  g0 = g1 = in[0]; ++in;
  for (k=1; k < N; ++k) {
    h0 = *in++; h1 = *in++;
    g1 = VSWAPHL(g1, h0);
    h0 = VSWAPHL(h0, h1);
    UNINTERLEAVE2(h0, g1, out[0], out[1]); out += out_stride;
    g1 = h1;
  }
  h0 = *in++; h1 = g0;
  g1 = VSWAPHL(g1, h0);
  h0 = VSWAPHL(h0, h1);
  UNINTERLEAVE2(h0, g1, out[0], out[1]);
}

void pffftd_zreorder(PFFFTD_Setup *setup, const double *in, double *out, pffftd_direction_t direction) {
  int k, N = setup->N, Ncvec = setup->Ncvec;
  const v4sd *vin = (const v4sd*)in;
  v4sd *vout = (v4sd*)out;
  assert(in != out);
  if (setup->transform == PFFFTD_REAL) {
    int k, dk = N/32;
    if (direction == PFFFTD_FORWARD) {
      for (k=0; k < dk; ++k) {
        INTERLEAVE2(vin[k*8 + 0], vin[k*8 + 1], vout[2*(0*dk + k) + 0], vout[2*(0*dk + k) + 1]);
        INTERLEAVE2(vin[k*8 + 4], vin[k*8 + 5], vout[2*(2*dk + k) + 0], vout[2*(2*dk + k) + 1]);
      }
      reversed_copy(dk, vin+2, 8, (v4sd*)(out + N/2));
      reversed_copy(dk, vin+6, 8, (v4sd*)(out + N));
    } else {
      for (k=0; k < dk; ++k) {
        UNINTERLEAVE2(vin[2*(0*dk + k) + 0], vin[2*(0*dk + k) + 1], vout[k*8 + 0], vout[k*8 + 1]);
        UNINTERLEAVE2(vin[2*(2*dk + k) + 0], vin[2*(2*dk + k) + 1], vout[k*8 + 4], vout[k*8 + 5]);
      }
      unreversed_copy(dk, (v4sd*)(in + N/4), (v4sd*)(out + N - 6*SIMD_SZ), -8);
      unreversed_copy(dk, (v4sd*)(in + 3*N/4), (v4sd*)(out + N - 2*SIMD_SZ), -8);
    }
  } else {
    if (direction == PFFFTD_FORWARD) {
      for (k=0; k < Ncvec; ++k) { 
        int kk = (k/4) + (k%4)*(Ncvec/4);
        INTERLEAVE2(vin[k*2], vin[k*2+1], vout[kk*2], vout[kk*2+1]);
      }
    } else {
      for (k=0; k < Ncvec; ++k) { 
        int kk = (k/4) + (k%4)*(Ncvec/4);
        UNINTERLEAVE2(vin[kk*2], vin[kk*2+1], vout[k*2], vout[k*2+1]);
      }
    }
  }
}

void pffftd_cplx_finalize(int Ncvec, const v4sd *in, v4sd *out, const v4sd *e) {
  int k, dk = Ncvec/SIMD_SZ; // number of 4x4 matrix blocks
  v4sd r0, i0, r1, i1, r2, i2, r3, i3;
  v4sd sr0, dr0, sr1, dr1, si0, di0, si1, di1;
  assert(in != out);
  for (k=0; k < dk; ++k) {    
    r0 = in[8*k+0]; i0 = in[8*k+1];
    r1 = in[8*k+2]; i1 = in[8*k+3];
    r2 = in[8*k+4]; i2 = in[8*k+5];
    r3 = in[8*k+6]; i3 = in[8*k+7];
    VTRANSPOSE4(r0,r1,r2,r3);
    VTRANSPOSE4(i0,i1,i2,i3);
    VCPLXMUL(r1,i1,e[k*6+0],e[k*6+1]);
    VCPLXMUL(r2,i2,e[k*6+2],e[k*6+3]);
    VCPLXMUL(r3,i3,e[k*6+4],e[k*6+5]);

    sr0 = VADD(r0,r2); dr0 = VSUB(r0, r2);
    sr1 = VADD(r1,r3); dr1 = VSUB(r1, r3);
    si0 = VADD(i0,i2); di0 = VSUB(i0, i2);
    si1 = VADD(i1,i3); di1 = VSUB(i1, i3);

    /*
      transformation for each column is:
      
      [1   1   1   1   0   0   0   0]   [r0]
      [1   0  -1   0   0  -1   0   1]   [r1]
      [1  -1   1  -1   0   0   0   0]   [r2]
      [1   0  -1   0   0   1   0  -1]   [r3]
      [0   0   0   0   1   1   1   1] * [i0]
      [0   1   0  -1   1   0  -1   0]   [i1]
      [0   0   0   0   1  -1   1  -1]   [i2]
      [0  -1   0   1   1   0  -1   0]   [i3]    
    */
    
    r0 = VADD(sr0, sr1); i0 = VADD(si0, si1);
    r1 = VADD(dr0, di1); i1 = VSUB(di0, dr1);
    r2 = VSUB(sr0, sr1); i2 = VSUB(si0, si1);
    r3 = VSUB(dr0, di1); i3 = VADD(di0, dr1);
  
    *out++ = r0; *out++ = i0; *out++ = r1; *out++ = i1;
    *out++ = r2; *out++ = i2; *out++ = r3; *out++ = i3;
  }
}

void pffftd_cplx_preprocess(int Ncvec, const v4sd *in, v4sd *out, const v4sd *e) {
  int k, dk = Ncvec/SIMD_SZ; // number of 4x4 matrix blocks
  v4sd r0, i0, r1, i1, r2, i2, r3, i3;
  v4sd sr0, dr0, sr1, dr1, si0, di0, si1, di1;
  assert(in != out);
  for (k=0; k < dk; ++k) {    
    r0 = in[8*k+0]; i0 = in[8*k+1];
    r1 = in[8*k+2]; i1 = in[8*k+3];
    r2 = in[8*k+4]; i2 = in[8*k+5];
    r3 = in[8*k+6]; i3 = in[8*k+7];

    sr0 = VADD(r0,r2); dr0 = VSUB(r0, r2);
    sr1 = VADD(r1,r3); dr1 = VSUB(r1, r3);
    si0 = VADD(i0,i2); di0 = VSUB(i0, i2);
    si1 = VADD(i1,i3); di1 = VSUB(i1, i3);

    r0 = VADD(sr0, sr1); i0 = VADD(si0, si1);
    r1 = VSUB(dr0, di1); i1 = VADD(di0, dr1);
    r2 = VSUB(sr0, sr1); i2 = VSUB(si0, si1);
    r3 = VADD(dr0, di1); i3 = VSUB(di0, dr1);

    VCPLXMULCONJ(r1,i1,e[k*6+0],e[k*6+1]);
    VCPLXMULCONJ(r2,i2,e[k*6+2],e[k*6+3]);
    VCPLXMULCONJ(r3,i3,e[k*6+4],e[k*6+5]);

    VTRANSPOSE4(r0,r1,r2,r3);
    VTRANSPOSE4(i0,i1,i2,i3);

    *out++ = r0; *out++ = i0; *out++ = r1; *out++ = i1;
    *out++ = r2; *out++ = i2; *out++ = r3; *out++ = i3;
  }
}


static ALWAYS_INLINE(void) pffftd_real_finalize_4x4(const v4sd *in0, const v4sd *in1, const v4sd *in,
                            const v4sd *e, v4sd *out) {
  v4sd r0, i0, r1, i1, r2, i2, r3, i3;
  v4sd sr0, dr0, sr1, dr1, si0, di0, si1, di1;
  r0 = *in0; i0 = *in1;
  r1 = *in++; i1 = *in++; r2 = *in++; i2 = *in++; r3 = *in++; i3 = *in++;
  VTRANSPOSE4(r0,r1,r2,r3);
  VTRANSPOSE4(i0,i1,i2,i3);
 
  /*
    transformation for each column is:

    [1   1   1   1   0   0   0   0]   [r0]
    [1   0  -1   0   0  -1   0   1]   [r1]
    [1   0  -1   0   0   1   0  -1]   [r2]
    [1  -1   1  -1   0   0   0   0]   [r3]
    [0   0   0   0   1   1   1   1] * [i0]
    [0  -1   0   1  -1   0   1   0]   [i1]
    [0  -1   0   1   1   0  -1   0]   [i2]
    [0   0   0   0  -1   1  -1   1]   [i3]    
  */
  
  //cerr << "matrix initial, before e , REAL:\n 1: " << r0 << "\n 1: " << r1 << "\n 1: " << r2 << "\n 1: " << r3 << "\n";
  //cerr << "matrix initial, before e, IMAG :\n 1: " << i0 << "\n 1: " << i1 << "\n 1: " << i2 << "\n 1: " << i3 << "\n";

  VCPLXMUL(r1,i1,e[0],e[1]);
  VCPLXMUL(r2,i2,e[2],e[3]);
  VCPLXMUL(r3,i3,e[4],e[5]);

  //cerr << "matrix initial, real part:\n 1: " << r0 << "\n 1: " << r1 << "\n 1: " << r2 << "\n 1: " << r3 << "\n";
  //cerr << "matrix initial, imag part:\n 1: " << i0 << "\n 1: " << i1 << "\n 1: " << i2 << "\n 1: " << i3 << "\n";

  sr0 = VADD(r0,r2); dr0 = VSUB(r0,r2); 
  sr1 = VADD(r1,r3); dr1 = VSUB(r3,r1);
  si0 = VADD(i0,i2); di0 = VSUB(i0,i2); 
  si1 = VADD(i1,i3); di1 = VSUB(i3,i1);

  r0 = VADD(sr0, sr1);
  r3 = VSUB(sr0, sr1);
  i0 = VADD(si0, si1);
  i3 = VSUB(si1, si0);
  r1 = VADD(dr0, di1);
  r2 = VSUB(dr0, di1);
  i1 = VSUB(dr1, di0);
  i2 = VADD(dr1, di0);

  *out++ = r0;
  *out++ = i0;
  *out++ = r1;
  *out++ = i1;
  *out++ = r2;
  *out++ = i2;
  *out++ = r3;
  *out++ = i3;

}

static NEVER_INLINE(void) pffftd_real_finalize(int Ncvec, const v4sd *in, v4sd *out, const v4sd *e) {
  int k, dk = Ncvec/SIMD_SZ; // number of 4x4 matrix blocks
  /* fftpack order is f0r f1r f1i f2r f2i ... f(n-1)r f(n-1)i f(n)r */

  v4sd_union cr, ci, *uout = (v4sd_union*)out;
  v4sd save = in[7], zero=VZERO();
  double xr0, xi0, xr1, xi1, xr2, xi2, xr3, xi3;
  static const double s = M_SQRT2/2;

  cr.v = in[0]; ci.v = in[Ncvec*2-1];
  assert(in != out);
  pffftd_real_finalize_4x4(&zero, &zero, in+1, e, out);

  /*
    [cr0 cr1 cr2 cr3 ci0 ci1 ci2 ci3]

    [Xr(1)]  ] [1   1   1   1   0   0   0   0]
    [Xr(N/4) ] [0   0   0   0   1   s   0  -s]
    [Xr(N/2) ] [1   0  -1   0   0   0   0   0]
    [Xr(3N/4)] [0   0   0   0   1  -s   0   s]
    [Xi(1)   ] [1  -1   1  -1   0   0   0   0]
    [Xi(N/4) ] [0   0   0   0   0  -s  -1  -s]
    [Xi(N/2) ] [0  -1   0   1   0   0   0   0]
    [Xi(3N/4)] [0   0   0   0   0  -s   1  -s]
  */

  xr0=(cr.f[0]+cr.f[2]) + (cr.f[1]+cr.f[3]); uout[0].f[0] = xr0;
  xi0=(cr.f[0]+cr.f[2]) - (cr.f[1]+cr.f[3]); uout[1].f[0] = xi0;
  xr2=(cr.f[0]-cr.f[2]);                     uout[4].f[0] = xr2;
  xi2=(cr.f[3]-cr.f[1]);                     uout[5].f[0] = xi2;
  xr1= ci.f[0] + s*(ci.f[1]-ci.f[3]);        uout[2].f[0] = xr1;
  xi1=-ci.f[2] - s*(ci.f[1]+ci.f[3]);        uout[3].f[0] = xi1;
  xr3= ci.f[0] - s*(ci.f[1]-ci.f[3]);        uout[6].f[0] = xr3;
  xi3= ci.f[2] - s*(ci.f[1]+ci.f[3]);        uout[7].f[0] = xi3; 

  for (k=1; k < dk; ++k) {
    v4sd save_next = in[8*k+7];
    pffftd_real_finalize_4x4(&save, &in[8*k+0], in + 8*k+1,
                           e + k*6, out + k*8);
    save = save_next;
  }

}

static ALWAYS_INLINE(void) pffftd_real_preprocess_4x4(const v4sd *in, 
                                             const v4sd *e, v4sd *out, int first) {
  v4sd r0=in[0], i0=in[1], r1=in[2], i1=in[3], r2=in[4], i2=in[5], r3=in[6], i3=in[7];
  /*
    transformation for each column is:

    [1   1   1   1   0   0   0   0]   [r0]
    [1   0   0  -1   0  -1  -1   0]   [r1]
    [1  -1  -1   1   0   0   0   0]   [r2]
    [1   0   0  -1   0   1   1   0]   [r3]
    [0   0   0   0   1  -1   1  -1] * [i0]
    [0  -1   1   0   1   0   0   1]   [i1]
    [0   0   0   0   1   1  -1  -1]   [i2]
    [0   1  -1   0   1   0   0   1]   [i3]    
  */

  v4sd sr0 = VADD(r0,r3), dr0 = VSUB(r0,r3); 
  v4sd sr1 = VADD(r1,r2), dr1 = VSUB(r1,r2);
  v4sd si0 = VADD(i0,i3), di0 = VSUB(i0,i3); 
  v4sd si1 = VADD(i1,i2), di1 = VSUB(i1,i2);

  r0 = VADD(sr0, sr1);
  r2 = VSUB(sr0, sr1);
  r1 = VSUB(dr0, si1);
  r3 = VADD(dr0, si1);
  i0 = VSUB(di0, di1);
  i2 = VADD(di0, di1);
  i1 = VSUB(si0, dr1);
  i3 = VADD(si0, dr1);

  VCPLXMULCONJ(r1,i1,e[0],e[1]);
  VCPLXMULCONJ(r2,i2,e[2],e[3]);
  VCPLXMULCONJ(r3,i3,e[4],e[5]);

  VTRANSPOSE4(r0,r1,r2,r3);
  VTRANSPOSE4(i0,i1,i2,i3);

  if (!first) {
    *out++ = r0;
    *out++ = i0;
  }
  *out++ = r1;
  *out++ = i1;
  *out++ = r2;
  *out++ = i2;
  *out++ = r3;
  *out++ = i3;
}

static NEVER_INLINE(void) pffftd_real_preprocess(int Ncvec, const v4sd *in, v4sd *out, const v4sd *e) {
  int k, dk = Ncvec/SIMD_SZ; // number of 4x4 matrix blocks
  /* fftpack order is f0r f1r f1i f2r f2i ... f(n-1)r f(n-1)i f(n)r */

  v4sd_union Xr, Xi, *uout = (v4sd_union*)out;
  double cr0, ci0, cr1, ci1, cr2, ci2, cr3, ci3;
  static const double s = M_SQRT2;
  assert(in != out);
  for (k=0; k < 4; ++k) {
    Xr.f[k] = ((double*)in)[8*k];
    Xi.f[k] = ((double*)in)[8*k+4];
  }

  pffftd_real_preprocess_4x4(in, e, out+1, 1); // will write only 6 values

  /*
    [Xr0 Xr1 Xr2 Xr3 Xi0 Xi1 Xi2 Xi3]

    [cr0] [1   0   2   0   1   0   0   0]
    [cr1] [1   0   0   0  -1   0  -2   0]
    [cr2] [1   0  -2   0   1   0   0   0]
    [cr3] [1   0   0   0  -1   0   2   0]
    [ci0] [0   2   0   2   0   0   0   0]
    [ci1] [0   s   0  -s   0  -s   0  -s]
    [ci2] [0   0   0   0   0  -2   0   2]
    [ci3] [0  -s   0   s   0  -s   0  -s]
  */
  for (k=1; k < dk; ++k) {    
    pffftd_real_preprocess_4x4(in+8*k, e + k*6, out-1+k*8, 0);
  }

  cr0=(Xr.f[0]+Xi.f[0]) + 2*Xr.f[2]; uout[0].f[0] = cr0;
  cr1=(Xr.f[0]-Xi.f[0]) - 2*Xi.f[2]; uout[0].f[1] = cr1;
  cr2=(Xr.f[0]+Xi.f[0]) - 2*Xr.f[2]; uout[0].f[2] = cr2;
  cr3=(Xr.f[0]-Xi.f[0]) + 2*Xi.f[2]; uout[0].f[3] = cr3;
  ci0= 2*(Xr.f[1]+Xr.f[3]);                       uout[2*Ncvec-1].f[0] = ci0;
  ci1= s*(Xr.f[1]-Xr.f[3]) - s*(Xi.f[1]+Xi.f[3]); uout[2*Ncvec-1].f[1] = ci1;
  ci2= 2*(Xi.f[3]-Xi.f[1]);                       uout[2*Ncvec-1].f[2] = ci2;
  ci3=-s*(Xr.f[1]-Xr.f[3]) - s*(Xi.f[1]+Xi.f[3]); uout[2*Ncvec-1].f[3] = ci3;
}


void pffftd_transform_internal(PFFFTD_Setup *setup, const double *finput, double *foutput, v4sd *scratch,
                             pffftd_direction_t direction, int ordered) {
  int k, Ncvec   = setup->Ncvec;
  int nf_odd = (setup->ifac[1] & 1);

  // temporary buffer is allocated on the stack if the scratch pointer is NULL
  int stack_allocate = (scratch == 0 ? Ncvec*2 : 1);
  VLA_ARRAY_ON_STACK(v4sd, scratch_on_stack, stack_allocate);

  const v4sd *vinput = (const v4sd*)finput;
  v4sd *voutput      = (v4sd*)foutput;
  v4sd *buff[2]      = { voutput, scratch ? scratch : scratch_on_stack };
  int ib = (nf_odd ^ ordered ? 1 : 0);

  assert(VALIGNED(finput) && VALIGNED(foutput));

  //assert(finput != foutput);
  if (direction == PFFFTD_FORWARD) {
    ib = !ib;
    if (setup->transform == PFFFTD_REAL) { 
      ib = (rfftf1_ps(Ncvec*2, vinput, buff[ib], buff[!ib],
                      setup->twiddle, &setup->ifac[0]) == buff[0] ? 0 : 1);      
      pffftd_real_finalize(Ncvec, buff[ib], buff[!ib], (v4sd*)setup->e);
    } else {
      v4sd *tmp = buff[ib];
      for (k=0; k < Ncvec; ++k) {
        UNINTERLEAVE2(vinput[k*2], vinput[k*2+1], tmp[k*2], tmp[k*2+1]);
      }
      ib = (cfftf1_ps(Ncvec, buff[ib], buff[!ib], buff[ib], 
                      setup->twiddle, &setup->ifac[0], -1) == buff[0] ? 0 : 1);
      pffftd_cplx_finalize(Ncvec, buff[ib], buff[!ib], (v4sd*)setup->e);
    }
    if (ordered) {
      pffftd_zreorder(setup, (double*)buff[!ib], (double*)buff[ib], PFFFTD_FORWARD);       
    } else ib = !ib;
  } else {
    if (vinput == buff[ib]) { 
      ib = !ib; // may happen when finput == foutput
    }
    if (ordered) {
      pffftd_zreorder(setup, (double*)vinput, (double*)buff[ib], PFFFTD_BACKWARD); 
      vinput = buff[ib]; ib = !ib;
    }
    if (setup->transform == PFFFTD_REAL) {
      pffftd_real_preprocess(Ncvec, vinput, buff[ib], (v4sd*)setup->e);
      ib = (rfftb1_ps(Ncvec*2, buff[ib], buff[0], buff[1], 
                      setup->twiddle, &setup->ifac[0]) == buff[0] ? 0 : 1);
    } else {
      pffftd_cplx_preprocess(Ncvec, vinput, buff[ib], (v4sd*)setup->e);
      ib = (cfftf1_ps(Ncvec, buff[ib], buff[0], buff[1], 
                      setup->twiddle, &setup->ifac[0], +1) == buff[0] ? 0 : 1);
      for (k=0; k < Ncvec; ++k) {
        INTERLEAVE2(buff[ib][k*2], buff[ib][k*2+1], buff[ib][k*2], buff[ib][k*2+1]);
      }
    }
  }
  
  if (buff[ib] != voutput) {
    /* extra copy required -- this situation should only happen when finput == foutput */
    assert(finput==foutput);
    for (k=0; k < Ncvec; ++k) {
      v4sd a = buff[ib][2*k], b = buff[ib][2*k+1];
      voutput[2*k] = a; voutput[2*k+1] = b;
    }
    ib = !ib;
  }
  assert(buff[ib] == voutput);
}

void pffftd_zconvolve_accumulate(PFFFTD_Setup *s, const double *a, const double *b, double *ab, double scaling) {
  int Ncvec = s->Ncvec;
  const v4sd * RESTRICT va = (const v4sd*)a;
  const v4sd * RESTRICT vb = (const v4sd*)b;
  v4sd * RESTRICT vab = (v4sd*)ab;

#ifdef __arm__
  __builtin_prefetch(va);
  __builtin_prefetch(vb);
  __builtin_prefetch(vab);
  __builtin_prefetch(va+2);
  __builtin_prefetch(vb+2);
  __builtin_prefetch(vab+2);
  __builtin_prefetch(va+4);
  __builtin_prefetch(vb+4);
  __builtin_prefetch(vab+4);
  __builtin_prefetch(va+6);
  __builtin_prefetch(vb+6);
  __builtin_prefetch(vab+6);
# ifndef __clang__
#   define ZCONVOLVE_USING_INLINE_NEON_ASM
# endif
#endif

  double ar, ai, br, bi, abr, abi;
#ifndef ZCONVOLVE_USING_INLINE_ASM
  v4sd vscal = LD_PS1(scaling);
  int i;
#endif

  assert(VALIGNED(a) && VALIGNED(b) && VALIGNED(ab));
  ar = ((v4sd_union*)va)[0].f[0];
  ai = ((v4sd_union*)va)[1].f[0];
  br = ((v4sd_union*)vb)[0].f[0];
  bi = ((v4sd_union*)vb)[1].f[0];
  abr = ((v4sd_union*)vab)[0].f[0];
  abi = ((v4sd_union*)vab)[1].f[0];
 
#ifdef ZCONVOLVE_USING_INLINE_ASM // inline asm version, unfortunately miscompiled by clang 3.2, at least on ubuntu.. so this will be restricted to gcc
  const double *a_ = a, *b_ = b; double *ab_ = ab;
  int N = Ncvec;
  asm volatile("mov         r8, %2                  \n"
               "vdup.f32    q15, %4                 \n"
               "1:                                  \n"
               "pld         [%0,#64]                \n"
               "pld         [%1,#64]                \n"
               "pld         [%2,#64]                \n"
               "pld         [%0,#96]                \n"
               "pld         [%1,#96]                \n"
               "pld         [%2,#96]                \n"
               "vld1.f32    {q0,q1},   [%0,:128]!         \n"
               "vld1.f32    {q4,q5},   [%1,:128]!         \n"
               "vld1.f32    {q2,q3},   [%0,:128]!         \n"
               "vld1.f32    {q6,q7},   [%1,:128]!         \n"
               "vld1.f32    {q8,q9},   [r8,:128]!          \n"
               
               "vmul.f32    q10, q0, q4             \n"
               "vmul.f32    q11, q0, q5             \n"
               "vmul.f32    q12, q2, q6             \n" 
               "vmul.f32    q13, q2, q7             \n"                 
               "vmls.f32    q10, q1, q5             \n"
               "vmla.f32    q11, q1, q4             \n"
               "vld1.f32    {q0,q1}, [r8,:128]!     \n"
               "vmls.f32    q12, q3, q7             \n"
               "vmla.f32    q13, q3, q6             \n"
               "vmla.f32    q8, q10, q15            \n"
               "vmla.f32    q9, q11, q15            \n"
               "vmla.f32    q0, q12, q15            \n"
               "vmla.f32    q1, q13, q15            \n"
               "vst1.f32    {q8,q9},[%2,:128]!    \n"
               "vst1.f32    {q0,q1},[%2,:128]!    \n"
               "subs        %3, #2                  \n"
               "bne         1b                      \n"
               : "+r"(a_), "+r"(b_), "+r"(ab_), "+r"(N) : "r"(scaling) : "r8", "q0","q1","q2","q3","q4","q5","q6","q7","q8","q9", "q10","q11","q12","q13","q15","memory");
#else // default routine, works fine for non-arm cpus with current compilers
  for (i=0; i < Ncvec; i += 2) {
    v4sd ar, ai, br, bi;
    ar = va[2*i+0]; ai = va[2*i+1];
    br = vb[2*i+0]; bi = vb[2*i+1];
    VCPLXMUL(ar, ai, br, bi);
    vab[2*i+0] = VMADD(ar, vscal, vab[2*i+0]);
    vab[2*i+1] = VMADD(ai, vscal, vab[2*i+1]);
    ar = va[2*i+2]; ai = va[2*i+3];
    br = vb[2*i+2]; bi = vb[2*i+3];
    VCPLXMUL(ar, ai, br, bi);
    vab[2*i+2] = VMADD(ar, vscal, vab[2*i+2]);
    vab[2*i+3] = VMADD(ai, vscal, vab[2*i+3]);
  }
#endif
  if (s->transform == PFFFTD_REAL) {
    ((v4sd_union*)vab)[0].f[0] = abr + ar*br*scaling;
    ((v4sd_union*)vab)[1].f[0] = abi + ai*bi*scaling;
  }
}


#else // defined(PFFFTD_SIMD_DISABLE)

// standard routine using scalar floats, without SIMD stuff.

#define pffftd_zreorder_nosimd pffftd_zreorder
void pffftd_zreorder_nosimd(PFFFTD_Setup *setup, const double *in, double *out, pffftd_direction_t direction) {
  int k, N = setup->N;
  if (setup->transform == PFFFTD_COMPLEX) {
    for (k=0; k < 2*N; ++k) out[k] = in[k];
    return;
  }
  else if (direction == PFFFTD_FORWARD) {
    double x_N = in[N-1];
    for (k=N-1; k > 1; --k) out[k] = in[k-1]; 
    out[0] = in[0];
    out[1] = x_N;
  } else {
    double x_N = in[1];
    for (k=1; k < N-1; ++k) out[k] = in[k+1]; 
    out[0] = in[0];
    out[N-1] = x_N;
  }
}

#define pffftd_transform_internal_nosimd pffftd_transform_internal
void pffftd_transform_internal_nosimd(PFFFTD_Setup *setup, const double *input, double *output, double *scratch,
                                    pffftd_direction_t direction, int ordered) {
  int Ncvec   = setup->Ncvec;
  int nf_odd = (setup->ifac[1] & 1);

  // temporary buffer is allocated on the stack if the scratch pointer is NULL
  int stack_allocate = (scratch == 0 ? Ncvec*2 : 1);
  VLA_ARRAY_ON_STACK(v4sd, scratch_on_stack, stack_allocate);
  double *buff[2];
  int ib;
  if (scratch == 0) scratch = scratch_on_stack;
  buff[0] = output; buff[1] = scratch;

  if (setup->transform == PFFFTD_COMPLEX) ordered = 0; // it is always ordered.
  ib = (nf_odd ^ ordered ? 1 : 0);

  if (direction == PFFFTD_FORWARD) {
    if (setup->transform == PFFFTD_REAL) { 
      ib = (rfftf1_ps(Ncvec*2, input, buff[ib], buff[!ib],
                      setup->twiddle, &setup->ifac[0]) == buff[0] ? 0 : 1);      
    } else {
      ib = (cfftf1_ps(Ncvec, input, buff[ib], buff[!ib], 
                      setup->twiddle, &setup->ifac[0], -1) == buff[0] ? 0 : 1);
    }
    if (ordered) {
      pffftd_zreorder(setup, buff[ib], buff[!ib], PFFFTD_FORWARD); ib = !ib;
    }
  } else {    
    if (input == buff[ib]) { 
      ib = !ib; // may happen when finput == foutput
    }
    if (ordered) {
      pffftd_zreorder(setup, input, buff[!ib], PFFFTD_BACKWARD); 
      input = buff[!ib];
    }
    if (setup->transform == PFFFTD_REAL) {
      ib = (rfftb1_ps(Ncvec*2, input, buff[ib], buff[!ib], 
                      setup->twiddle, &setup->ifac[0]) == buff[0] ? 0 : 1);
    } else {
      ib = (cfftf1_ps(Ncvec, input, buff[ib], buff[!ib], 
                      setup->twiddle, &setup->ifac[0], +1) == buff[0] ? 0 : 1);
    }
  }
  if (buff[ib] != output) {
    int k;
    // extra copy required -- this situation should happens only when finput == foutput
    assert(input==output);
    for (k=0; k < Ncvec; ++k) {
      double a = buff[ib][2*k], b = buff[ib][2*k+1];
      output[2*k] = a; output[2*k+1] = b;
    }
    ib = !ib;
  }
  assert(buff[ib] == output);
}

#define pffftd_zconvolve_accumulate_nosimd pffftd_zconvolve_accumulate
void pffftd_zconvolve_accumulate_nosimd(PFFFTD_Setup *s, const double *a, const double *b,
                                       double *ab, double scaling) {
  int i, Ncvec = s->Ncvec;

  if (s->transform == PFFFTD_REAL) {
    // take care of the fftpack ordering
    ab[0] += a[0]*b[0]*scaling;
    ab[2*Ncvec-1] += a[2*Ncvec-1]*b[2*Ncvec-1]*scaling;
    ++ab; ++a; ++b; --Ncvec;
  }
  for (i=0; i < Ncvec; ++i) {
    double ar, ai, br, bi;
    ar = a[2*i+0]; ai = a[2*i+1];
    br = b[2*i+0]; bi = b[2*i+1];
    VCPLXMUL(ar, ai, br, bi);
    ab[2*i+0] += ar*scaling;
    ab[2*i+1] += ai*scaling;
  }
}

#endif // defined(PFFFTD_SIMD_DISABLE)

void pffftd_transform(PFFFTD_Setup *setup, const double *input, double *output, double *work, pffftd_direction_t direction) {
  pffftd_transform_internal(setup, input, output, (v4sd*)work, direction, 0);
}

void pffftd_transform_ordered(PFFFTD_Setup *setup, const double *input, double *output, double *work, pffftd_direction_t direction) {
  pffftd_transform_internal(setup, input, output, (v4sd*)work, direction, 1);
}
