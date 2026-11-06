
/* Copyright (c) 2026  Michael Neuling ( mikey@neuling.org )

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
*/

#ifndef PF_RVV_DBL_H
#define PF_RVV_DBL_H

/*
  RISC-V Vector (RVV 1.0) support macros for double precision.

  pffft hard-codes SIMD_SZ=4 for double precision, matching AVX 256-bit
  vectors. We use 256-bit fixed-length vectors via GCC's vector_size
  attribute. With VLEN >= 256 the compiler lowers these to single RVV
  instructions with VL=4, e64, m1 (e.g. the Spacemit X-class cores); with
  VLEN=128 it splits each operation in two.
*/
#if !defined(SIMD_SZ) && !defined(PFFFT_SIMD_DISABLE) && defined(PFFFT_ENABLE_RVV) && defined(__riscv_vector) && (__riscv_v_elen >= 64)
#pragma message( __FILE__ ": RISC-V Vector (RVV) double macros are defined" )

typedef double   v4sf __attribute__((vector_size(32), aligned(8)));
typedef long int v4si __attribute__((vector_size(32), aligned(8)));

#  define SIMD_SZ 4

typedef union v4sf_union {
  v4sf   v;
  double f[SIMD_SZ];
} v4sf_union;

/* Both clang and recent GCC accept __builtin_shufflevector; older GCC
 * only ships __builtin_shuffle. Pick whichever is available. GCC < 10 has
 * no __has_builtin at all, so it must not appear in the same #if as the
 * defined() check. */
#if defined(__has_builtin)
#  if __has_builtin(__builtin_shufflevector)
#    define PF_RVV_HAVE_SHUFFLEVECTOR 1
#  endif
#endif
#ifdef PF_RVV_HAVE_SHUFFLEVECTOR
#  define PF_RVV_SHUF4(a, b, i0, i1, i2, i3) \
     __builtin_shufflevector((a), (b), i0, i1, i2, i3)
#else
#  define PF_RVV_SHUF4(a, b, i0, i1, i2, i3) \
     __builtin_shuffle((a), (b), (v4si){i0, i1, i2, i3})
#endif

#  define VARCH "RVV"
#  define VREQUIRES_ALIGN 0
#  define VZERO()       ((v4sf){0.0, 0.0, 0.0, 0.0})
#  define VMUL(a,b)     ((a) * (b))
#  define VADD(a,b)     ((a) + (b))
#  define VSUB(a,b)     ((a) - (b))
#  define VMADD(a,b,c)  ((a) * (b) + (c))
#  define VMSUB(a,b,c)  ((c) - (a) * (b))
#  define LD_PS1(p)     ((v4sf){(p), (p), (p), (p)})
#  define VLOAD_UNALIGNED(ptr)  (*(const v4sf *)(ptr))
#  define VLOAD_ALIGNED(ptr)    (*(const v4sf *)(ptr))

/* INTERLEAVE2: {a0,a1,a2,a3},{b0,b1,b2,b3} -> {a0,b0,a1,b1},{a2,b2,a3,b3} */
#  define INTERLEAVE2(in1, in2, out1, out2) do {                        \
    v4sf in1__ = (in1), in2__ = (in2);                                  \
    out1 = PF_RVV_SHUF4(in1__, in2__, 0, 4, 1, 5);                      \
    out2 = PF_RVV_SHUF4(in1__, in2__, 2, 6, 3, 7);                      \
  } while (0)

/* UNINTERLEAVE2: {a0,a1,a2,a3},{b0,b1,b2,b3} -> {a0,a2,b0,b2},{a1,a3,b1,b3} */
#  define UNINTERLEAVE2(in1, in2, out1, out2) do {                      \
    v4sf in1__ = (in1), in2__ = (in2);                                  \
    out1 = PF_RVV_SHUF4(in1__, in2__, 0, 2, 4, 6);                      \
    out2 = PF_RVV_SHUF4(in1__, in2__, 1, 3, 5, 7);                      \
  } while (0)

#  define VTRANSPOSE4(x0, x1, x2, x3) do {                              \
    v4sf r0__ = (x0), r1__ = (x1), r2__ = (x2), r3__ = (x3);            \
    v4sf t0__ = PF_RVV_SHUF4(r0__, r1__, 0, 4, 1, 5);                   \
    v4sf t1__ = PF_RVV_SHUF4(r0__, r1__, 2, 6, 3, 7);                   \
    v4sf t2__ = PF_RVV_SHUF4(r2__, r3__, 0, 4, 1, 5);                   \
    v4sf t3__ = PF_RVV_SHUF4(r2__, r3__, 2, 6, 3, 7);                   \
    x0 = PF_RVV_SHUF4(t0__, t2__, 0, 1, 4, 5);                          \
    x1 = PF_RVV_SHUF4(t0__, t2__, 2, 3, 6, 7);                          \
    x2 = PF_RVV_SHUF4(t1__, t3__, 0, 1, 4, 5);                          \
    x3 = PF_RVV_SHUF4(t1__, t3__, 2, 3, 6, 7);                          \
  } while (0)

/* VSWAPHL(a,b) -> {b0, b1, a2, a3} */
#  define VSWAPHL(a, b)  PF_RVV_SHUF4((b), (a), 0, 1, 6, 7)

/* reverse/flip all doubles */
#  define VREV_S(a)      PF_RVV_SHUF4((a), (a), 3, 2, 1, 0)
/* reverse/flip complex doubles */
#  define VREV_C(a)      PF_RVV_SHUF4((a), (a), 2, 3, 0, 1)

/* Both load paths happen to be safe for any 8-byte aligned pointer
 * (the typedef requests only aligned(8)). We still report v4sf-sized
 * alignment here so the --test-simd self-test sees a mix of aligned
 * and unaligned addresses, matching the SSE2/AVX double convention. */
#  define VALIGNED(ptr) ((((uintptr_t)(ptr)) & 0x1F) == 0)

#endif

#endif /* PF_RVV_DBL_H */
