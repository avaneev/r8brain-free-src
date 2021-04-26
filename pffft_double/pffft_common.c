
#include "pffft.h"

#include <stdlib.h>

/* SSE and co like 16-bytes aligned pointers
 * with a 64-byte alignment, we are even aligned on L2 cache lines... */
#define MALLOC_V4SF_ALIGNMENT 64

static void * Valigned_malloc(size_t nb_bytes) {
  void *p, *p0 = malloc(nb_bytes + MALLOC_V4SF_ALIGNMENT);
  if (!p0) return (void *) 0;
  p = (void *) (((size_t) p0 + MALLOC_V4SF_ALIGNMENT) & (~((size_t) (MALLOC_V4SF_ALIGNMENT-1))));
  *((void **) p - 1) = p0;
  return p;
}

static void Valigned_free(void *p) {
  if (p) free(*((void **) p - 1));
}


static int next_power_of_two(int N) {
  /* https://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2 */
  /* compute the next highest power of 2 of 32-bit v */
  unsigned v = N;
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v++;
  return v;
}

static int is_power_of_two(int N) {
  /* https://graphics.stanford.edu/~seander/bithacks.html#DetermineIfPowerOf2 */
  int f = N && !(N & (N - 1));
  return f;
}

static int min_fft_size(pffft_transform_t transform) {
  /* unfortunately, the fft size must be a multiple of 16 for complex FFTs
     and 32 for real FFTs -- a lot of stuff would need to be rewritten to
     handle other cases (or maybe just switch to a scalar fft, I don't know..) */
  int simdSz = pffft_simd_size();
  if (transform == PFFFT_REAL)
    return ( 2 * simdSz * simdSz );
  else if (transform == PFFFT_COMPLEX)
    return ( simdSz * simdSz );
  else
    return 1;
}


void *pffft_aligned_malloc(size_t nb_bytes) { return Valigned_malloc(nb_bytes); }
void pffft_aligned_free(void *p) { Valigned_free(p); }
int pffft_next_power_of_two(int N) { return next_power_of_two(N); }
int pffft_is_power_of_two(int N) { return is_power_of_two(N); }
int pffft_min_fft_size(pffft_transform_t transform) { return min_fft_size(transform); }

void *pffftd_aligned_malloc(size_t nb_bytes) { return Valigned_malloc(nb_bytes); }
void pffftd_aligned_free(void *p) { Valigned_free(p); }
int pffftd_next_power_of_two(int N) { return next_power_of_two(N); }
int pffftd_is_power_of_two(int N) { return is_power_of_two(N); }
int pffftd_min_fft_size(pffft_transform_t transform) { return min_fft_size(transform); }

