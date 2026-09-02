# r8brain-free-src v7.5 — Undefined Behavior, C++ Compliance & Coverage Evaluation

**Code under test:** r8brain-free-src 7.5 (`R8B_VERSION "7.5"`, © 2013–2026 Aleksey Vaneev) — 28 uploaded files, ~23,000 LOC total: the header-only r8brain core (`CDSP*.h`, `r8b*.h`, two `_inc.h` half-band kernel files), the bundled Ooura FFT (`fft4g.h`), and the vendored PFFFT library (float `pffft.c`, double `pffft_double.c`, SIMD macro headers `pf_*.h`).

**Toolchain:** GCC 12.2.0 (Debian 12.2.0-14), x86-64 Linux. Sanitizers: AddressSanitizer + UndefinedBehaviorSanitizer (`-fno-sanitize-recover=all`) and ThreadSanitizer. Coverage: gcov. Static analysis: `-fanalyzer`, `-Wall -Wextra -Wpedantic`.

**Method note:** the uploaded `*_inc.h` files were renamed back to their original `.inc` names for the build (they are `#include`d by name from the parent headers). No external code was fetched; `example.cpp` cannot be built standalone because it depends on the library's unbundled `CWaveFile` helper class — everything else was exercised directly.

---

## 1. Test harness

A purpose-built ~1,200-line harness (`r8brain_test_harness.cpp`, delivered alongside this report) drives **26,245–26,249 assertions per run** across these groups:

| Test group | What it exercises |
|---|---|
| Base utilities | `CFixedBuffer`, `align_ptr`, `getBitOccupancy`, refcount/keeper classes |
| Sinc filter gen | window/band/Hilbert/fractional-delay generation, all 3 window types × power on/off, frac delays incl. exact 0.0 and 1.0 |
| FFT keeper | round-trip forward/inverse for LenBits 1–20 (5–20 for PFFFT), max error bounds |
| r8butil | filter-response scans, −3 dB level search, group delay |
| HB up/downsamplers | steepness indices 0–5 × third-band on/off × **attenuation sweep 42–220 dB** (selects every per-length SSE2 convolve kernel in the `.inc` files) |
| Frac interpolator (direct) | whole-stepping and interpolated modes, 8 rate pairs × third on/off × 3 attenuations |
| Block convolver (direct) | up factors {1,2,3,4,5,6} × down factors {1,2,3,4,6} × latency-consume on/off |
| Sine fidelity | passband amplitude/RMS and stopband attenuation measurements |
| Phase response | linear- vs minimum-phase chains, min-phase filter generation |
| API edges | `MaxInLen=1` sample-by-sample streaming, zero-length `process()`, `clear()`+reuse determinism, `oneshot()`, `getInLenBeforeOutPos()` consistency |
| Ratio matrix | 34 sample-rate ratios incl. 1:1, 2:1/1:2 … 16:1, third-band chains, whole-stepping ratios, fractional ratios, 480:1 / 1:480 extreme ratios, sub-Hz rates |
| Atten × transition-band matrix | 4 attenuations × 4 transition widths × 2 chain shapes |
| Threads (C++11+) | 4 concurrent threads, each with own resampler, hammering the global FFT/filter caches |

Input is randomized (xorshift PRNG) plus deterministic sine sweeps; every output block is scanned for NaN/Inf; streaming totals are checked against the theoretical `totalIn · dst/src` count.

---

## 2. Dynamic analysis results — no defects in public-API usage

| Configuration | Std | Checks | Failures | Sanitizer findings |
|---|---|---|---|---|
| Ooura FFT (default) | C++17 | 26,249 | 0 | none |
| Ooura FFT | **C++03** | 26,248 | 0 | none |
| `R8B_PFFFT=1` (float32 FFT) | C++17 | 26,245 | 0 | none |
| `R8B_PFFFT_DOUBLE=1` | C++17 | 26,245 | 0 | none |
| `R8B_EXTFFT=1`, `R8B_FLTTEST=1` (coverage runs) | C++17 | 26,249 | 0 | — |
| Ooura, ThreadSanitizer | C++17 | 26,249 | 0 | **no data races** |

- **ASan+UBSan** (`-fsanitize=address,undefined`, leak detection on): zero heap/stack/global overflows, zero use-after-free, zero leaks, zero signed-overflow/shift/alignment/null violations across all four build configurations.
- **TSan**: zero races — the README's "the code is thread-safe [with] a separate resampler object … per stream" claim holds for the tested scenario, including the shared global FFT/filter caches (mutex-protected).
- Extreme configurations (1:480 upsample producing 19.2 M output samples; 480:1 downsample with ~918 k-sample internal latency; `MaxInLen=1` streaming) all behave correctly.

---

## 3. Findings (library side)

No undefined behavior and no memory-safety defect was found in any usage reachable through the public `CDSPResampler` API. The following are **latent robustness gaps**, all in internal machinery or out-of-contract usage:

### 3.1 Internal `CDSPBlockConvolver`: unguarded-by-default precondition → heap overflow under PFFFT backends (Low)

`CDSPBlockConvolver` requires that power-of-2 `UpFactor > 1` and power-of-2 `DownFactor > 1` never co-occur (the code comments say "This case never happens in practice due to mutual exclusion"). If a direct user constructs e.g. up=4/down=4 anyway:

- **Ooura build:** silently works (spectrum lives in the caller-sized double buffer).
- **`R8B_PFFFT` build:** heap-buffer-overflow — `mirrorInputSpectrum()` expands the spectrum to `BlockLen2` elements, but `process()` handed it the smaller `fftin` work buffer (both keepers share the same FFT length here). Reproduced under ASan: 2,048-byte `memcpy` past a 512-float buffer at `CDSPBlockConvolver.h:630`.
- **`R8B_PFFFT_DOUBLE` build:** same overflow at `CDSPBlockConvolver.h:619`.

The precondition *is* asserted (`R8BASSERT(UpShift == 0)` at `CDSPBlockConvolver.h:120`, verified to trap when enabled), but `R8BASSERT` is a **no-op by default** (`r8bconf.h`). Not reachable via `CDSPResampler` (its `CommonRatios` are reduced fractions {1,2},{1,3},{2,3},{3,2},{3,4}; all internal convolvers use pure-up or pure-down factors).

### 3.2 `CDSPRealFFTKeeper` direct use requires SIMD-aligned buffers — undocumented (Low)

With `R8B_PFFFT_DOUBLE`, `CDSPRealFFT::forward/inverse` pass the caller's buffer straight into `pffftd_transform`, which requires **32-byte alignment** (`assert(VALIGNED(...))`; misaligned SSE2 `movapd` faults in release builds). The library's own buffers are 64-byte-aligned via `CFixedBuffer`, so internal use is safe; but a direct user passing a `std::vector<double>` (16-byte-aligned) buffer trips the assert — observed non-deterministically (allocation-layout dependent). Not documented on the public-in-header `CDSPRealFFTKeeper` interface.

### 3.3 No runtime guard against `process()` input length > `MaxInLen` (Informational)

`CDSPResampler::process()` documents that `l` "should not exceed the `MaxInLen` supplied in the constructor", but only asserts `l >= 0`. Violating the contract produces a genuine heap overflow deep in `CDSPBlockConvolver::copyToOutput` (reproduced under ASan during harness development). This is documented API contract, so it is not a defect; adding `R8BASSERT(l <= MaxInLen)` to `CDSPResampler::process` and intermediate processors would turn silent corruption into a debuggable trap.

### 3.4 `CDSPSincFilterGen`: init/generate window pairing is a silent-mismatch footgun (Design note)

`initBand/initFrac(wftKaiser, …)` configures generator state but `generateBand/generateFrac` default their `wfunc` argument to `calcWindowBlackman`. Calling them without the matching `calcWindowKaiser/calcWindowGaussian` pointer silently yields a wrong (near-zero-energy) filter — wrong results, not UB. The library's internal caller (`CDSPFracInterpolator.h:112`) pairs them correctly. A stored accessor would remove the trap.

### 3.5 Aggressive-but-isolated type punning in PFFFT float path (Note)

`CDSPRealFFT::forward` (R8B_PFFFT) reuses a live `double` buffer as `float` scratch via placement-new (`construct_ptr<float>`) and later restores doubles the same way. Formally this rides the edge of the C++ object-lifetime rules (fully sanctioned only in C++20 implicit-lifetime terms); it is the standard DSP-library pattern, isolated to two functions, and works on all mainstream compilers. No action needed; worth knowing if the code is ever ported to an exotic toolchain.

---

## 4. C++ compliance

| Standard | `-Wall -Wextra -Wpedantic` result |
|---|---|
| C++03 | one warning: `r8bconf.h:41` — anonymous variadic macro `R8BCONSOLE(...)` (variadic macros are C++11). Accepted by GCC/Clang/MSVC as an extension; strict `-pedantic-errors` C++03 builds must predefine `R8BCONSOLE` (the file is the designated user-editable config header) |
| C++11 | clean |
| C++17 | clean |
| C++20 | clean |

- No C++17-removed constructs (`register`, dynamic exception specs, `auto_ptr`, …) anywhere.
- All 14 core headers are individually self-contained (compile standalone).
- `-m32` syntax check passes.
- Global caches use function-local statics (thread-safe magic statics in C++11+; the README already conditions pre-C++11 use on linking a threading lib).
- `CDSPProcessor` has a virtual destructor — polymorphic deletion through the base is safe.
- `CFixedBuffer` provides 64-byte alignment and correct `realloc` clamping; allocator frees via the raw pointer (`Data0`), placement-new is used over `char[]` storage — lifetime-correct for the trivially-destructible element types used.

### Static analysis (`g++ -fanalyzer`)

15 warnings, **all false positives**: 12 × `-Wanalyzer-malloc-leak` (analyzer loses pointers escaping into `CFixedBuffer` member fields; ASan's leak checker reports zero leaks at exit), 3 × `-Wanalyzer-possible-null-dereference/-argument` (analyzer models `operator new` as nullable; r8b uses throwing `new`, never null).

---

## 5. Coverage (gcov, union over 5 build configs)

Configs: Ooura, `R8B_PFFFT`, `R8B_PFFFT_DOUBLE`, `R8B_EXTFFT`, `R8B_FLTTEST`. A line counts as covered if executed in any config.

| File | Covered / Total | % |
|---|---|---|
| CDSPResampler.h | 197 / 202 | 97.5 |
| CDSPBlockConvolver.h | 249 / 251 | 99.2 |
| CDSPFIRFilter.h | 228 / 229 | 99.6 |
| CDSPFracInterpolator.h | 327 / 371 | 88.1 |
| CDSPHBDownsampler.h | 92 / 100 | 92.0 |
| CDSPHBUpsampler.h | 281 / 298 | 94.3 |
| CDSPHBDownsampler.inc / CDSPHBUpsampler.inc (SSE2 kernel tables) | 317 / 467 | **100 / 100** |
| CDSPProcessor.h | 5 / 5 | 100.0 |
| CDSPRealFFT.h | 223 / 233 | 95.7 |
| CDSPSincFilterGen.h | 226 / 230 | 98.3 |
| r8bbase.h | 228 / 244 | 93.4 |
| r8butil.h | 102 / 105 | 97.1 |
| fft4g.h (Ooura) | 507 / 507 | 100.0 |
| **r8brain core total** | **3,132 / 3,242** | **96.6** |
| pffft.c / pffft_double.c / pffft_priv_impl.h (vendored) | 1,084 / 2,607 | 41.6 |

**Uncovered core lines** are confined to: latency-consumption early-return branches that require non-zero `PrevLatency` chains; the `ElementSize == 2/4` filter-bank table layouts (not instantiated by any current caller); small-FFT special-case loops in `CDSPRealFFT.h`; `getBitOccupancy` high-word branch; the `R8B_FLTTEST` fraction-override assignment; and trivial `getLatency()` bodies on the `DoConsumeLatency=true` path. No uncovered line is on a data-processing hot path.

**PFFFT's 41.6%** reflects that CDSPRealFFT uses only PFFFT's *real, sequential-order* transform entry points (setup/destroy/transform/zreorder — all exercised); the uncovered remainder is unused vendored features (complex/interleaved transforms, unused radix specializations, scalar fallback code compiled but not selected).

---

## 6. Verdict

**The library is clean for its documented public API.** Across ~105k cumulative assertions on four sanitizer configurations spanning C++03–C++20, three FFT backends and extreme resampling ratios, no undefined behavior, memory error, leak, or data race was observed. Coverage of the core library reaches 96.6% (100% for the SSE2 kernel tables and the Ooura FFT), so this is not a shallow result.

The only genuine hazards found sit at internal/out-of-contract boundaries: the internal block convolver's power-of-2 up/down mutual-exclusion precondition (§3.1) and the FFT keeper's alignment requirement (§3.2) are both enforced only by assertions that are compiled out by default. Enabling `R8BASSERT` in debug builds (one line in `r8bconf.h`) mitigates both.

---

## Appendix — harness bugs found and fixed during testing (not library defects)

For transparency, these false alarms were root-caused to the test code, and each is a plausible user error worth documenting:

1. Feeding 4,096 samples to a resampler constructed with `MaxInLen=1024` → heap overflow (the §3.3 contract violation, caught by ASan).
2. Flush loop that kept zero-input blocks flowing past the drain point, corrupting output-count expectations.
3. `generateFrac/generateBand` called without the matching `calcWindow*` function after `init*(wftKaiser/wftGaussian)` (the §3.4 footgun) → near-zero filter sums and a −40 dB "stopband" measurement.
4. `findFIRFilterResponseMinLtoR` seeded inside the passband (contract requires a monotonically-decreasing region) and scanned only to 0.5 instead of 1.0.
5. `std::vector<double>` buffer passed to `CDSPRealFFTKeeper::forward` under PFFFT_DOUBLE → alignment assert (§3.2).
6. Shared PRNG state across harness threads → TSan race (harness-side only).
7. FFT round-trip tolerance 1e-9 applied to the float32 PFFFT backend → 16 false failures at ~1e-7; relaxed to 1e-5 for that backend.
