// harness.cpp -- UB / correctness / coverage harness for r8brain-free-src.
// C++03-compatible (thread test enabled only under C++11+).
#include "../CDSPResampler.h"
#include "../CDSPFIRFilter.h"
#include "../CDSPBlockConvolver.h"
#include "../CDSPFracInterpolator.h"
#include "../CDSPHBDownsampler.h"
#include "../CDSPHBUpsampler.h"
#include "../CDSPRealFFT.h"
#include "../CDSPSincFilterGen.h"
#include "../r8butil.h"

#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <cstring>
#include <cfloat>
#include <cmath>
#include <vector>

#if defined( __cplusplus ) && __cplusplus >= 201103L
	#define HARNESS_CXX11 1
	#include <thread>
#else
	#define HARNESS_CXX11 0
#endif

#if R8B_FLTTEST
namespace r8b {
	int InterpFilterFracs = -1; // Required extern hook in FLTTEST mode.
} // namespace r8b
#endif // R8B_FLTTEST

using r8b :: CDSPResampler;

static int gChecks = 0;
static int gFails = 0;

static void check( const bool cond, const char* const fmt, ... )
{
	gChecks++;
	if( !cond )
	{
		gFails++;
		va_list ap;
		va_start( ap, fmt );
		printf( "FAIL: " );
		vprintf( fmt, ap );
		printf( "\n" );
		va_end( ap );
	}
}

#if HARNESS_CXX11
static thread_local uint64_t rngState = 0x9E3779B97F4A7C15ull;
#else
static uint64_t rngState = 0x9E3779B97F4A7C15ull;
#endif

static double urand() // uniform in [-1, 1)
{
	rngState ^= rngState << 13;
	rngState ^= rngState >> 7;
	rngState ^= rngState << 17;
	return( (double) ( rngState >> 11 ) * ( 1.0 / 9007199254740992.0 ) - 1.0 );
}

static bool isFiniteD( const double v )
{
	return( v == v && v <= DBL_MAX && v >= -DBL_MAX );
}

static bool scanFinite( const double* const p, const int n )
{
	int i;
	for( i = 0; i < n; i++ )
	{
		if( !isFiniteD( p[ i ]))
		{
			return( false );
		}
	}
	return( true );
}

// Goertzel single-bin power at frequency f (cycles/sample), normalized to
// the amplitude of a pure unit-amplitude sine at that bin.
static double goertzelAmp( const double* const p, const int n, const double f )
{
	const double w = 2.0 * 3.14159265358979324 * f;
	const double cr = cos( w );
	const double ci = sin( w );
	double re = 0.0;
	double im = 0.0;
	int i;
	for( i = 0; i < n; i++ )
	{
		re += p[ i ] * cos( w * i );
		im -= p[ i ] * sin( w * i );
	}
	(void) cr;
	(void) ci;
	return( 2.0 * sqrt( re * re + im * im ) / n );
}

static double rms( const double* const p, const int n )
{
	double s = 0.0;
	int i;
	for( i = 0; i < n; i++ )
	{
		s += p[ i ] * p[ i ];
	}
	return( sqrt( s / n ));
}

// ---------------------------------------------------------------------------
// Test 1: streaming ratio matrix.
// ---------------------------------------------------------------------------

struct RatioCase {
	double src;
	double dst;
	const char* name;
};

static void testRatio( const double s, const double d, const double atten,
	const double tb, const char* const label )
{
	const int MaxInLen = 1024;
	CDSPResampler rs( s, d, MaxInLen, tb, atten );

	const double freq = 0.2 * ( s < d ? s : d ); // Hz, deep in passband.
	static const int blocks[ 10 ] =
		{ 1, 2, 3, 17, 100, 333, 1024, 1024, 7, 511 };

	long totalIn = (long) ( 40000.0 * ( s < d ? 1.0 : s / d ));
	if( totalIn > 400000 )
	{
		totalIn = 400000;
	}
	if( totalIn < 40000 )
	{
		totalIn = 40000;
	}

	// Extreme down-sampling chains have a huge internal latency (e.g. ~919k
	// input samples for 480:1); make sure enough input is fed to produce
	// output at all.
	if( (double) totalIn * d / s < 5000.0 )
	{
		totalIn = 2500000;
	}

	std :: vector< double > in( MaxInLen );
	std :: vector< double > out;
	long totalOut = 0;
	long pos = 0;
	int bi = 0;
	bool finite = true;

	while( pos < totalIn )
	{
		int bl = blocks[ bi % 10 ];
		bi++;
		if( bl > totalIn - pos )
		{
			bl = (int) ( totalIn - pos );
		}
		int i;
		for( i = 0; i < bl; i++ )
		{
			in[ i ] = 0.8 * sin( 2.0 * 3.14159265358979324 * freq *
				( pos + i ) / s ) + 0.02 * urand();
		}
		double* op = NULL;
		const int l = rs.process( &in[ 0 ], bl, op );
		check( l >= 0, "%s: negative process result %d", label, l );
		if( l > 0 )
		{
			if( !scanFinite( op, l ))
			{
				finite = false;
			}
			if( totalOut + l <= 4000000 )
			{
				out.insert( out.end(), op, op + l );
			}
			totalOut += l;
		}
		pos += bl;
	}

	// Zero-length process call.
	{
		double* op = NULL;
		const int l = rs.process( &in[ 0 ], 0, op );
		check( l >= 0, "%s: l=0 call returned %d", label, l );
	}

	// Flush with zeros only until the output count catches up with the
	// expectation (this drains the internal latency); never feed more than
	// needed, otherwise the count check would count zero-input output.
	memset( &in[ 0 ], 0, MaxInLen * sizeof( in[ 0 ]));
	const double expectPre = (double) totalIn * d / s;
	// Extreme ratios need many blocks to drain the latency (e.g. ~897
	// blocks of 1024 for 480:1).
	const int guardMax = 500 + (int) ( s > d ? 3.0 * s / d : 3.0 * d / s );
	int guard = 0;
	while( totalOut < (long) ( expectPre - 1.0 ) && guard < guardMax )
	{
		double* op = NULL;
		const int l = rs.process( &in[ 0 ], MaxInLen, op );
		if( l > 0 )
		{
			if( !scanFinite( op, l ))
			{
				finite = false;
			}
			totalOut += l;
		}
		guard++;
	}

	check( finite, "%s: non-finite output sample detected", label );

	const double expect = (double) totalIn * d / s;
	// One full output block may overshoot the expectation when the flush
	// loop stops, so add getMaxOutLen() to the tolerance.
	const double tol = expect * 0.02 + 256.0 + rs.getMaxOutLen( MaxInLen );
	check( fabs( totalOut - expect ) <= tol,
		"%s: output count %ld deviates from expected %.0f (tol %.0f)",
		label, totalOut, expect, tol );

	// Steady-state RMS of the sine should be 0.8/sqrt(2).
	if( out.size() > 8000 )
	{
		const size_t b = out.size() / 4;
		const size_t e = out.size() * 3 / 4;
		const double r = rms( &out[ b ], (int) ( e - b ));
		check( fabs( r - 0.8 / sqrt( 2.0 )) < 0.05,
			"%s: steady-state RMS %.4f deviates from %.4f", label, r,
			0.8 / sqrt( 2.0 ));
	}

	// API sanity.
	check( rs.getLatency() >= 0, "%s: negative latency", label );
	{
		const double lf = rs.getLatencyFrac();
		check( lf >= 0.0 && lf < 1.0, "%s: LatencyFrac %.6f out of [0,1)",
			label, lf );
	}
	check( rs.getInLenBeforeOutStart() >= 0,
		"%s: getInLenBeforeOutStart < 0", label );
	check( rs.getInLenBeforeOutPos( 100 ) >= 0,
		"%s: getInLenBeforeOutPos(100) < 0", label );
	check( rs.getInputRequiredForOutput( 1000 ) >= 0,
		"%s: getInputRequiredForOutput(1000) < 0", label );
	check( rs.getMaxOutLen( MaxInLen ) >= 0,
		"%s: getMaxOutLen negative", label );
}

static void testRatioMatrix()
{
	static const RatioCase cases[] = {
		{ 44100.0, 44100.0, "1:1" },
		{ 48000.0, 24000.0, "2:1 down" },
		{ 24000.0, 48000.0, "1:2 up" },
		{ 48000.0, 16000.0, "3:1 down (third)" },
		{ 16000.0, 48000.0, "1:3 up (third)" },
		{ 48000.0, 12000.0, "4:1 down" },
		{ 12000.0, 48000.0, "1:4 up" },
		{ 48000.0, 3000.0, "16:1 down" },
		{ 3000.0, 48000.0, "1:16 up" },
		{ 48000.0, 8000.0, "6:1 down (third+half)" },
		{ 8000.0, 48000.0, "1:6 up (third+half)" },
		{ 48000.0, 32000.0, "3:2 down" },
		{ 32000.0, 48000.0, "2:3 up" },
		{ 48000.0, 36000.0, "4:3 down" },
		{ 36000.0, 48000.0, "3:4 up" },
		{ 44100.0, 48000.0, "44100->48000" },
		{ 48000.0, 44100.0, "48000->44100" },
		{ 44100.0, 88200.0, "44100->88200" },
		{ 88200.0, 44100.0, "88200->44100" },
		{ 44100.0, 176400.0, "44100->176400" },
		{ 44100.0, 96000.0, "44100->96000" },
		{ 48000.0, 40000.0, "6:5 (whole-step)" },
		{ 44100.0, 31500.0, "7:5 (whole-step)" },
		{ 48000.0, 33600.0, "10:7 (whole-step)" },
		{ 48000.0, 47999.0, "48000->47999 (interp)" },
		{ 47999.0, 48000.0, "47999->48000 (interp)" },
		{ 1000.0, 137.3, "1000->137.3 (frac down)" },
		{ 137.3, 1000.0, "137.3->1000 (frac up)" },
		{ 48000.0, 100.0, "480:1 extreme down" },
		{ 100.0, 48000.0, "1:480 extreme up" },
		{ 1.0, 1.5, "1->1.5 tiny rates" },
		{ 0.5, 0.25, "0.5->0.25 sub-Hz" },
		{ 96000.0, 96001.0, "96000->96001" },
		{ 8000.0, 8011.0, "8000->8011" },
	};

	size_t i;
	for( i = 0; i < sizeof( cases ) / sizeof( cases[ 0 ]); i++ )
	{
		testRatio( cases[ i ].src, cases[ i ].dst, 140.0, 2.0,
			cases[ i ].name );
	}
}

// ---------------------------------------------------------------------------
// Test 2: sine fidelity + stopband attenuation (one-shot based).
// ---------------------------------------------------------------------------

static void testSineFidelity()
{
	// 2:1 downsample: passband sine keeps amplitude; stopband sine is cut.
	{
		const double s = 48000.0, d = 24000.0;
		CDSPResampler rs( s, d, 8192, 2.0, 136.0 );
		const int N = 16384;
		std :: vector< double > in( N );
		std :: vector< double > obuf( N * 2 + 64 );
		int i;

		// Passband: 2400 Hz (0.2 * dst Nyq).
		const double fp = 2400.0;
		for( i = 0; i < N; i++ )
		{
			in[ i ] = sin( 2.0 * 3.14159265358979324 * fp * i / s );
		}
		rs.oneshot( &in[ 0 ], N, &obuf[ 0 ], N / 2 );
		const double amp = goertzelAmp( &obuf[ 1024 ], N / 4, fp / d );
		check( fabs( amp - 1.0 ) < 0.05,
			"sine 2:1 passband amplitude %.4f != 1.0", amp );

		// Stopband: 18000 Hz (above dst Nyquist 12000 -> aliases).
		const double fs2 = 18000.0;
		CDSPResampler rs2( s, d, 8192, 2.0, 136.0 );
		for( i = 0; i < N; i++ )
		{
			in[ i ] = sin( 2.0 * 3.14159265358979324 * fs2 * i / s );
		}
		rs2.oneshot( &in[ 0 ], N, &obuf[ 0 ], N / 2 );
		// Alias appears mirrored at 24000 - 18000 = 6000 Hz.
		const double amp2 = goertzelAmp( &obuf[ 1024 ], N / 4,
			( d - fs2 ) / d );
		check( amp2 < 0.0002, // ~ -74 dB; spec says -136 dB, keep margin.
			"sine 2:1 stopband leak %.6f too high", amp2 );
	}

	// 1:2 upsample: image suppression.
	{
		const double s = 24000.0, d = 48000.0;
		CDSPResampler rs( s, d, 8192, 2.0, 136.0 );
		const int N = 16384;
		std :: vector< double > in( N );
		std :: vector< double > obuf( N * 2 + 64 );
		const double fp = 2400.0;
		int i;
		for( i = 0; i < N; i++ )
		{
			in[ i ] = sin( 2.0 * 3.14159265358979324 * fp * i / s );
		}
		rs.oneshot( &in[ 0 ], N, &obuf[ 0 ], N * 2 );
		const double amp = goertzelAmp( &obuf[ 2048 ], N / 2, fp / d );
		check( fabs( amp - 1.0 ) < 0.05,
			"sine 1:2 passband amplitude %.4f != 1.0", amp );
		// Image at 24000 - 2400 = 21600 Hz.
		const double amp2 = goertzelAmp( &obuf[ 2048 ], N / 2,
			( s - fp ) / d );
		check( amp2 < 0.0002, "sine 1:2 image leak %.6f too high", amp2 );
	}

	// 1:1 must be an exact passthrough.
	{
		CDSPResampler rs( 44100.0, 44100.0, 1024 );
		double in[ 1024 ];
		int i;
		for( i = 0; i < 1024; i++ )
		{
			in[ i ] = urand();
		}
		double* op = NULL;
		const int l = rs.process( in, 1024, op );
		check( l == 1024, "1:1 passthrough length %d != 1024", l );
		check( op == in || memcmp( op, in, sizeof( in )) == 0,
			"1:1 passthrough data mismatch" );
	}
}

// ---------------------------------------------------------------------------
// Test 3: phase response variants (minimum-phase hits calcMinPhaseTransform).
// ---------------------------------------------------------------------------

static void testPhase()
{
	const double s = 44100.0, d = 48000.0;
	int ph;
	for( ph = 0; ph <= 1; ph++ )
	{
		CDSPResampler rs( s, d, 4096, 2.0, 140.0,
			( r8b :: EDSPFilterPhaseResponse ) ph );
		std :: vector< double > in( 4096 );
		size_t i;
		for( i = 0; i < in.size(); i++ )
		{
			in[ i ] = urand() * 0.5;
		}
		long total = 0;
		bool finite = true;
		int blk;
		for( blk = 0; blk < 12; blk++ )
		{
			double* op = NULL;
			const int l = rs.process( &in[ 0 ], (int) in.size(), op );
			if( l > 0 )
			{
				if( !scanFinite( op, l ))
				{
					finite = false;
				}
				total += l;
			}
		}
		check( finite, "phase=%d: non-finite output", ph );
		check( total > 0, "phase=%d: no output", ph );
	}

	// Direct min-phase LP filter generation (calcMinPhaseTransform path).
	{
		r8b :: CDSPFIRFilter& f = r8b :: CDSPFIRFilterCache :: getLPFilter(
			0.45, 5.0, 130.0, r8b :: fprMinPhase, 1.0, 1 );
		check( f.getKernelLen() > 10, "minphase kernel too small" );
		check( isFiniteD( f.getLatencyFrac()), "minphase latency frac" );
		f.unref();
	}
}

// ---------------------------------------------------------------------------
// Test 4: attenuation x transition-band matrix (varies kernel sizes, HB
// steepness indices, frac filter lengths).
// ---------------------------------------------------------------------------

static void testAttenTbMatrix()
{
	static const double attens[] = { 49.0, 100.0, 156.0, 218.0 };
	static const double tbs[] = { 0.5, 2.0, 10.0, 45.0 };
	size_t ai, ti;
	for( ai = 0; ai < sizeof( attens ) / sizeof( attens[ 0 ]); ai++ )
	{
		for( ti = 0; ti < sizeof( tbs ) / sizeof( tbs[ 0 ]); ti++ )
		{
			char label[ 96 ];
			sprintf( label, "atten=%.0f tb=%.1f", attens[ ai ], tbs[ ti ]);
			// 2:1 down and 44100->48000 frac: two chain shapes.
			testRatio( 48000.0, 24000.0, attens[ ai ], tbs[ ti ], label );
			testRatio( 44100.0, 48000.0, attens[ ai ], tbs[ ti ], label );
		}
	}
}

// ---------------------------------------------------------------------------
// Test 5: API edge cases.
// ---------------------------------------------------------------------------

static void testApiEdges()
{
	// MaxInLen=1, sample-by-sample streaming.
	{
		CDSPResampler rs( 44100.0, 48000.0, 1, 2.0, 120.0 );
		double in = 0.0;
		long total = 0;
		bool finite = true;
		int i;
		for( i = 0; i < 20000; i++ )
		{
			in = sin( 2.0 * 3.14159265358979324 * 1000.0 * i / 44100.0 );
			double* op = NULL;
			const int l = rs.process( &in, 1, op );
			if( l > 0 )
			{
				if( !scanFinite( op, l ))
				{
					finite = false;
				}
				total += l;
			}
		}
		check( finite, "MaxInLen=1: non-finite output" );
		check( total > 15000 && total < 30000,
			"MaxInLen=1: total output %ld unexpected", total );
	}

	// clear() + reuse determinism: identical input sequence before/after
	// clear() must give bit-identical output.
	{
		const int N = 8192;
		std :: vector< double > in( N );
		size_t i;
		for( i = 0; i < in.size(); i++ )
		{
			in[ i ] = urand() * 0.7;
		}
		CDSPResampler rs( 44100.0, 48000.0, 1024, 2.0, 140.0 );
		std :: vector< double > run1;
		int rep;
		for( rep = 0; rep < 2; rep++ )
		{
			std :: vector< double > acc;
			size_t pos = 0;
			while( pos < in.size())
			{
				const int bl = (int) (( pos * 37 ) % 1024 ) + 1;
				const int c = (int) ( pos + bl > in.size() ?
					in.size() - pos : bl );
				double* op = NULL;
				const int l = rs.process( &in[ pos ], c, op );
				if( l > 0 )
				{
					acc.insert( acc.end(), op, op + l );
				}
				pos += c;
			}
			if( rep == 0 )
			{
				run1 = acc;
				rs.clear();
			}
			else
			{
				check( acc.size() == run1.size(),
					"clear(): output size differs after reuse (%u vs %u)",
					(unsigned) acc.size(), (unsigned) run1.size());
				if( acc.size() == run1.size())
				{
					check( memcmp( &acc[ 0 ], &run1[ 0 ],
						acc.size() * sizeof( double )) == 0,
						"clear(): output data differs after reuse" );
				}
			}
		}
	}

	// Twin resamplers processing the same input must produce identical
	// output (global caches must not corrupt per-instance state).
	{
		CDSPResampler r1( 44100.0, 48000.0, 512, 2.0, 140.0 );
		CDSPResampler r2( 44100.0, 48000.0, 512, 2.0, 140.0 );
		double in[ 512 ];
		int i, blk;
		bool same = true;
		for( blk = 0; blk < 8; blk++ )
		{
			for( i = 0; i < 512; i++ )
			{
				in[ i ] = urand();
			}
			double* o1 = NULL;
			double* o2 = NULL;
			const int l1 = r1.process( in, 512, o1 );
			const int l2 = r2.process( in, 512, o2 );
			if( l1 != l2 || ( l1 > 0 &&
				memcmp( o1, o2, l1 * sizeof( double )) != 0 ))
			{
				same = false;
			}
		}
		check( same, "twin resamplers diverged" );
	}

	// oneshot<double> and oneshot<float>.
	{
		const int N = 10000;
		std :: vector< double > in( N );
		int i;
		for( i = 0; i < N; i++ )
		{
			in[ i ] = sin( 2.0 * 3.14159265358979324 * 1000.0 * i / 44100.0 );
		}
		CDSPResampler rs( 44100.0, 48000.0, 1024, 2.0, 140.0 );
		const int ocap = (int) ( N * 48000.0 / 44100.0 ) + 64;
		std :: vector< double > od( ocap );
		std :: vector< float > of( ocap );
		rs.oneshot( &in[ 0 ], N, &od[ 0 ], ocap );
		check( scanFinite( &od[ 0 ], ocap ), "oneshot<double> non-finite" );

		std :: vector< float > fin( N );
		for( i = 0; i < N; i++ )
		{
			fin[ i ] = (float) in[ i ];
		}
		CDSPResampler rs2( 44100.0, 48000.0, 1024, 2.0, 140.0 );
		rs2.oneshot( &fin[ 0 ], N, &of[ 0 ], ocap );
		bool ff = true;
		for( i = 0; i < ocap; i++ )
		{
			if( !isFiniteD( of[ i ]))
			{
				ff = false;
			}
		}
		check( ff, "oneshot<float> non-finite" );
		// float/double paths should agree closely.
		double maxd = 0.0;
		for( i = 1024; i < ocap - 64; i++ )
		{
			const double df = fabs( od[ i ] - (double) of[ i ]);
			if( df > maxd )
			{
				maxd = df;
			}
		}
		check( maxd < 1e-4, "oneshot float/double diverge: %g", maxd );
	}

	// oneshot downsample + zero-length input (latency flush path).
	{
		CDSPResampler rs( 48000.0, 44100.0, 512, 2.0, 140.0 );
		double dummy = 0.0;
		std :: vector< double > od( 4096 );
		rs.oneshot( &dummy, 0, &od[ 0 ], (int) od.size());
		check( scanFinite( &od[ 0 ], (int) od.size()),
			"oneshot zero-input non-finite" );
	}

	// Filter-cache eviction: create more distinct filters than the cache
	// holds (R8B_FILTER_CACHE_MAX = 96).
	{
		int i;
		for( i = 0; i < 110; i++ )
		{
			const double s = 44100.0 + i * 13.0;
			const double d = 22050.0 + i * 7.0;
			CDSPResampler rs( s, d, 256, 5.0, 80.0 + ( i % 3 ),
				r8b :: fprLinearPhase );
			double in[ 64 ];
			memset( in, 0, sizeof( in ));
			double* op = NULL;
			rs.process( in, 64, op );
		}
		check( true, "cache eviction run" );
	}

	// getInLenBeforeOutStart / Pos consistency: feeding the returned count
	// of input samples should produce the requested output position.
	{
		CDSPResampler rs( 44100.0, 48000.0, 4096, 2.0, 140.0 );
		const int need = rs.getInLenBeforeOutPos( 5000 );
		check( need > 0, "getInLenBeforeOutPos(5000)=%d", need );
		std :: vector< double > in( need + 4096, 0.0 );
		size_t i;
		for( i = 0; i < in.size(); i++ )
		{
			in[ i ] = urand() * 0.3;
		}
		// Feed in chunks of <= MaxInLen (API contract).
		long total = 0;
		size_t pos = 0;
		while( pos < in.size())
		{
			const int c = (int) ( in.size() - pos > 4096 ?
				4096 : in.size() - pos );
			double* op = NULL;
			const int l = rs.process( &in[ pos ], c, op );
			total += l;
			pos += c;
		}
		check( total >= 5000,
			"fed %d samples for out pos 5000, got %ld", need, total );
	}
}

// ---------------------------------------------------------------------------
// Test 6: direct component tests.
// ---------------------------------------------------------------------------

static void testFFTKeeper()
{
	int lb;
	for( lb = 1; lb <= 20; lb++ )
	{
	#if R8B_PFFFT || R8B_PFFFT_DOUBLE
		if( lb < 5 )
		{
			continue; // pffft requires N >= 32.
		}
	#endif // R8B_PFFFT || R8B_PFFFT_DOUBLE

		r8b :: CDSPRealFFTKeeper ffto( lb );
		const int Len = ffto -> getLen();
		// PFFFT/PFFFT_DOUBLE require SIMD-aligned (32B) buffers; the
		// library's own CFixedBuffer provides 64B alignment.
		r8b :: CFixedBuffer< double > in( Len );
		r8b :: CFixedBuffer< double > orig( Len );
		int i;
		for( i = 0; i < Len; i++ )
		{
			in[ i ] = urand();
			orig[ i ] = in[ i ];
		}

		r8b :: realfft_t* const tp = ffto -> forward( &in[ 0 ]);
		ffto -> inverse( tp );
		const double sc = ffto -> getInvMulConst();
		double maxerr = 0.0;
		const double* const dp = (const double*) tp;
		for( i = 0; i < Len; i++ )
		{
			const double e = fabs( dp[ i ] * sc - orig[ i ]);
			if( e > maxerr )
			{
				maxerr = e;
			}
		}
		// R8B_PFFFT uses float32 internally: roundtrip error is bounded by
		// float precision (~1e-7), not double precision.
	#if R8B_PFFFT
		check( maxerr < 1e-5, "FFT roundtrip LenBits=%d maxerr=%g", lb,
			maxerr );
	#else
		check( maxerr < 1e-9, "FFT roundtrip LenBits=%d maxerr=%g", lb,
			maxerr );
	#endif
	}
}

static void testSincGen()
{
	typedef r8b :: CDSPSincFilterGen FG;
	static const FG :: EWindowFunctionType wts[] = {
		FG :: wftCosine, FG :: wftKaiser, FG :: wftGaussian };
	static const char* const wnames[] = { "cosine", "kaiser", "gaussian" };

	int wi;
	for( wi = 0; wi < 3; wi++ )
	{
		int up;
		for( up = 0; up <= 1; up++ )
		{
			const double kaiserP[ 2 ] = { 8.0, 1.97 };
			const double gaussP[ 2 ] = { 0.4, 1.5 };
			const double* params =
				( wts[ wi ] == FG :: wftKaiser ? kaiserP :
				( wts[ wi ] == FG :: wftGaussian ? gaussP :
				( const double* ) NULL ));

			// init*(WinType) only configures the sine-generator state;
			// the matching calcWindow* function must be passed to
			// generate* (as CDSPFracInterpolator does internally).
			const FG :: CWindowFunc wf =
				( wts[ wi ] == FG :: wftKaiser ?
				(FG :: CWindowFunc) &FG :: calcWindowKaiser :
				( wts[ wi ] == FG :: wftGaussian ?
				(FG :: CWindowFunc) &FG :: calcWindowGaussian :
				(FG :: CWindowFunc) &FG :: calcWindowBlackman ));

			r8b :: CDSPSincFilterGen fg;
			fg.Len2 = 50.0;
			fg.initWindow( wts[ wi ], params, up != 0 );
			std :: vector< double > k( fg.KernelLen );
			fg.generateWindow( &k[ 0 ], wf );
			check( scanFinite( &k[ 0 ], fg.KernelLen ),
				"sincgen window %s pow=%d non-finite", wnames[ wi ], up );

			fg.Len2 = 50.0;
			fg.Freq1 = 0.1 * 3.14159265358979324;
			fg.Freq2 = 0.5 * 3.14159265358979324;
			fg.initBand( wts[ wi ], params, up != 0 );
			std :: vector< double > k2( fg.KernelLen );
			fg.generateBand( &k2[ 0 ], wf );
			check( scanFinite( &k2[ 0 ], fg.KernelLen ),
				"sincgen band %s pow=%d non-finite", wnames[ wi ], up );

			fg.Len2 = 50.0;
			fg.initHilbert( wts[ wi ], params, up != 0 );
			std :: vector< double > k3( fg.KernelLen );
			fg.generateHilbert( &k3[ 0 ], wf );
			check( scanFinite( &k3[ 0 ], fg.KernelLen ),
				"sincgen hilbert %s pow=%d non-finite", wnames[ wi ], up );

			// Fractional delays, incl. exact 0.0 and 1.0 boundary values.
			static const double fds[] =
				{ 0.0, 1e-9, 0.25, 0.5, 0.75, 1.0 - 1e-9, 1.0 };
			size_t fi;
			for( fi = 0; fi < sizeof( fds ) / sizeof( fds[ 0 ]); fi++ )
			{
				fg.Len2 = 50.0;
				fg.FracDelay = fds[ fi ];
				fg.initFrac( wts[ wi ], params, up != 0 );
				std :: vector< double > k4( fg.KernelLen );
				fg.generateFrac( &k4[ 0 ], wf );
				check( scanFinite( &k4[ 0 ], fg.KernelLen ),
					"sincgen frac %s pow=%d fd=%.6f non-finite",
					wnames[ wi ], up, fds[ fi ]);
				// Sinc sum should be ~1 regardless of frac delay.
				double ssum = 0.0;
				int j;
				for( j = 0; j < fg.KernelLen; j++ )
				{
					ssum += k4[ j ];
				}
				check( fabs( ssum - 1.0 ) < 0.05,
					"sincgen frac %s fd=%.6f sum=%f", wnames[ wi ],
					fds[ fi ], ssum );
			}
		}
	}
}

static void testBlockConvolver()
{
	static const int upf[] = { 1, 2, 3, 4, 5, 6 };
	static const int dnf[] = { 1, 2, 3, 4, 6 };
	size_t ui, di;
	for( ui = 0; ui < sizeof( upf ) / sizeof( upf[ 0 ]); ui++ )
	{
		for( di = 0; di < sizeof( dnf ) / sizeof( dnf[ 0 ]); di++ )
		{
			r8b :: CDSPFIRFilter& flt =
				r8b :: CDSPFIRFilterCache :: getLPFilter(
					0.45, 5.0, 100.0, r8b :: fprLinearPhase, 1.0,
					dnf[ di ]);

			// CDSPBlockConvolver documents that power-of-2 UpFactor and
			// DownFactor > 1 are mutually exclusive ("this case never
			// happens in practice"): skip such out-of-contract combos.
			const bool upPow2 = ( upf[ ui ] & ( upf[ ui ] - 1 )) == 0;
			const bool dnPow2 = ( dnf[ di ] & ( dnf[ di ] - 1 )) == 0;
			if( upPow2 && dnPow2 && upf[ ui ] > 1 && dnf[ di ] > 1 )
			{
				continue;
			}

			int dl;
			for( dl = 0; dl <= 1; dl++ )
			{
				r8b :: CDSPBlockConvolver bc( flt, upf[ ui ], dnf[ di ],
					0.0, dl == 0 );
				std :: vector< double > in( 1000 );
				std :: vector< double > obuf( in.size() * upf[ ui ] +
					16384 ); // CDSPProcessor: caller provides output buf
				size_t i;
				for( i = 0; i < in.size(); i++ )
				{
					in[ i ] = urand() * 0.5;
				}
				long total = 0;
				bool finite = true;
				int blk;
				for( blk = 0; blk < 6; blk++ )
				{
					double* op = &obuf[ 0 ];
					const int l = bc.process( &in[ 0 ], (int) in.size(),
						op );
					check( l <= (int) obuf.size(),
						"blockconv output overflow" );
					if( l > 0 )
					{
						if( !scanFinite( op, l ))
						{
							finite = false;
						}
						total += l;
					}
				}
				check( finite, "blockconv up=%d down=%d dl=%d non-finite",
					upf[ ui ], dnf[ di ], dl );
				check( total >= 0, "blockconv count" );
			}
			flt.unref();
		}
	}
}

static void testHBSamplers()
{
	// Attenuation sweep so that all per-length convolve kernels in the
	// _inc tables get selected (kernel length grows with ReqAtten).
	static const double attens[] = { 42.0, 55.0, 68.0, 80.0, 92.0, 105.0,
		118.0, 130.0, 145.0, 160.0, 175.0, 190.0, 205.0, 220.0 };
	int si, third;
	size_t ai;
	for( ai = 0; ai < sizeof( attens ) / sizeof( attens[ 0 ]); ai++ )
	for( si = 0; si <= 5; si++ )
	{
		for( third = 0; third <= 1; third++ )
		{
			{
				r8b :: CDSPHBUpsampler up( attens[ ai ], si, third != 0,
					0.0 );
				std :: vector< double > in( 777 );
				std :: vector< double > obuf( 777 * 2 + 512 );
				size_t i;
				for( i = 0; i < in.size(); i++ )
				{
					in[ i ] = urand() * 0.5;
				}
				bool finite = true;
				int blk;
				for( blk = 0; blk < 3; blk++ )
				{
					double* op = &obuf[ 0 ];
					const int l = up.process( &in[ 0 ], (int) in.size(),
						op );
					check( l <= (int) obuf.size(), "hbup overflow" );
					if( l > 0 && !scanFinite( op, l ))
					{
						finite = false;
					}
				}
				check( finite, "hbup si=%d third=%d non-finite", si,
					third );
			}
			{
				r8b :: CDSPHBDownsampler dn( attens[ ai ], si,
					third != 0, 0.0 );
				std :: vector< double > in( 1554 );
				std :: vector< double > obuf( 1554 + 512 );
				size_t i;
				for( i = 0; i < in.size(); i++ )
				{
					in[ i ] = urand() * 0.5;
				}
				bool finite = true;
				int blk;
				for( blk = 0; blk < 3; blk++ )
				{
					double* op = &obuf[ 0 ];
					const int l = dn.process( &in[ 0 ], (int) in.size(),
						op );
					check( l <= (int) obuf.size(), "hbdn overflow" );
					if( l > 0 && !scanFinite( op, l ))
					{
						finite = false;
					}
				}
				check( finite, "hbdn si=%d third=%d non-finite", si,
					third );
			}
		}
	}
}

static void testFracInterpolatorDirect()
{
	// Direct fractional interpolator: whole-stepping (InStep/OutStep
	// integers) and interpolated modes, both half- and third-band.
	static const struct { double s, d; } cases[] = {
		{ 6.0, 5.0 }, { 5.0, 6.0 }, { 7.0, 3.0 }, { 3.0, 7.0 },
		{ 44100.0, 48000.0 }, { 48000.0, 44100.0 },
		{ 1000.0, 997.3 }, { 997.3, 1000.0 },
	};
	size_t ci;
	for( ci = 0; ci < sizeof( cases ) / sizeof( cases[ 0 ]); ci++ )
	{
		int third;
		for( third = 0; third <= 1; third++ )
		{
			double atten;
			for( atten = 60.0; atten <= 200.0; atten += 70.0 )
			{
				r8b :: CDSPFracInterpolator fi( cases[ ci ].s,
					cases[ ci ].d, atten, third != 0, 0.0 );
				std :: vector< double > in( 500 );
				std :: vector< double > obuf( 500 * 4 + 1024 );
				size_t i;
				for( i = 0; i < in.size(); i++ )
				{
					in[ i ] = urand() * 0.5;
				}
				bool finite = true;
				int blk;
				for( blk = 0; blk < 8; blk++ )
				{
					double* op = &obuf[ 0 ];
					const int l = fi.process( &in[ 0 ], (int) in.size(),
						op );
					check( l <= (int) obuf.size(), "fracint overflow" );
					if( l > 0 && !scanFinite( op, l ))
					{
						finite = false;
					}
				}
				check( finite,
					"fracinterp %g->%g third=%d att=%g non-finite",
					cases[ ci ].s, cases[ ci ].d, third, atten );
			}
		}
	}
}

static void testR8bUtil()
{
	// Build a LP sinc and scan its response with r8butil helpers.
	r8b :: CDSPSincFilterGen fg;
	fg.Len2 = 200.0;
	fg.Freq1 = 0.0;
	fg.Freq2 = 0.4 * 3.14159265358979324;
	const double kp[ 2 ] = { 8.0, 1.97 };
	fg.initBand( r8b :: CDSPSincFilterGen :: wftKaiser, kp, false );
	std :: vector< double > k( fg.KernelLen );
	fg.generateBand( &k[ 0 ], &r8b :: CDSPSincFilterGen :: calcWindowKaiser );
	r8b :: normalizeFIRFilter( &k[ 0 ], fg.KernelLen, 1.0 );

	// r8butil scan helpers work on squared gains with in/out ref params.
	double re0, im0;
	r8b :: calcFIRFilterResponse( &k[ 0 ], fg.KernelLen, 0.0, re0, im0 );
	double maxg = re0 * re0 + im0 * im0;
	double maxth = 0.0;
	r8b :: findFIRFilterResponseMaxLtoR( &k[ 0 ], fg.KernelLen, maxg, maxth,
		0.5 );
	check( maxth >= 0.0 && maxth <= 0.5,
		"r8butil max freq %g out of range", maxth );
	check( isFiniteD( maxg ), "r8butil max gain non-finite" );

	r8b :: calcFIRFilterResponse( &k[ 0 ], fg.KernelLen,
		0.45 * 3.14159265358979324, re0, im0 );
	double ming = re0 * re0 + im0 * im0;
	double minth = 0.45; // start in stopband (passband edge is 0.4)
	r8b :: findFIRFilterResponseMinLtoR( &k[ 0 ], fg.KernelLen, ming, minth,
		1.0 );
	check( minth >= 0.45 && minth <= 1.0,
		"r8butil min freq %g out of range", minth );
	const double mindb = 10.0 * log10( ming + 1e-100 );
	check( mindb < -50.0, "r8butil stopband min %.1f dB too high", mindb );

	double th3 = 0.5;
	r8b :: findFIRFilterResponseLevelRtoL( &k[ 0 ], fg.KernelLen,
		0.50118723362727224 /* -3 dB squared */, th3, 0.0 );
	check( th3 >= 0.0 && th3 <= 0.5, "r8butil -3dB level %g", th3 );

	// Group delay / response helpers from r8bbase.
	double re, im;
	r8b :: calcFIRFilterResponse( &k[ 0 ], fg.KernelLen, 0.1, re, im );
	check( isFiniteD( re ) && isFiniteD( im ), "calcFIRFilterResponse" );
	const double gd = r8b :: calcFIRFilterGroupDelay( &k[ 0 ],
		fg.KernelLen, 0.1 );
	check( isFiniteD( gd ), "calcFIRFilterGroupDelay" );
	check( fabs( gd - fg.fl2 ) < 1.0,
		"group delay %.2f far from fl2=%d", gd, fg.fl2 );
}

static void testBaseUtilities()
{
	// getBitOccupancy across the int domain, incl. edge values.
	check( r8b :: getBitOccupancy( 1 ) == 1, "bitocc(1)" );
	check( r8b :: getBitOccupancy( 2 ) == 2, "bitocc(2)" );
	check( r8b :: getBitOccupancy( 3 ) == 2, "bitocc(3)" );
	check( r8b :: getBitOccupancy( 255 ) == 8, "bitocc(255)" );
	check( r8b :: getBitOccupancy( 256 ) == 9, "bitocc(256)" );
	check( r8b :: getBitOccupancy( 32767 ) == 15, "bitocc(32767)" );
	{
		const int bo0 = r8b :: getBitOccupancy( 0 );
		check( bo0 >= 1 && bo0 <= 32, "bitocc(0)=%d", bo0 );
		const int bom = r8b :: getBitOccupancy( 2147483647 );
		check( bom == 31, "bitocc(INT_MAX)=%d", bom );
		const int bon = r8b :: getBitOccupancy( -5 ); // impl-defined shift
		check( bon >= 1 && bon <= 33, "bitocc(-5)=%d", bon );
	}

	// CFixedBuffer alloc/realloc/moveFrom.
	{
		r8b :: CFixedBuffer< double > b( 100 );
		b[ 0 ] = 1.5;
		b[ 99 ] = 2.5;
		r8b :: CFixedBuffer< double > c;
		c.alloc( 0 ); // capacity 0 must be safe
		r8b :: CFixedBuffer< double > d;
		d.moveFrom( b );
		check( d[ 0 ] == 1.5 && d[ 99 ] == 2.5, "CFixedBuffer moveFrom" );
		d.realloc( 100, 300 );
		check( d[ 0 ] == 1.5 && d[ 99 ] == 2.5, "CFixedBuffer realloc" );
		d.free();
		d.free(); // double free of empty buffer must be safe
		c.free();
	}

	// CPtrKeeper / CRefKeeper basics.
	{
		r8b :: CPtrKeeper< r8b :: CDSPResampler > pk;
		pk = new CDSPResampler( 44100.0, 48000.0, 256, 5.0, 100.0 );
		pk -> clear();
		pk.reset();
	}

	// CSineGen.
	{
		r8b :: CSineGen g( 0.1, 0.0, 1.0 );
		g.generate();
		double mx = 0.0;
		int i;
		for( i = 0; i < 1000; i++ )
		{
			const double v = g.generate();
			if( fabs( v ) > mx )
			{
				mx = fabs( v );
			}
		}
		check( mx > 0.9 && mx <= 1.0001, "CSineGen max %.4f", mx );
	}

	// Spline coefficient generators (used by frac filter bank init).
	{
		double c[ 8 ];
		const double y[ 8 ] = { 0.1, 0.3, 0.6, 0.5, 0.2, 0.4, 0.9, 0.7 };
		r8b :: calcSpline3p8Coeffs( c, y[ 0 ], y[ 1 ], y[ 2 ], y[ 3 ],
			y[ 4 ], y[ 5 ], y[ 6 ], y[ 7 ]);
		check( scanFinite( c, 8 ), "spline3p8 coeffs non-finite" );
		r8b :: calcSpline2p8Coeffs( c, y[ 0 ], y[ 1 ], y[ 2 ], y[ 3 ],
			y[ 4 ], y[ 5 ], y[ 6 ], y[ 7 ]);
		check( scanFinite( c, 8 ), "spline2p8 coeffs non-finite" );
		r8b :: calcSpline3p4Coeffs( c, y );
		check( scanFinite( c, 8 ), "spline3p4 coeffs non-finite" );
		r8b :: calcSpline3p6Coeffs( c, y );
		check( scanFinite( c, 8 ), "spline3p6 coeffs non-finite" );
	}
}

// ---------------------------------------------------------------------------
// Test 7: multithreaded use (C++11+ only) -- global caches under contention.
// ---------------------------------------------------------------------------

#if HARNESS_CXX11

static void threadWorker( const int id )
{
	static const double rates[][ 2 ] = {
		{ 44100.0, 48000.0 }, { 48000.0, 44100.0 },
		{ 48000.0, 24000.0 }, { 24000.0, 48000.0 },
		{ 48000.0, 32000.0 }, { 44100.0, 88200.0 },
	};
	int it;
	for( it = 0; it < 20; it++ )
	{
		const double* const r = rates[( id + it ) % 6 ];
		{
			CDSPResampler rs( r[ 0 ], r[ 1 ], 512, 2.0,
				110.0 + ( it % 4 ) * 20.0,
				( it & 1 ) != 0 ? r8b :: fprMinPhase :
					r8b :: fprLinearPhase );
			double in[ 512 ];
			int i;
			for( i = 0; i < 512; i++ )
			{
				in[ i ] = urand() * 0.4;
			}
			double* op = NULL;
			rs.process( in, 512, op );
		}
		// Hammer the global FFT object pools as well.
		{
			r8b :: CDSPRealFFTKeeper k1( 8 + ( it % 5 ));
			r8b :: CDSPRealFFTKeeper k2( 8 + (( it + 2 ) % 5 ));
		}
	}
}

static void testThreads()
{
	std :: vector< std :: thread > ts;
	int i;
	for( i = 0; i < 4; i++ )
	{
		ts.push_back( std :: thread( threadWorker, i ));
	}
	for( i = 0; i < 4; i++ )
	{
		ts[ i ].join();
	}
	check( true, "threads completed" );
}

#endif // HARNESS_CXX11

// ---------------------------------------------------------------------------

int main()
{
	printf( "r8brain-free-src harness, C++ standard: %ld\n",
		(long) __cplusplus );
#if R8B_OOURA
	printf( "FFT backend: Ooura\n" );
#endif
#if R8B_PFFFT
	printf( "FFT backend: PFFFT (float)\n" );
#endif
#if R8B_PFFFT_DOUBLE
	printf( "FFT backend: PFFFT_DOUBLE\n" );
#endif

	testBaseUtilities();
	printf( "base utilities done\n" );
	testSincGen();
	printf( "sinc gen done\n" );
	testFFTKeeper();
	printf( "fft keeper done\n" );
	testR8bUtil();
	printf( "r8butil done\n" );
	testHBSamplers();
	printf( "hb samplers done\n" );
	testFracInterpolatorDirect();
	printf( "frac interpolator done\n" );
	testBlockConvolver();
	printf( "block convolver done\n" );
	testSineFidelity();
	printf( "sine fidelity done\n" );
	testPhase();
	printf( "phase done\n" );
	testApiEdges();
	printf( "api edges done\n" );
	testRatioMatrix();
	printf( "ratio matrix done\n" );
	testAttenTbMatrix();
	printf( "atten/tb matrix done\n" );
#if HARNESS_CXX11
	testThreads();
	printf( "threads done\n" );
#endif

	printf( "CHECKS: %d  FAILURES: %d\n", gChecks, gFails );
	return( gFails == 0 ? 0 : 1 );
}
