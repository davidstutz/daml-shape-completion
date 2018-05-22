#ifndef REMESHER_DEFINES_HEADER
#define REMESHER_DEFINES_HEADER

/*
 * Library compile-time configuration.
 */

/* Checks for NAN values. This is not very expensive. */
#define REMESHER_NAN_CHECKS 1
#if REMESHER_NAN_CHECKS
#   ifdef WIN32
namespace std
{
    template <typename T>
    inline bool isnan (T const& x)
    { return (x != x); }
}
#   endif
#   define REMESHER_NAN_CHECK(x) if (std::isnan(x)) { \
    std::cout << "NAN error in " << __FILE__ \
        << ":" << __LINE__ << std::endl; }
#else
#   define REMESHER_NAN_CHECK(x) /* x */
#endif

/* Caching of patches improves performance A LOT. */
#define REMESHER_PATCH_CACHING 1
#define REMESHER_PATCH_CACHE_SIZE (1 * 1024 * 1024 * 1024) // 1GB

/* Fast cache lookup requires more memory. Essential for parallelization. */
#define REMESHER_FAST_PATCH_LOOKUP 1

/* Parallelization of the Lloyd relaxation. */
#define REMESHER_PARALLELIZATION 1 // 0 = no threads, 1 = threads
#define REMESHER_RELAXATION_THREADS 8

/* The namespace definition the library uses. */
#define REMESHER_NAMESPACE_BEGIN namespace Remesher {
#define REMESHER_NAMESPACE_END }

/* Error checks. */
#if REMESHER_PARALLELIZATION \
    && REMESHER_PATCH_CACHING \
    && !REMESHER_FAST_PATCH_LOOKUP
# error "Parallelization without fast patch lookup performs very bad!"
#endif

#if REMESHER_PARALLELIZATION && REMESHER_RELAXATION_THREADS == 1
# error "Parallelization is active for a single thread!"
#endif

/* Parametrization with free boundary or not? */
#define REMESHER_FB_PARAMETRIZATION 0

/*
 * Sexy math.
 * Note that arguments for MY_MAX and MY_MIN should NOT be an
 * expensive expression, because either x or y is evaluated twice.
 */
#define MY_MAX(x, y) (x < y ? y : x)
#define MY_MIN(x, y) (x < y ? x : y)
#define MY_PI 3.14159265358979323846
#define MY_2PI (2.0 * MY_PI)
#define MY_PI2 (MY_PI / 2.0)
#define MY_RAD2DEG(x) (x*360.0/(2.0*MY_PI))
#define MY_DEG2RAD(x) (x*(2.0*MY_PI/360.0))
#define MY_FLT_ROUND(x) ((float)((int)(x + 0.5f)))
#define MY_DBL_ROUND(x) ((double)((int)(x + 0.5)))

/* Floating point epsilon comparisons. */
#define MACHINE_FLT_EPSILON 1.1920928955078125e-07
#define MACHINE_DBL_EPSILON 2.2204460492503131e-16
#define MY_FLT_EPS (MACHINE_FLT_EPSILON * 10.0f)
#define MY_DBL_EPS (MACHINE_DBL_EPSILON * 10.0)
#define EPSILON_EQ(x,v,eps) (((v - eps) <= x) && (x <= (v + eps)))
#define FLOAT_EQ(x,v) EPSILON_EQ(x,v,MY_FLT_EPS)
#define DOUBLE_EQ(x,v) EPSILON_EQ(x,v,MY_DBL_EPS)

/* Other useful defines and constants. */
#define MAX_SIZE_T ((std::size_t)-1)
#define UNUSED(arg) /* arg */

#endif /* REMESHER_DEFINES_HEADER */
