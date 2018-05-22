#ifndef REMESHER_ELAPSED_TIMER_HEADER
#define REMESHER_ELAPSED_TIMER_HEADER

#include <ctime>

#include "defines.h"

REMESHER_NAMESPACE_BEGIN

/**
 * Very simple timer class to take execution times. The reported
 * float values are in seconds, the integer values are in milli seconds.
 * The functions that provide milli seconds should be preferred.
 *
 * This class should not be used for timings that rely on the
 * actual real world time but rather for computational timings.
 * The timings here are pure processing times which may vary from
 * real time if the application is scheduled or sleeps.
 *
 * FIXME: What is the impact of I/O on the timings?
 */
class ElapsedTimer
{
private:
    std::size_t start;

public:
    ElapsedTimer (void);
    void reset (void);

    static float now_sec (void);
    static std::size_t now (void);

    std::size_t get_elapsed (void) const;
    float get_elapsed_sec (void) const;
};

/* ---------------------------------------------------------------- */

inline
ElapsedTimer::ElapsedTimer (void)
{
    this->reset();
}

inline void
ElapsedTimer::reset (void)
{
    this->start = ElapsedTimer::now();
}

inline float
ElapsedTimer::now_sec (void)
{
    return (float)std::clock() / (float)CLOCKS_PER_SEC;
}

inline std::size_t
ElapsedTimer::now (void)
{
    return ((std::size_t)(std::clock()) * 1000) / (std::size_t)CLOCKS_PER_SEC;
}

inline float
ElapsedTimer::get_elapsed_sec (void) const
{
    return (1.0f / 1000.0f) * (float)this->get_elapsed();
}

inline std::size_t
ElapsedTimer::get_elapsed (void) const
{
    return ElapsedTimer::now() - start;
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_ELAPSED_TIMER_HEADER */
