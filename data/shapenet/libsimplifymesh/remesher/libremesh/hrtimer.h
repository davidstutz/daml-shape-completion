#ifndef REMESHER_HR_TIMER_HEADER
#define REMESHER_HR_TIMER_HEADER

#ifdef WIN32
#   define NOMINMAX
#   include "windows.h"
#else
#   include <sys/time.h>
#endif

#include "defines.h"

REMESHER_NAMESPACE_BEGIN

/**
 * Cross-platform high-resolution real-time timer.
 * This implementation returns milli seconds as smallest unit.
 * The windows implementation uses functions with poor precision (~15ms).
 */
class HRTimer
{
private:
#ifdef WIN32
    std::size_t start;
#else
    struct timeval start;
#endif

public:
    HRTimer (void);

    /** Returns the milli seconds timer creation. */
    std::size_t get_elapsed (void) const;

    /** Returns the seconds since timer creation. */
    float get_elapsed_sec (void) const;
};

/* ---------------------------------------------------------------- */

inline
HRTimer::HRTimer (void)
{
#ifdef WIN32
    // FIXME: ::GetTickCount has poor precision (~10 - 16ms)
    start = ::GetTickCount();
#else
    ::gettimeofday(&start, 0);
#endif
}

inline std::size_t
HRTimer::get_elapsed (void) const
{
#ifdef WIN32
    return ::GetTickCount() - start;
#else
    struct timeval curTime;
    ::gettimeofday(&curTime, 0);
    std::size_t ret = (curTime.tv_sec - start.tv_sec) * 1000;
    std::size_t cur_ms = curTime.tv_usec / 1000;
    std::size_t start_ms = start.tv_usec / 1000;
    if (cur_ms >= start_ms)
        ret += (cur_ms - start_ms);
    else
        ret -= (start_ms - cur_ms);
    return ret;
#endif
}

inline float
HRTimer::get_elapsed_sec (void) const
{
    return (1.0f / 1000.0f) * (float)this->get_elapsed();
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_TIMER_HEADER */
