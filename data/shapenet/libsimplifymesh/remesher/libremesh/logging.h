#ifndef REMESHER_LOGGING_HEADER
#define REMESHER_LOGGING_HEADER

#include <sstream>

#include "defines.h"

REMESHER_NAMESPACE_BEGIN

/*
 * A pretty simple logging facility. Logging is done with a simple call:
 *
 *   LOG(LOG_ERROR) << "Some number: " << 12;
 *
 * A newline is automatically inserted. Log instructions above the
 * REPORTING_LEVEL are not printed and instructions are completely
 * erased at compile time if the specified level is a compile time
 * constant, thus make sure you don't call functions with side effects
 * in a logging statement.
 */

enum LogLevel
{
  LOG_ERROR,
  LOG_WARNING,
  LOG_INFO,
  LOG_DEBUG,
  LOG_DEBUG1,
  LOG_DEBUG2,
  LOG_DEBUG3,
  LOG_DEBUG4
};

/* ---------------------------------------------------------------- */

#define REPORTING_LEVEL LOG_DEBUG2
#define LOG(level) if (level <= REPORTING_LEVEL) Log(level).get()

/* ---------------------------------------------------------------- */

class Log
{
  protected:
    LogLevel level;

  protected:
    char const* level_str (void);

  public:
    std::stringstream out;

  public:
    Log (LogLevel level);
    virtual ~Log (void);
    std::stringstream& get (void);
};

/* ---------------------------------------------------------------- */

inline
Log::Log (LogLevel level)
  : level(level)
{
   this->out << this->level_str();
}

inline
Log::~Log (void)
{
  std::cout << this->out.str() << std::endl << std::flush;
}

inline std::stringstream&
Log::get (void)
{
  return this->out;
}

inline char const*
Log::level_str (void)
{
  switch (this->level)
  {
    case LOG_ERROR:   return "Error:   ";
    case LOG_WARNING: return "Warning: ";
    case LOG_INFO:    return "Info:    ";
    case LOG_DEBUG:   return "Debug:   ";
    case LOG_DEBUG1:  return "Debug1:  ";
    case LOG_DEBUG2:  return "Debug2:  ";
    case LOG_DEBUG3:  return "Debug3:  ";
    case LOG_DEBUG4:  return "Debug4:  ";
    default:          return "Message: ";
  }
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_LOGGING_HEADER */
