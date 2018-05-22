#ifndef REMESHER_EXCEPTION_HEADER
#define REMESHER_EXCEPTION_HEADER

#include <string>
#include <stdexcept>

#include "defines.h"

REMESHER_NAMESPACE_BEGIN

/* Universal simple exception class. */
class Exception : public std::exception, public std::string
{
  public:
    Exception (void) throw()
    { }

    Exception (std::string const& msg) throw() : std::string(msg)
    { }

    virtual ~Exception (void) throw()
    { }

    virtual const char* what (void) const throw()
    { return this->c_str(); }
};

REMESHER_NAMESPACE_END

#endif /* REMESHER_EXCEPTION_HEADER */
