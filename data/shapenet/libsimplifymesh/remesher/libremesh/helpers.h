#ifndef REMESHER_HELPERS_HEADER
#define REMESHER_HELPERS_HEADER

#include <string>
#include <vector>

#include "defines.h"

REMESHER_NAMESPACE_BEGIN

typedef std::vector<std::string> StringVector;

class Helpers
{
  public:
    static std::string get_string (int value);
    static std::string get_string (unsigned int value);
    static std::string get_string (float value, int digits);
    static std::string get_string (double value, int digits);
#if __WORDSIZE==64
    static std::string get_string (std::size_t value);
#endif

    static unsigned int get_uint_from_string (std::string const& value);
    static std::size_t get_sizet_from_string (std::string const& value);
    static int get_int_from_string (std::string const& value);
    static float get_float_from_string (std::string const& value);
    static double get_double_from_string (std::string const& value);

    static std::string get_dotted_str (int value);
    static std::string get_dotted_str (unsigned int value);
#if __WORDSIZE==64
    static std::string get_dotted_str (std::size_t value);
#endif
    static std::string get_dotted_str (std::string const& str);

    static StringVector split_string (std::string const& str, char delim);
    static void chop_string (std::string& str);
    static void clip_string (std::string& str);
    static void normalize_string (std::string& str);

    template <class T>
    static void swap_endianess (T& subject);
};

REMESHER_NAMESPACE_END

#endif /* REMESHER_HELPERS_HEADER */
