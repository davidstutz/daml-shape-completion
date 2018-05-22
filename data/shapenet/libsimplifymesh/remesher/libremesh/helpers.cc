#include <iomanip>
#include <sstream>

#include "exception.h"
#include "helpers.h"

REMESHER_NAMESPACE_BEGIN

std::string
Helpers::get_string (int value)
{
  std::stringstream ss;
  ss << value;
  return ss.str();
}

/* ---------------------------------------------------------------- */

std::string
Helpers::get_string (unsigned int value)
{
  std::stringstream ss;
  ss << value;
  return ss.str();
}

/* ---------------------------------------------------------------- */

std::string
Helpers::get_string (float value, int digits)
{
  std::stringstream ss;
  ss << std::fixed << std::setprecision(digits) << value;
  return ss.str();

}

/* ---------------------------------------------------------------- */

#if __WORDSIZE==64
std::string
Helpers::get_string (std::size_t value)
{
  std::stringstream ss;
  ss << value;
  return ss.str();

}
#endif

/* ---------------------------------------------------------------- */

std::string
Helpers::get_string (double value, int digits)
{
  std::stringstream ss;
  ss << std::fixed << std::setprecision(digits) << value;
  return ss.str();

}

/* ---------------------------------------------------------------- */

int
Helpers::get_int_from_string (std::string const& value)
{
  std::stringstream ss(value);
  int ret;
  ss >> ret;
  return ret;
}

/* ---------------------------------------------------------------- */

unsigned int
Helpers::get_uint_from_string (std::string const& value)
{
  std::stringstream ss(value);
  unsigned int ret;
  ss >> ret;
  return ret;
}

/* ---------------------------------------------------------------- */

std::size_t
Helpers::get_sizet_from_string (std::string const& value)
{
  std::stringstream ss(value);
  std::size_t ret;
  ss >> ret;
  return ret;
}

/* ---------------------------------------------------------------- */

float
Helpers::get_float_from_string (std::string const& value)
{
  std::stringstream ss(value);
  float ret;
  ss >> ret;
  return ret;
}

/* ---------------------------------------------------------------- */

double
Helpers::get_double_from_string (std::string const& value)
{
  std::stringstream ss(value);
  double ret;
  ss >> ret;
  return ret;
}

/* ---------------------------------------------------------------- */

std::string
Helpers::get_dotted_str (int value)
{
  std::stringstream ss;
  ss << value;
  return Helpers::get_dotted_str(ss.str());
}

/* ---------------------------------------------------------------- */

std::string
Helpers::get_dotted_str (unsigned int value)
{
  std::stringstream ss;
  ss << value;
  return Helpers::get_dotted_str(ss.str());
}

/* ---------------------------------------------------------------- */

#if __WORDSIZE==64
std::string
Helpers::get_dotted_str (std::size_t value)
{
  std::stringstream ss;
  ss << value;
  return Helpers::get_dotted_str(ss.str());
}
#endif

/* ---------------------------------------------------------------- */

std::string
Helpers::get_dotted_str (std::string const& str)
{
  std::string ret;

  int cnt = 0;
  for (int i = (int)str.size() - 1; i >= 0; --i)
  {
    if (cnt % 3 == 0 && cnt > 0)
      ret.insert(ret.begin(), 1, ',');
    ret.insert(ret.begin(), 1, str[i]);
    cnt += 1;
  }
  return ret;
}

/* ---------------------------------------------------------------- */

StringVector
Helpers::split_string (const std::string& str, char delim)
{
  StringVector parts;

  std::size_t last = 0;
  std::size_t cur = 0;
  for (; cur < str.size(); ++cur)
    if (str[cur] == delim)
    {
      parts.push_back(str.substr(last, cur - last));
      last = cur + 1;
    }

  if (last < str.size())
    parts.push_back(str.substr(last));

  return parts;
}

/* ---------------------------------------------------------------- */

void
Helpers::chop_string (std::string& str)
{
  while (!str.empty() && (str[str.size() - 1] == '\r'
      || str[str.size() - 1] == '\n'))
    str.resize(str.size() - 1);
}

/* ---------------------------------------------------------------- */

void
Helpers::clip_string (std::string& str)
{
  while (!str.empty() && (str[str.size() - 1] == ' '
      || str[str.size() - 1] == '\t'))
    str.erase(--str.end());

  while (!str.empty() && (str[0] == ' ' || str[0] == '\t'))
    str.erase(str.begin());
}

/* ---------------------------------------------------------------- */

void
Helpers::normalize_string (std::string& str)
{
  std::size_t iter = 0;
  bool was_whitespace = false;
  while (iter < str.size())
  {
    if (str[iter] == '\t')
      str[iter] = ' ';

    if (str[iter] == ' ')
    {
      if (was_whitespace)
      {
        str.erase(str.begin() + iter);
        was_whitespace = true;
        continue;
      }
      was_whitespace = true;
    }
    else
    {
      was_whitespace = false;
    }

    iter += 1;
  }
}

/* ---------------------------------------------------------------- */

template <>
void
Helpers::swap_endianess (float& subject)
{
  char* data = (char*)&subject;
  std::swap(data[0], data[3]);
  std::swap(data[1], data[2]);
}

/* ---------------------------------------------------------------- */

template <>
void
Helpers::swap_endianess (unsigned int& subject)
{
  char* data = (char*)&subject;
  std::swap(data[0], data[3]);
  std::swap(data[1], data[2]);
}

/* ---------------------------------------------------------------- */

template <>
void
Helpers::swap_endianess (int& subject)
{
  char* data = (char*)&subject;
  std::swap(data[0], data[3]);
  std::swap(data[1], data[2]);
}

/* ---------------------------------------------------------------- */

template <class T>
void
Helpers::swap_endianess (T& subject)
{
  throw Exception("No template secialization found");
}

REMESHER_NAMESPACE_END
