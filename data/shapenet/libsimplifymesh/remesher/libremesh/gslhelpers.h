#ifndef REMESHER_GSL_HELPERS_HEADER
#define REMESHER_GSL_HELPERS_HEADER

#include <iostream>
#include <iomanip>
#include <string>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

namespace GSL
{
  void
  pretty_print_matrix (gsl_matrix* m, std::string const& name)
  {
    std::cout << name << " = [";
    for (std::size_t i = 0; i < m->size2; ++i)
    {
      if (i != 0) std::cout << ";" << std::endl << "     ";
      for (std::size_t j = 0; j < m->size1; ++j)
      {
        if (j != 0) std::cout << ", ";
        std::cout << std::setw(6) << std::fixed << std::setprecision(2)
            << gsl_matrix_get(m, i, j);
      }
    }
    std::cout << "];" << std::endl;
  }

  void
  pretty_print_vector (gsl_vector* v, std::string const& name)
  {
    std::cout << name << " = [";
    for (std::size_t i = 0; i < v->size; ++i)
    {
      if (i != 0) std::cout << "; ";
      std::cout << std::setw(6) << std::fixed << std::setprecision(2)
            << gsl_vector_get(v, i);
    }
    std::cout << "];" << std::endl;
  }

}

#endif /* REMESHER_GSL_HELPERS_HEADER */
