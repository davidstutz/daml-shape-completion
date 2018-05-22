#include "exception.h"
#include "averageplane.h"

#if AVG_PLANE_USE_GMM
#   define GETFEM_PARA_LEVEL 0
#   include "libgmm/gmm_matrix.h"
#   include "libgmm/gmm_dense_qr.h"
#else
#   include <gsl/gsl_matrix.h>
#   include <gsl/gsl_vector.h>
#   include <gsl/gsl_linalg.h>
#endif

REMESHER_NAMESPACE_BEGIN

Plane3f
AveragePlane::get_average (std::vector<Vec3f> const& points)
{
  if (points.size() < 3)
    throw Exception("AveragePlane: Expected at least 3 points!");

  /* Calculate centroid. */
  Vec3f c(0.0f, 0.0f, 0.0f);
  for (std::size_t i = 0; i < points.size(); ++i)
    c += points[i];
  c /= (float)points.size();

#if AVG_PLANE_USE_GMM

  /*
   * This version always creates a 3x3 Matrix and solves for eigenvalues.
   */

  float a00 = 0.0f;
  float a11 = 0.0f;
  float a22 = 0.0f;
  float a01 = 0.0f;
  float a02 = 0.0f;
  float a12 = 0.0f;

  for (std::size_t i = 0; i < points.size(); ++i)
  {
    a00 += (points[i][0] - c[0]) * (points[i][0] - c[0]);
    a11 += (points[i][1] - c[1]) * (points[i][1] - c[1]);
    a22 += (points[i][2] - c[2]) * (points[i][2] - c[2]);
    a01 += (points[i][0] - c[0]) * (points[i][1] - c[1]);
    a02 += (points[i][0] - c[0]) * (points[i][2] - c[2]);
    a12 += (points[i][1] - c[1]) * (points[i][2] - c[2]);
  }

  /* Setup matrix. */
  gmm::dense_matrix<float> matrix_m(3, 3);
  matrix_m(0, 0) = a00;  matrix_m(0, 1) = a01;  matrix_m(0, 2) = a02;
  matrix_m(1, 0) = a01;  matrix_m(1, 1) = a11;  matrix_m(1, 2) = a12;
  matrix_m(2, 0) = a02;  matrix_m(2, 1) = a12;  matrix_m(2, 2) = a22;

  /* Solve for eigenvalues and eigenvectors. */
  std::vector<float> vec_ev(3);
  gmm::dense_matrix<float> matrix_ev(3, 3);
  gmm::symmetric_qr_algorithm(matrix_m, vec_ev, matrix_ev, MY_FLT_EPS);

  int min_elem = 0;
  for (int i = 1; i < 3; ++i)
    if (std::abs(vec_ev[i]) < std::abs(vec_ev[min_elem]))
        min_elem = i;

  Vec3f normal(matrix_ev(0, min_elem),
      matrix_ev(1, min_elem),
      matrix_ev(2, min_elem));

#else /* Use GSL */

  /*
   * This version always creates a Nx3 Matrix and performes SVD.
   */

  /* Create Matrix M */
  gsl_matrix* matrix_m = gsl_matrix_alloc(points.size(), 3);
  for (std::size_t i = 0; i < points.size(); ++i)
  {
    Vec3f p = points[i] - c;
    //std::cout << p[0] << " " << p[1] << " " << p[2] << std::endl;
    gsl_matrix_set(matrix_m, i, 0, p[0]);
    gsl_matrix_set(matrix_m, i, 1, p[1]);
    gsl_matrix_set(matrix_m, i, 2, p[2]);
  }

  /* Singular Value Decomposition of M = U S V^T. */
  gsl_matrix* matrix_v = gsl_matrix_alloc(3, 3);
  gsl_vector* vector_s = gsl_vector_alloc(3);
  gsl_vector* gsl_work = gsl_vector_alloc(3);

  /* Singular value decomposition, M = USV^T. M is replaced by U. */
  gsl_linalg_SV_decomp(matrix_m, matrix_v, vector_s, gsl_work);

  //gsl_matrix_fprintf(stdout, matrix_v, "%f");
  //std::cout << std::endl;
  //gsl_vector_fprintf(stdout, vector_s, "%f");

  gsl_matrix_free(matrix_m);
  gsl_vector_free(gsl_work);

  /* Retrieve plane normal. */
  Vec3f normal;

  {
    /* Find singular vector for smallest singular value. */
    std::size_t sidx = 0;
    double sval = 0.0;
    for (std::size_t i = 0; i < 3; ++i)
    {
      double si = gsl_vector_get(vector_s, i);
      if (i == 0 || sval > si)
      {
        sidx = i;
        sval = si;
      }
    }

    normal[0] = (float)gsl_matrix_get(matrix_v, 0, sidx);
    normal[1] = (float)gsl_matrix_get(matrix_v, 1, sidx);
    normal[2] = (float)gsl_matrix_get(matrix_v, 2, sidx);
  }

  gsl_matrix_free(matrix_v);
  gsl_vector_free(vector_s);

#endif

  /* Set plane from centroid and normal. */
  return Plane3f(normal, c);
}

REMESHER_NAMESPACE_END
