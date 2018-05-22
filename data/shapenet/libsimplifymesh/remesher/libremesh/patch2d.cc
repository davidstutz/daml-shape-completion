#include <iomanip> // REMOVE
#include <iostream>
#include <algorithm>

#define GETFEM_PARA_LEVEL 0
#include "libgmm/gmm_dense_Householder.h"
#include "libgmm/gmm_matrix.h"
#include "libgmm/gmm_iter_solvers.h"

#include "elapsedtimer.h"
#include "exception.h"
#include "modelwriter.h"
#include "patch2d.h"

REMESHER_NAMESPACE_BEGIN

void
Patch2d::create_patch (Patch3dPtr patch)
{
  /* Initialize members. */
  this->mesh3d = patch->get_mesh();
  this->mesh2d = TriangleMesh::create();
  this->create_index_mappings(patch);

  /* Create the patch mesh and parametrize. */
  this->create_patch_mesh();
  this->create_parametrization();

  /* Release reference to the reference model. */
  this->mesh3d.reset();
}

/* ---------------------------------------------------------------- */

void
Patch2d::create_patch (TriangleMeshPtr mesh)
{
  /* Initialize members. */
  this->mesh3d = mesh;
  this->mesh2d = TriangleMesh::create();
  this->create_identity_index_mappings();

  /* Create the patch mesh and parametrize. */
  this->create_patch_mesh();
  this->create_parametrization();

  /* Release reference to the reference model. */
  this->mesh3d.reset();
}

/* ---------------------------------------------------------------- */

void
Patch2d::create_index_mappings (Patch3dPtr patch)
{
  this->vmap.clear();
  this->fmap.clear();

  MeshFaceList const& faces = patch->get_mesh()->get_faces();

  /* Create mappings from new vertex and face indices to indices
   * on the reference mesh. Since indices for the reference mesh
   * are sorted, reverse lookup takes logarithmic time. */
  Patch3d::MeshFaceSet const& fset = patch->get_face_set();

  typedef std::set<std::size_t> VertexIndexSet;
  VertexIndexSet vset;
  for (Patch3d::MeshFaceSet::const_iterator iter = fset.begin();
      iter != fset.end(); ++iter)
  {
    this->fmap.push_back(*iter);
    for (int v = 0; v < 3; ++v)
      vset.insert(faces[*iter * 3 + v]);
  }

  for (VertexIndexSet::iterator iter = vset.begin();
      iter != vset.end(); ++iter)
    this->vmap.push_back(*iter);
}

/* ---------------------------------------------------------------- */

void
Patch2d::create_identity_index_mappings (void)
{
  this->vmap.clear();
  this->fmap.clear();

  std::size_t vsize = this->mesh3d->get_vertices().size();
  std::size_t fsize = this->mesh3d->get_faces().size() / 3;

  this->vmap.resize(vsize);
  for (std::size_t i = 0; i < vsize; ++i)
    this->vmap[i] = i;

  this->fmap.resize(fsize);
  for (std::size_t i = 0; i < fsize; ++i)
    this->fmap[i] = i;
}

/* ---------------------------------------------------------------- */

void
Patch2d::create_patch_mesh (void)
{
  MeshFaceList const& faces = this->mesh3d->get_faces();
  MeshVertexList const& verts = this->mesh3d->get_vertices();

  /* Create a new triangle mesh from the set of vertices. */
  MeshVertexList& verts2d = this->mesh2d->get_vertices();
  for (std::size_t i = 0; i < this->vmap.size(); ++i)
    verts2d.push_back(verts[this->vmap[i]]);

  MeshFaceList& faces2d = this->mesh2d->get_faces();
  for (std::size_t i = 0; i < this->fmap.size(); ++i)
  {
    std::size_t v1idx = this->lookup_vertex(faces[this->fmap[i] * 3 + 0]);
    std::size_t v2idx = this->lookup_vertex(faces[this->fmap[i] * 3 + 1]);
    std::size_t v3idx = this->lookup_vertex(faces[this->fmap[i] * 3 + 2]);

    faces2d.push_back((MeshVIndex)v1idx);
    faces2d.push_back((MeshVIndex)v2idx);
    faces2d.push_back((MeshVIndex)v3idx);
  }
}

/* ---------------------------------------------------------------- */

void
Patch2d::create_parametrization (void)
{
  /* Create vertex information for the 3D patch. */
  VertexInfoListPtr vinfo = VertexInfoList::create(this->mesh2d);

  /* Reorders vertices of the mesh in a way such that all boundary
   * vertices are moved to the end of the vertex list. Data structures
   * are updated appropriately. */
  this->reorder_boundary_vertices(vinfo);

#if CREATE_DEBUG_3D_MESHES
  ModelWriter::save_off_model(DEBUG_MESH_3D_PATCH_FILE, this->mesh2d);
#endif

  PatchBoundary boundary;
  this->create_patch_boundary(boundary, vinfo);

#if 0
  if (this->mesh2d->get_vertices().size() - boundary.size() == 2845)
    ModelWriter::save_off_model(DEBUG_MESH_3D_PATCH_FILE, this->mesh2d);
#endif

  this->create_embedding(boundary, vinfo);

#if CREATE_DEBUG_2D_MESHES

// Platform agnostic "wait for enter": #include <iostream> #include <limits> void pause() { using namespace std; cin.clear(); for ( char c; cin.readsome(&c, 1); ) {} cin.ignore(numeric_limits<streamsize>::max(), '\n'); }

  ModelWriter::save_off_model(DEBUG_MESH_2D_PATCH_FILE, this->mesh2d);
  std::cout << "2D patch saved... press a key" << std::endl;
  std::string str;
  std::cin >> str;
#endif
}

/* ---------------------------------------------------------------- */

void
Patch2d::reorder_boundary_vertices (VertexInfoListPtr vinfo)
{
  MeshVertexList& verts = this->mesh2d->get_vertices();
  MeshFaceList& faces = this->mesh2d->get_faces();
  std::size_t i = 0;
  std::size_t j = verts.size() - 1;

  while (i < j)
  {
    /* Find a border vertex. */
    while (i < j && vinfo[i].vclass != VERTEX_CLASS_BORDER)
      i += 1;

    /* Find a non-border vertex. */
    while (i < j && vinfo[j].vclass == VERTEX_CLASS_BORDER)
      j -= 1;

    if (i >= j)
      return;

    /* Swap vertices and fix face indices. */
    std::swap(verts[i], verts[j]);
    std::swap(vinfo[i], vinfo[j]);
    std::swap(vmap[i], vmap[j]);

    /* Walk over all faces and fix the vertex reference. */
    /* FIXME Inefficient. */
    for (std::size_t k = 0; k < faces.size(); ++k)
    {
      if (faces[k] == i)
        faces[k] = (MeshVIndex)j;
      else if (faces[k] == j)
        faces[k] = (MeshVIndex)i;
    }
  }
}

/* ---------------------------------------------------------------- */

void
Patch2d::create_patch_boundary (PatchBoundary& b, VertexInfoListPtr vinfo)
{
  MeshFaceList const& faces = this->mesh2d->get_faces();

  /* Search a border vertex and use it as seed. Count all border vertices
   * to check if the mesh is topological equivalent to a disc. */
  std::size_t border_seed = MAX_SIZE_T;
  std::size_t border_vamount = 0;
  for (std::size_t i = 0; i < vinfo->size(); ++i)
  {
    if (vinfo[i].vclass == VERTEX_CLASS_BORDER)
    {
      border_vamount += 1;
      border_seed = i;
    }
    else if (vinfo[i].vclass == VERTEX_CLASS_COMPLEX)
      throw Exception("Complex vertex in 3D patch");
  }

  if (border_seed == MAX_SIZE_T)
    throw Exception("There is no border (the patch is closed)");

  /* Use the border seed point to build the patch boundary. */
  std::size_t prev = border_seed;
  std::size_t current = border_seed;
  while (prev == border_seed || current != border_seed)
  {
    std::size_t adj_tri = vinfo[current].adj_faces[0];
    std::size_t adj_border_vertex = MAX_SIZE_T;
    for (int i = 0; i < 3; ++i)
      if (faces[adj_tri * 3 + i] == current)
        adj_border_vertex = faces[adj_tri * 3 + (i + 1) % 3];

    if (adj_border_vertex != MAX_SIZE_T
        && vinfo[adj_border_vertex].vclass == VERTEX_CLASS_BORDER
        && adj_border_vertex != prev && adj_border_vertex != current)
    {
      prev = current;
      current = adj_border_vertex;
      b.push_back(adj_border_vertex);
    }
    else
    {
      //ModelWriter::save_off_model(DEBUG_MESH_3D_PATCH_FILE, this->mesh2d);
      //std::cout << "Cannot find new border vertex" << std::endl;
      //std::string mystr;
      //std::getline(std::cin, mystr);
      throw Exception("Error finding patch boundary (SHOULD NOT HAPPEN)");
    }
  }

  if (b.size() != border_vamount)
  {
    //std::cout << "Mesh is not isomorphic to a disc" << std::endl;
    throw Exception("Patch has multiple boundaries (not disc-like)");
  }

  #if 0
  std::cout << "Boundary: ";
  for (std::size_t i = 0; i < b.size(); ++i)
    std::cout << b[i] << " ";
  std::cout << std::endl;
  #endif
}

/* ---------------------------------------------------------------- */

void
Patch2d::create_embedding (PatchBoundary& b, VertexInfoListPtr vinfo)
{
  MeshVertexList& verts2d = this->mesh2d->get_vertices();
  MeshVertexList& verts3d = this->mesh3d->get_vertices();

  /*
   * The boundary is embedded to 2D space now.
   */
  {
    float total_length = 0.0f;
    for (std::size_t i = 0; i < b.size(); ++i)
    {
      std::size_t ip1 = (i + 1) % b.size();
      total_length += (verts2d[b[i]] - verts2d[b[ip1]]).length();
    }

    float len_so_far = 0.0f;
    for (std::size_t i = 0; i < b.size(); ++i)
    {
      std::size_t ip1 = (i + 1) % b.size();
      float angle = len_so_far * (float)MY_2PI / total_length;
      len_so_far += (verts2d[b[i]] - verts2d[b[ip1]]).length();
      /* Modify vertex to have 2D coordinates. */
      verts2d[b[i]] = Vec3f(std::cos(angle), std::sin(angle), 0.0f);
    }
  }

  /*
   * Create matrix A for the linear system Ax = b. Also create the
   * right hand side "b" two times for x and y coordinates.
   */
  //ElapsedTimer timer;
  std::size_t matrix_size = verts2d.size() - b.size();

  //std::cout << "Solving linear system: " << matrix_size
  //    << "x" << matrix_size << " ... " << std::flush;
  typedef gmm::row_matrix<gmm::rsvector<double> > GMMMatrix;
  GMMMatrix matrix_a(matrix_size, matrix_size);
  std::vector<double> vector_bx(matrix_size);
  std::vector<double> vector_by(matrix_size);

  for (std::size_t i = 0; i < matrix_size; ++i)
  {
    /* Retrieve the adjacent vertices. */
    VertexInfo::VertexList vlist;
    vinfo->get_adjacent_vertices(i, vlist);
    std::size_t avs = vlist.size();

    /* Calculate the weights for the adjacent vertices. */
    std::vector<float> weights;
    weights.resize(avs);
    float total_weight = 0.0f;
    float border_weight = 0.0f;
    float inner_weight = 0.0f;

    for (std::size_t j = 0; j < weights.size(); ++j)
    {
      Vec3f const& v0 = verts3d[this->vmap[i]];
      Vec3f const& v1 = verts3d[this->vmap[vlist[(j + avs - 1) % avs]]];
      Vec3f const& v2 = verts3d[this->vmap[vlist[j]]];
      Vec3f const& v3 = verts3d[this->vmap[vlist[(j + 1) % avs]]];

      float d2len = (v2 - v0).length();
      Vec3f d1 = (v1 - v0).norm();
      Vec3f d2 = (v2 - v0) / d2len;
      Vec3f d3 = (v3 - v0).norm();
      float scalar1 = d1.scalar(d2);
      float scalar2 = d2.scalar(d3);
      scalar1 = MY_MIN(1.0f, scalar1);
      scalar1 = MY_MAX(-1.0f, scalar1);
      scalar2 = MY_MIN(1.0f, scalar2);
      scalar2 = MY_MAX(-1.0f, scalar2);
      float weight = (std::tan(std::acos(scalar1) / 2.0f)
          + std::tan(std::acos(scalar2) / 2.0f)) / d2len;

#if REMESHER_NAN_CHECKS
      if (std::isnan(weight))
      {
        std::cout << "NAN weight detected!" << std::endl;
        std::cout << std::fixed << std::setprecision(9);
        std::cout << "Scalar1: " << scalar1 << std::endl;
        std::cout << "Scalar2: " << scalar2 << std::endl;
        std::cout << "Acos1: " << std::acos(scalar1) << std::endl;
        std::cout << "Acos2: " << std::acos(scalar2) << std::endl;
        std::cout << "Tan1: " << std::tan(std::acos(scalar1)) << std::endl;
        std::cout << "Tan2: " << std::tan(std::acos(scalar2)) << std::endl;
        std::cout << "D2len: " << d2len << std::endl;
      }
#endif

      weights[j] = weight;
      total_weight += weight;
      if (vlist[j] >= matrix_size)
        border_weight += weight;
      else
        inner_weight += weight;
    }

    float border_weight_x = 0.0f;
    float border_weight_y = 0.0f;

    /* Set matrix weights for inner vertices, calc border
     * weight for boundary vertices. */
    for (std::size_t j = 0; j < weights.size(); ++j)
    {
      if (vlist[j] < matrix_size)
        matrix_a(i, vlist[j]) = -weights[j] / total_weight;
      else
      {
        border_weight_x += verts2d[vlist[j]][0] * weights[j] / total_weight;
        border_weight_y += verts2d[vlist[j]][1] * weights[j] / total_weight;
      }
    }

    /* Set diagonal element and right-hand side vector elements. */
    matrix_a(i, i) = 1.0;
    vector_bx[i] = border_weight_x;
    vector_by[i] = border_weight_y;
  }

  //std::cout << "Building matrix took " << timer.get_elapsed()
  //    << " ms" << std::endl;

  /*
   * Solve the linear system. The linear system needs
   * to be solved separately for the x and y coordinates.
   */

  /*
   * Since the linear system to be solved is sparse, we use
   * special data structures and solvers for sparse matrices.
   * GMRES is a solver for non-symmetric, positive definit matrices.
   * The preconditioner efficently reduces calculation time for
   * huge matrices.
   */
  std::vector<double> vector_xx(matrix_size);
  std::vector<double> vector_xy(matrix_size);
  //gmm::identity_matrix matrix_precond;
  gmm::ilut_precond<GMMMatrix> matrix_precond(matrix_a, 5, 1.0e-4);
  gmm::iteration iter(1.0e-8);
  iter.set_maxiter(1000);

  gmm::gmres(matrix_a, vector_xx, vector_bx, matrix_precond, 30, iter);
  gmm::gmres(matrix_a, vector_xy, vector_by, matrix_precond, 30, iter);
  //gmm::bicgstab(matrix_a, vector_xx, vector_bx, matrix_precond, iter);
  //gmm::bicgstab(matrix_a, vector_xy, vector_by, matrix_precond, iter);

  if (iter.converged())
  {
    //std::cout << "[" << matrix_size << "x" << matrix_size
    //    << ":" << iter.get_iteration() << "] " << std::flush;
    //std::cout << iter.get_iteration() << " iterations" << std::endl;
  }
  else
  {
    //std::cout << "ERROR: Solving LSE did not converge!" << std::endl;
    //ModelWriter::save_off_model(DEBUG_MESH_2D_PATCH_FILE, this->mesh2d);
    //::exit(1); // Sanity exit
    throw Exception("Linear system solver did not converge");
  }

  /* Store the coordinates back to the mesh. */
  for (std::size_t i = 0; i < matrix_size; ++i)
  {
    verts2d[i][0] = (float)vector_xx[i];
    verts2d[i][1] = (float)vector_xy[i];
    verts2d[i][2] = 0.0f;
  }

  //std::cout << "Parametrization took " << timer.get_elapsed()
  //    << " ms" << std::endl;
}

/* ---------------------------------------------------------------- */

std::size_t
Patch2d::lookup_intern (IndexMapping const& map, std::size_t i,
    std::size_t a, std::size_t b) const
{
  if (a > b)
    return MAX_SIZE_T;

  if (a == b)
  {
    if (map[a] == i)
      return a;
    return MAX_SIZE_T;
  }

  std::size_t mid = (a + b) / 2;
  if (map[mid] == i)
    return mid;

  if (map[mid] > i && mid > a)
    return this->lookup_intern(map, i, a, mid - 1);

  if (map[mid] < i && mid < b)
    return this->lookup_intern(map, i, mid + 1, b);

  return MAX_SIZE_T;
}

/* ---------------------------------------------------------------- */

std::size_t
Patch2d::get_memory_usage (void) const
{
  std::size_t s_vmap = this->vmap.capacity() * sizeof(std::size_t);
  std::size_t s_fmap = this->fmap.capacity() * sizeof(std::size_t);
  std::size_t s_mesh2d = this->mesh2d->get_memory_usage();

  return s_vmap + s_fmap + s_mesh2d;
}

REMESHER_NAMESPACE_END
