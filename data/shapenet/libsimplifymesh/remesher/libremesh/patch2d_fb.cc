#include <iomanip> // REMOVE
#include <iostream>
#include <algorithm>
#include <stack>

#define GETFEM_PARA_LEVEL 0
#include "libgmm/gmm_dense_Householder.h"
#include "libgmm/gmm_matrix.h"
#include "libgmm/gmm_iter_solvers.h"
//#include "libnl/nl.h"

#include "elapsedtimer.h"
#include "exception.h"
#include "modelwriter.h"
#include "patch2d_fb.h"

REMESHER_NAMESPACE_BEGIN

void
Patch2dFB::create_patch (Patch3dPtr patch)
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
Patch2dFB::create_patch (TriangleMeshPtr mesh)
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
Patch2dFB::create_index_mappings (Patch3dPtr patch)
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

#if 0
    std::cout << std::endl << "Vertex mapping" << std::endl;
    for (std::size_t i = 0; i < this->vmap.size(); ++i)
        std::cout << i << " -> " << this->vmap[i] << std::endl;

    std::cout << std::endl << "Face mapping" << std::endl;
    for (std::size_t i = 0; i < this->fmap.size(); ++i)
        std::cout << i << " -> " << this->fmap[i] << std::endl;
#endif
}

/* ---------------------------------------------------------------- */

void
Patch2dFB::create_identity_index_mappings (void)
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
Patch2dFB::create_patch_mesh (void)
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
Patch2dFB::create_parametrization (void)
{
  /* Create vertex information for the 3D patch. */
  VertexInfoListPtr vinfo = VertexInfoList::create(this->mesh2d);

  /* Reorders vertices of the mesh in a way such that all boundary
   * vertices are moved to the end of the vertex list. Data structures
   * are updated appropriately. */
  this->reorder_boundary_vertices(vinfo);

#if CREATE_DEBUG_3D_MESHES
  std::cout << "Saving model with " << this->mesh2d->get_vertices().size()
    << " vertices and " << (this->mesh2d->get_faces().size() / 3)
    << " faces." << std::endl;
  ModelWriter::save_off_model(DEBUG_MESH_3D_PATCH_FILE, this->mesh2d);
#endif

  /* This is only to check if there is exactly one boundary. */
  {
    PatchBoundary boundary;
    this->create_patch_boundary(boundary, vinfo);
  }
  
  this->create_embedding(vinfo);

#if CREATE_DEBUG_2D_MESHES

// Platform agnostic "wait for enter": #include <iostream> #include <limits> void pause() { using namespace std; cin.clear(); for ( char c; cin.readsome(&c, 1); ) {} cin.ignore(numeric_limits<streamsize>::max(), '\n'); }

  ModelWriter::save_off_model(DEBUG_MESH_2D_PATCH_FILE, this->mesh2d);
  std::cout << "Patches saved... press a key" << std::endl;
  std::string str;
  std::getline(std::cin, str);
#endif
}

/* ---------------------------------------------------------------- */

void
Patch2dFB::reorder_boundary_vertices (VertexInfoListPtr vinfo)
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
  }

  if (border_seed == MAX_SIZE_T)
    throw Exception("There is no border (the patch is closed)");

  /* Use the border seed point to build the patch boundary. */
  std::size_t prev = border_seed;
  std::size_t current = border_seed;
  while (prev == border_seed || current != border_seed)
  {
    std::size_t adj_tri = vinfo[current].adj_faces[0];
    bool found_new = false;
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
      found_new = true;
    }
    else
    //if (!found_new)
    {
      //ModelWriter::save_off_model(DEBUG_MESH_3D_PATCH_FILE, this->mesh2d);
      //std::cout << "Cannot find new border vertex" << std::endl;
      throw Exception("Error finding patch boundary");
    }
  }

  if (b.size() != border_vamount)
  {
    //std::cout << "Mesh is not isomorphic to a disc" << std::endl;
    throw Exception("Mesh is not isomorphic to a disc");
  }
}

/* ---------------------------------------------------------------- */

void
Patch2dFB::anglesToMesh(VertexInfoListPtr vinfo, std::size_t v1, std::size_t v2)
{
    MeshVertexList & verts3d = this->mesh3d->get_vertices();
    MeshVertexList & verts2d = this->mesh2d->get_vertices();
    MeshFaceList & faces2d = this->mesh2d->get_faces();
     
    verts2d.assign(verts2d.size(), Vec3f(0.0f, 0.0f, 1.0f));

    float edge_length = (verts3d[vmap[v1]] - verts3d[vmap[v2]]).length();
    verts2d[v1] = Vec3f(0.0f, 0.0f, 0.0f);
    verts2d[v2] = Vec3f(edge_length, 0.0f, 0.0f);

//    std::cout << "Starting vertices: " << verts3d[v1] << ", " << verts3d[v2] << std::endl;
//    std::cout << "Projected: " << verts2d[v1] << ", " << verts2d[v2] << std::endl;

    std::vector<bool> faces(faces2d.size()/3, false);
    std::stack<std::pair<std::size_t, std::size_t> > edge_stack;

    edge_stack.push(std::pair<std::size_t, std::size_t>(v1, v2));

    while (!edge_stack.empty())
    {
    std::size_t vA = edge_stack.top().first;
    std::size_t vB = edge_stack.top().second;
    edge_stack.pop();
    //std::cout << "POP " << vA << ", " << vB << std::endl;

    for (std::size_t i=0; i<vinfo[vA].adj_faces.size(); i++)
    {
        for (std::size_t j=0; j<vinfo[vB].adj_faces.size(); j++)
        {
            // find faces containing this edge
            std::cout << vinfo[vA].adj_faces[i]
                << " " << vinfo[vB].adj_faces[j] << std::endl;
            if (vinfo[vA].adj_faces[i] != vinfo[vB].adj_faces[j])
                continue;
#if 0
            std::cout << "Processing face " << vinfo[vA].adj_faces[i]
                << " (" << faces2d[3*vinfo[vA].adj_faces[i] + 0]
                << ", " << faces2d[3*vinfo[vA].adj_faces[i] + 1]
                << ", " << faces2d[3*vinfo[vA].adj_faces[i] + 2]
                << ")." << std::endl;
#endif
            // is this face marked as set?
            if (faces[vinfo[vA].adj_faces[i]]) {
            //std::cout << "Face already marked." << std::endl;
                continue;
            }

            //std::cout << "Face not yet marked." << std::endl;

            // find indices of vA and vB within the face
            std::size_t vAindex = 0;
            for (vAindex=0; vAindex<3; vAindex++)
            {
                //std::cout << "vAindex: " << vAindex << " vertID: "
                //    << faces2d[3*vinfo[vA].adj_faces[i] + vAindex]
                //    << std::endl;
                if (faces2d[3*vinfo[vA].adj_faces[i] + vAindex] == vA)
                    break;
            }

            std::size_t vBindex=0;
            for (vBindex=0; vBindex<3; vBindex++)
                if (faces2d[3*vinfo[vB].adj_faces[j] + vBindex] == vB)
                    break;

            std::size_t vCindex = 0;
            for (std::size_t k=0; k<3; k++)
            {
                if ((k == vBindex) || (k == vAindex))
                    continue;
                vCindex = k;
            }

            //std::cout << "indices within face: " << vAindex
            //    << ", " << vBindex << ", " << vCindex << std::endl;

            std::size_t vC = faces2d[3*vinfo[vA].adj_faces[i] + vCindex];
            if ((verts2d[vC])[2] == 0.0f)
            {
                // already projected
                // std::cout << "Vertex " << vC << "(" << vinfo[vA].adj_faces[i] << ", " << vCindex << ") already projected." << std::endl;
                faces[vinfo[vA].adj_faces[i]] = true;
                continue;
            }

            std::size_t vFaceIndex = vinfo[vA].adj_faces[i];
            std::size_t vATmp = vA;
            std::size_t vBTmp = vB;
            if ((vAindex + 1)%3 != vBindex)
            {
                std::swap(vATmp, vBTmp);
                std::swap(vAindex, vBindex);
            }

            // compute vC in 2d
            Vec3f hVec(verts2d[vBTmp] - verts2d[vATmp]);
            double b = hVec.length()
                    * std::sin(alphas[3*vFaceIndex + vBindex])
                    / std::sin(M_PI - alphas[3*vFaceIndex + vBindex]
                    - alphas[3*vFaceIndex + vAindex]);

            hVec.norm_self();
            verts2d[vC] = verts2d[vATmp] + hVec * (float)(b * std::cos(alphas[3*vFaceIndex + vAindex]));

            float tmp = hVec[0];
            hVec[0] = -hVec[1];
            hVec[1] = tmp;

            verts2d[vC] += hVec * (float)(b * sin(alphas[3*vFaceIndex + vAindex]));

//		std::cout << "verts2d[" << vC << "] : " << verts2d[vC] << std::endl;

/*
		// compute vC in 2d
		Vec3f hVec(verts2d[vB] - verts2d[vA]);
		if (vBindex < vAindex) {
		    hVec *= -1.0f;
		}
		float b = hVec.length() * sin(alphas[3*vinfo[vB].adj_faces[j] + vBindex]) / sin(M_PI - alphas[3*vinfo[vB].adj_faces[j] + vBindex] - alphas[3*vinfo[vA].adj_faces[i] + vAindex]);

		hVec.norm_self();
		verts2d[vC] = verts2d[vA] + hVec * b * cos(alphas[3*vinfo[vA].adj_faces[i] + vAindex]);

		if (vCindex == 1) {
		    hVec *= -1.0f;
		}
		float tmp = hVec[0];
		hVec[0] = -hVec[1];
		hVec[1] = tmp;

		verts2d[vC] += hVec * b * sin(alphas[3*vinfo[vA].adj_faces[i] + vAindex]);
		std::cout << "verts2d[" << vC << "] : " << verts2d[vC] << std::endl;
*/

            faces[vFaceIndex] = true;
            edge_stack.push(std::pair<std::size_t, std::size_t>(vA, vC));
            edge_stack.push(std::pair<std::size_t, std::size_t>(vB, vC));
            // std::cout << "PUSH " << vA << ", " << vC << std::endl;
            // std::cout << "PUSH " << vB << ", " << vC << std::endl;
        }
    }
    } // while
}

/* ---------------------------------------------------------------- */
void
Patch2dFB::linearAnglesToMesh()
{
  ElapsedTimer timer;

  MeshVertexList& verts2d = this->mesh2d->get_vertices();
  MeshFaceList& faces2d = this->mesh2d->get_faces();

  std::size_t numVerts = verts2d.size();
  std::size_t numFaces = faces2d.size() / 3;

	std::size_t matrix_cols = 2*numVerts;

#if 1
	/* gmm Version */
  std::size_t matrix_rows = 2*numFaces+4;
  typedef gmm::row_matrix<gmm::rsvector<double> > GMMMatrix;
	GMMMatrix matrix_a(matrix_rows,matrix_cols);
  std::vector<double> vector_b(matrix_rows, 0.0);
	
	matrix_a(2*numFaces,2*faces2d[0]) = 1;
	matrix_a(2*numFaces+1,2*faces2d[0]+1) = 1;
	matrix_a(2*numFaces+2,2*faces2d[1]) = 1;
	matrix_a(2*numFaces+3,2*faces2d[1]+1) = 1;
	vector_b[2*numFaces+0] = 0.0f;
	vector_b[2*numFaces+1] = 0.0f;
	vector_b[2*numFaces+2] = 1.0f;
	vector_b[2*numFaces+3] = 0.0f;
	
	for (std::size_t f=0 ; f < numFaces ; ++f)
	{
		double a1 = alphas[f*3];
		double a2 = alphas[f*3+1];
		double a3 = alphas[f*3+2];
		double sina1 = std::sin(a1);
		double sina2 = std::sin(a2);
		double sina3 = std::sin(a3);

	  std::size_t v1 = faces2d[f*3];
		std::size_t v2 = faces2d[f*3 + 1];
		std::size_t v3 = faces2d[f*3 + 2];

	  double ratio = (sina3==0.0f) ? 1.0f : sina2/sina3;
		double cosine = std::cos(a1)*ratio;
		double sine = sina1*ratio;

		matrix_a(2*f, 2*v1) = cosine - 1.0f;
		matrix_a(2*f, 2*v1+1) = -sine;
		matrix_a(2*f, 2*v2) =  -cosine;
		matrix_a(2*f, 2*v2+1) = sine;
		matrix_a(2*f, 2*v3) = 1.0f;

		matrix_a(2*f+1, 2*v1) = sine;
		matrix_a(2*f+1, 2*v1+1) = cosine - 1.0f;
		matrix_a(2*f+1, 2*v2) =  -sine;
		matrix_a(2*f+1, 2*v2+1) = -cosine;
		matrix_a(2*f+1, 2*v3+1) = 1.0f;
	}
  
	std::vector<double> vector_x (matrix_cols, 1.0f);
  gmm::iteration iter(MY_FLT_EPS); //1.0e-8
  iter.set_maxiter(5*numVerts);
	gmm::least_squares_cg(matrix_a, vector_x, vector_b, iter);
	for (std::size_t v=0 ; v < numVerts ; ++v)
	{
			Vec3f vec((float)vector_x[2*v], (float)vector_x[2*v+1], 0.0f);
			verts2d[v] = vec;
	}
	
#else

	/* OpenNL Version */
	std::size_t fix1 = faces2d[0];
	std::size_t fix2 = faces2d[1];
	nlNewContext(); 
  nlSolverParameteri(NL_SOLVER, NL_CG);
  nlSolverParameteri(NL_NB_VARIABLES, matrix_cols);
  nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
  nlSolverParameteri(NL_MAX_ITERATIONS, 5*numVerts);
  nlSolverParameterd(NL_THRESHOLD, MY_FLT_EPS); //e-10
  nlBegin(NL_SYSTEM);
	nlSetVariable(2*fix1, 0.0f);
	nlSetVariable(2*fix1+1, 0.0f);
	nlSetVariable(2*fix2, 1.0f);
	nlSetVariable(2*fix2+1, 0.0f);
	nlLockVariable(2*fix1);
	nlLockVariable(2*fix1+1);
	nlLockVariable(2*fix2);
	nlLockVariable(2*fix2+1);

  nlBegin(NL_MATRIX);
	for (std::size_t f=0 ; f < numFaces ; ++f)
	{
		double a1 = alphas[f*3];
		double a2 = alphas[f*3+1];
		double a3 = alphas[f*3+2];
		double sina1 = std::sin(a1);
		double sina2 = std::sin(a2);
		double sina3 = std::sin(a3);

	  std::size_t v1 = faces2d[f*3];
		std::size_t v2 = faces2d[f*3 + 1];
		std::size_t v3 = faces2d[f*3 + 2];

	  double ratio = (sina3==0.0f) ? 1.0f : sina2/sina3;
		double cosine = std::cos(a1)*ratio;
		double sine = sina1*ratio;

	 	nlBegin(NL_ROW);
		nlCoefficient(2*v1,   cosine - 1.0);
		nlCoefficient(2*v1+1, -sine);
		nlCoefficient(2*v2,   -cosine);
		nlCoefficient(2*v2+1, sine);
		nlCoefficient(2*v3,   1.0);
		nlEnd(NL_ROW);

		nlBegin(NL_ROW);
		nlCoefficient(2*v1,   sine);
		nlCoefficient(2*v1+1, cosine - 1.0);
		nlCoefficient(2*v2,   -sine);
		nlCoefficient(2*v2+1, -cosine);
		nlCoefficient(2*v3+1, 1.0);
		nlEnd(NL_ROW);
	}
  
	nlEnd(NL_MATRIX);
	nlEnd(NL_SYSTEM);

	nlSolve();

	for (std::size_t v=0 ; v < numVerts ; ++v)
	{
			Vec3f vec (nlGetVariable(2*v), nlGetVariable(2*v+1), 0.0f);
			verts2d[v] = vec;
	}
#endif

  //std::cout << "Angles to mesh took " << timer.get_elapsed()
  //    << " ms" << std::endl;
}
/* ---------------------------------------------------------------- */
float
Patch2dFB::getAngle3D(std::size_t t, std::size_t k)
{
    MeshVertexList& verts3d = this->mesh3d->get_vertices();
    MeshFaceList& faces3d = this->mesh3d->get_faces();

#if 0
    std::size_t id1 = fmap[t] * 3 + (k + 0) % 3;
    std::size_t id2 = fmap[t] * 3 + (k + 1) % 3;
    std::size_t id3 = fmap[t] * 3 + (k + 2) % 3;
    std::cout << "Faces: " << faces3d.size()
        << ", IDs: " << id1 << "," << id2 << "," << id3 << std::endl;
#endif

    std::size_t v1 = faces3d[fmap[t] * 3 + (k + 0) % 3];
    std::size_t v2 = faces3d[fmap[t] * 3 + (k + 1) % 3];
    std::size_t v3 = faces3d[fmap[t] * 3 + (k + 2) % 3];

    // ??? Simon
    //v1 = vmap[v1];
    //v2 = vmap[v2];
    //v3 = vmap[v3];

    float angle = ((verts3d[v2] - verts3d[v1]).norm())
        .scalar((verts3d[v3] - verts3d[v1]).norm());

//    float angle =
//	((verts3d[faces3d[t * 3 + (k + 1) % 3]] -
//	  verts3d[faces3d[t * 3 + k]]).norm()).
//	scalar((verts3d[faces3d[t * 3 + (k + 2) % 3]] -
//		verts3d[faces3d[t * 3 + k]]).norm());
//    std::cout << "Computing angle " << k << " in face " << t << ": " << acos(angle) << std::endl;
//    std::cout << "Vertices: " << verts3d[faces3d[t * 3 + k]] << ", " << verts3d[faces3d[t * 3 + (k + 1) % 3]] << ", " << verts3d[faces3d[t * 3 + (k + 2) % 3]] << std::endl;
//      angle /= (verts3d[this->fmap[i]*3 + (k+1)%2] - verts3d[this->fmap[i]*3 + k]).length() * (verts3d[this->fmap[i]*3 + (k+2)%2] - verts3d[this->fmap[i]*3 + k]).length();


    angle = MY_MIN(1.0f, angle);
    angle = MY_MAX(-1.0f, angle);
    angle = std::acos(angle);
    REMESHER_NAN_CHECK(angle);

    return angle;
}

/* ---------------------------------------------------------------- */

float Patch2dFB::getAngleSum3d(VertexInfoListPtr vinfo, std::size_t vId)
{
    MeshVertexList & verts3d = this->mesh3d->get_vertices();

    // compute weights for optimal angles
    float angleSum = 0.0f;

    VertexInfo::VertexList vlist;
    vinfo->get_adjacent_vertices(vId, vlist);
    for (std::size_t i = 0; i < vlist.size(); ++i)
	vlist[i] = vmap[vlist[i]];
    vId = vmap[vId];

    Vec3f const& center(verts3d[vId]);
    for (std::size_t i = 0; i < vlist.size(); ++i)
    {
	std::size_t ip1 = (i + 1) % vlist.size();
	Vec3f const& v1(verts3d[vlist[i]]);
	Vec3f const& v2(verts3d[vlist[ip1]]);
	Vec3f e1((v1 - center).norm());
	Vec3f e2((v2 - center).norm());
	float cosangle = e1.scalar(e2);

        cosangle = MY_MIN(1.0f, cosangle);
	cosangle = MY_MAX(-1.0f, cosangle);

//	std::cout << "angle: " << std::acos(cosangle) << std::endl;
	angleSum += std::acos(cosangle);
    }
    //std::cout << "final sum: " << angleSum << std::endl;


    return angleSum;
}

/* ---------------------------------------------------------------- */

void
Patch2dFB::compute_optimalAngles(VertexInfoListPtr vinfo)
{
    MeshFaceList & faces2d = this->mesh2d->get_faces();

    std::size_t totalFaces = faces2d.size() / 3;

    this->optimalAngles.clear();
    this->optimalAngles.reserve(3 * totalFaces);

    for (std::size_t t = 0; t < totalFaces; ++t)
    {
		// compute angles in all 3d faces
		for (std::size_t k = 0; k < 3; k++)
		{
	    	float angle = this->getAngle3D(t, k);
	    	if (vinfo[faces2d[t * 3 + k]].vclass != VERTEX_CLASS_BORDER) {
					float angleSum =  this->getAngleSum3d(vinfo, faces2d[t * 3 + k]);
					angle *= 2.0f * (float)(MY_PI / angleSum);
	    	}
	    	this->optimalAngles.push_back(angle);
		}
    }
}

/* ---------------------------------------------------------------- */

void
Patch2dFB::create_embedding (VertexInfoListPtr vinfo)
{
  MeshVertexList& verts2d = this->mesh2d->get_vertices();
  MeshFaceList& faces2d = this->mesh2d->get_faces();
  
  ElapsedTimer timer;

  //std::cout << "Computing optimal angles .... " << std::endl;
  this->compute_optimalAngles(vinfo);
	
  std::size_t totalFaces = faces2d.size() / 3;

	std::size_t internalVerts = 0;
	for (std::size_t t = 0; t < verts2d.size(); ++t)
	{
	    if (vinfo[t].vclass != VERTEX_CLASS_SIMPLE)
				break;
			internalVerts++;
	}
  //std::cout << optimalAngles << std::endl;
	std::size_t matrix_rows = internalVerts*2 + totalFaces;
	std::size_t matrix_cols = this->optimalAngles.size();
  //std::cout << "Matrix A -> rows: " << matrix_rows << " - cols: " << matrix_cols << " ... " << std::endl;

  typedef gmm::row_matrix<gmm::rsvector<double> > GMMMatrix;
  GMMMatrix matrix_a(matrix_rows, matrix_cols);
  std::vector<double> vector_b(matrix_rows);
	
  //std::cout << "Adding constrains ... " << std::endl;
	/* add Vertex and Wheel consistency constrain */
	for (std::size_t v = 0 ; v < internalVerts ; ++v)
	{
		double sum_b_u = 0.0f;
		double sum_b_l = 0.0f;
		std::vector<std::size_t> afcs = vinfo[v].adj_faces;
		
		for (std::size_t f = 0 ; f < afcs.size() ; ++f)
		{
			for (std::size_t fv = 0 ; fv < 3 ; ++fv)
			{
				if(faces2d[afcs[f]*3 +fv] == v)
				{
					double cangle = optimalAngles[afcs[f]*3 + fv];
					matrix_a(v, afcs[f]*3 + fv) = 1;
					sum_b_u += cangle;

					double langle = optimalAngles[afcs[f]*3 + ((fv+1) % 3)];
					matrix_a(internalVerts + v, afcs[f]*3 + ((fv+1) % 3)) = 1/tan(langle); 
					double rangle = optimalAngles[afcs[f]*3 + ((fv+2) % 3)];
					matrix_a(internalVerts + v, afcs[f]*3 + ((fv+2) % 3)) = -1/tan(rangle);
					sum_b_l += log(sin(rangle)) - log(sin(langle));
				}
			}
		}
		vector_b[v] = 2.0f * M_PI - sum_b_u;
		vector_b[internalVerts + v] = sum_b_l;
	}


  /* add Triangle consistency constrain */
	for (std::size_t f = 0 ; f < totalFaces ; ++f)
	{
		double sum_b = 0.0f;
		for (std::size_t v = 0 ; v < 3 ; ++v)
		{
			sum_b += optimalAngles[f*3 + v];
			matrix_a(internalVerts*2 + f, f*3 + v) = 1;
		}
		vector_b[internalVerts*2 + f] = M_PI - sum_b;
	}
	
  //std::cout << " A : " << matrix_a << std::endl;
  //std::cout << "Computing matrix c ... " << std::endl;
	/* compute matrix c */
  GMMMatrix matrix_d(matrix_cols, matrix_cols);
  GMMMatrix matrix_c(matrix_rows, matrix_cols);
  for (std::size_t i=0 ; i < matrix_cols ; ++i)
	{
		matrix_d(i,i) = optimalAngles[i];
	}
	gmm::mult(matrix_a, matrix_d, matrix_c);
	

	/*compute C*C^t */
  //std::cout << "Computing matrix C^T ... " << std::endl;
	gmm::transposed_row_ref<GMMMatrix*> matrix_ct = gmm::transposed(matrix_c);

	GMMMatrix matrix_cct(matrix_rows, matrix_rows);

  //std::cout<< "C : " << matrix_c << std::endl << "Ct : " << matrix_ct << std::endl;
	gmm::mult (matrix_c, matrix_ct, matrix_cct);

  //std::cout << "Building matrix took " << timer.get_elapsed()
  //    << " ms" << std::endl;
	
  std::vector<double> vector_x(matrix_rows);
  gmm::ilut_precond<GMMMatrix> matrix_precond(matrix_cct, 5, 1.0e-4);
  //gmm::ilu_precond<GMMMatrix> matrix_precond(matrix_cct);
  gmm::iteration iter(1.0e-8);
  iter.set_maxiter(1000);
  //std::cout << "Solving linear system (C C^t) x = b... " << std::endl;
  gmm::gmres(matrix_cct, vector_x, vector_b, matrix_precond, 30, iter);
  //gmm::least_squares_cg(matrix_cct, vector_x, vector_b, iter);
  //std::cout << "Took " << timer.get_elapsed()
  //    << " ms" << std::endl;

	/* compute estimation error */
  //std::cout << "Compute result ... " << std::endl;
  std::vector<double> vector_eerror(matrix_cols);
  std::vector<double> vector_r(matrix_cols);
  GMMMatrix matrix_e(matrix_cols, matrix_rows);
  gmm::mult(matrix_ct, vector_x, vector_r);

  gmm::mult(matrix_d, vector_r, vector_eerror);

  /* compute solution */
  std::vector<double> vector_solution(matrix_cols);
  alphas.clear();
  for (std::size_t i=0 ; i < matrix_cols ; ++i)
	{
		vector_solution[i] = optimalAngles[i] + vector_eerror[i];
		alphas.push_back(vector_solution[i]);
		//alphas.push_back(optimalAngles[i]);
	}

  //std::cout << "Mesh flattening " << timer.get_elapsed()
  //    << " ms" << std::endl;

  // std::cout << "Final angles: " << alphas << std::endl;

  this->linearAnglesToMesh();
  //this->anglesToMesh(vinfo,faces2d[0],faces2d[1]);
}

/* ---------------------------------------------------------------- */

std::size_t
Patch2dFB::lookup_intern (IndexMapping const& map, std::size_t i,
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
Patch2dFB::get_memory_usage (void) const
{
  std::size_t s_vmap = this->vmap.capacity() * sizeof(std::size_t);
  std::size_t s_fmap = this->fmap.capacity() * sizeof(std::size_t);
  std::size_t s_mesh2d = this->mesh2d->get_memory_usage();

  return s_vmap + s_fmap + s_mesh2d;
}

REMESHER_NAMESPACE_END
