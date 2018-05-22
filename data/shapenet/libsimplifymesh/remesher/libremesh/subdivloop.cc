#include <iostream>

#include "exception.h"
#include "subdivloop.h"

REMESHER_NAMESPACE_BEGIN

void
SubdivLoop::subdiv_impl (void)
{
  MeshFaceList& faces = this->mesh->get_faces();
  MeshVertexList& verts = this->mesh->get_vertices();

  /*
   * First, all vertice are classified into SMOOTH, DART, CREASE and CORNER.
   * Each CREASE is further classified into REGULAR and NON-REGULAR crease.
   */
  this->classify_vertices();

  /* Every edge in the mesh is subdivided and a new vertex in inserted.
   * The indices of the vertices are stored in this temp list. */
  NewVertexIndices newvidx;
  newvidx.resize(faces.size(), MAX_SIZE_T);

  /* Walk over all triangles. Get the three neighboring triangles.
   * Subdivide the three edges if not done yet. Register the
   * new vertices in the temporary list to mark the edge as divided. */
  std::size_t vertex_amount = verts.size();
  std::size_t face_amount = faces.size() / 3;
  for (std::size_t i = 0; i < face_amount; ++i)
  {
    std::size_t fidx = i * 3;

    for (std::size_t j = 0; j < 3; ++j)
    {
      std::size_t jp1 = (j + 1) % 3;
      std::size_t v1idx = faces[fidx + j];
      std::size_t v2idx = faces[fidx + jp1];
      LoopVClass const& vc1 = this->vclasses[v1idx];
      LoopVClass const& vc2 = this->vclasses[v2idx];

      /* Get more information about the edge. */
      std::size_t face1, face2, v0idx, v3idx;
      bool unique = this->vinfo->get_edge_info(v1idx, v2idx,
          face1, face2, v0idx, v3idx);

      if (face1 == MAX_SIZE_T)
        throw Exception("Subdivision: Face1 was not identified");

      if (!unique)
      {
        // FIXME
        std::cout << "Warning: Non-unique constellation!" << std::endl;
        continue;
      }

      /* Check if the vertex has already been created. */
      if (newvidx[face1 * 3 + j] != MAX_SIZE_T)
        continue;

      /* The edge to be subdivided is v1idx -> v2idx. */
      bool is_smooth = true;
      if (!this->features->empty())
        is_smooth = !this->features->is_feature_edge(v1idx, v2idx);

      /*
       * Now classify the edges into one of:
       * -> smooth edge: get adjacent faces, apply smooth edge mask
       * -> reg crease: apply regular crease subdivision
       * -> non-reg crease: apply non-regular crease subdivision
       */

      if (is_smooth || vc1 == LVC_DART || vc2 == LVC_DART)
      {
        /* Smooth edge subdivision. */
        this->smooth_subdiv(v0idx, v1idx, v2idx, v3idx);
      }
      else if ((vc1 == LVC_CREASE_NONREG && vc2 == LVC_CREASE_NONREG)
          || (vc1 == LVC_CREASE_REG && vc2 == LVC_CREASE_REG)
          || (vc1 == LVC_CREASE_NONREG && vc2 == LVC_CORNER)
          || (vc1 == LVC_CORNER && vc2 == LVC_CREASE_NONREG)
          || (vc1 == LVC_CORNER && vc2 == LVC_CORNER))
      {
        /* Regular crease subdivision. */
        this->regular_crease_subdiv(v1idx, v2idx);
      }
      else if (vc1 == LVC_CREASE_REG
          && (vc2 == LVC_CREASE_NONREG || vc2 == LVC_CORNER))
      {
        /* Non-regular crease subdivision (vc1 -> vc2). */
        this->nonregular_crease_subdiv(v1idx, v2idx);
      }
      else if (vc2 == LVC_CREASE_REG
          && (vc1 == LVC_CREASE_NONREG || vc1 == LVC_CORNER))
      {
        /* Non-regular crease subdivision (vc2 -> vc1). */
        this->nonregular_crease_subdiv(v2idx, v1idx);
      }
      else
      {
        std::cout << "Error: Subdivision rule not determined for classes: "
            << "vc1 = " << vc1 << ", vc2 = " << vc2 << ", smooth: "
            << is_smooth << std::endl;
        throw Exception("Cannot find subdivision rule!");
      }

      /* Now insert the new vertex to the list of new vertices. */
      if (face2 != MAX_SIZE_T)
      {
        std::size_t jf2 = MAX_SIZE_T;
        for (std::size_t k = 0; k < 3; ++k)
          if (faces[face2 * 3 + k] == v2idx)
            jf2 = k;

        if (jf2 == MAX_SIZE_T)
          throw Exception("Subdivision: Error finding face2 index");

        if (newvidx[face1 * 3 + j] != newvidx[face2 * 3 + jf2])
          throw Exception("Subdivision: Edge info does not match");

        newvidx[face2 * 3 + jf2] = verts.size() - 1;
      }

      newvidx[face1 * 3 + j] = verts.size() - 1;
    }
  }

  /* Now calc new position for the original vertices. */
  for (std::size_t i = 0; i < vertex_amount; ++i)
  {
    //if (this->vinfo[i].vclass == VERTEX_CLASS_BORDER)
    //  continue;

    LoopVClass const& vc = this->vclasses[i];

    if (vc == LVC_SMOOTH || vc == LVC_DART)
    {
      VertexInfo::VertexList vlist;
      this->vinfo->get_adjacent_vertices(i, vlist);

      float vn = (float)vlist.size();
      float aterm =  0.375f + 0.25f * std::cos((float)MY_2PI / vn);
      float alpha = 0.625f - aterm * aterm;

      Vec3f allv(0.0f, 0.0f, 0.0f);
      for (std::size_t j = 0; j < vlist.size(); ++j)
        allv += verts[vlist[j]];

      Vec3f& v = verts[i];
      v = v * (1.0f - alpha) + allv * (alpha / vn);
    }
    else if (vc == LVC_CREASE_REG || vc == LVC_CREASE_NONREG)
    {
      if (!this->features->empty())
      {
        FeatureVertexEdges e = this->features->get_edges_for_vertex(i);
        Vec3f const& ve1 = verts[e[0]];
        Vec3f const& ve2 = verts[e[1]];
        Vec3f& v = verts[i];
        v = (v * 6.0f + ve1 + ve2) / 8.0f;
      }
      else
      {
        VertexInfo::VertexList vlist;
        this->vinfo->get_adjacent_vertices(i, vlist);
        Vec3f const& ve1 = verts[vlist[0]];
        Vec3f const& ve2 = verts[vlist[vlist.size() - 1]];
        Vec3f& v = verts[i];
        v = (v * 6.0f + ve1 + ve2) / 8.0f;
      }
    }
    else if (vc == LVC_CORNER)
    {
      /* Nothing to do for a corner, leave this vertex. */
    }
  }

  /* Some data is invalidated. Clear it to save memory. */
  this->vinfo->clear();
  this->vclasses.clear();

  /* Features are also refined using the inserted vertices. */
  if (!this->features->empty())
  {
    this->features->resize(verts.size());
    for (std::size_t i = 0; i < face_amount; ++i)
    {
      std::size_t fidx = i * 3;

      for (std::size_t j = 0; j < 3; ++j)
      {
        std::size_t jp1 = (j + 1) % 3;
        std::size_t v1idx = faces[fidx + j];
        std::size_t v2idx = faces[fidx + jp1];
        std::size_t vidx = newvidx[fidx + j];

        if (vidx == MAX_SIZE_T)
          continue;

        if (!this->features->is_feature_edge(v1idx, v2idx))
          continue;

        this->features->rm_feature_edge(v1idx, v2idx);
        this->features->add_feature_edge(v1idx, vidx);
        this->features->add_feature_edge(v2idx, vidx);
      }
    }
  }

  this->insert_faces(newvidx);
}

/* ---------------------------------------------------------------- */

void
SubdivLoop::classify_vertices (void)
{
  /*
   * Walk over all vertices and classifiy them. Classifications are:
   * - smooth vertex: no incident feature lines
   * - dart vertex: one incident feature line
   * - regular crease vertex: two incident features lines AND
   *   - EITHER simple, degree 6, two non-crease edges at each side
   *   - OR boundary, degree 4
   * - non-regular crease vertex: two incident features AND not regular
   * - corner vertex: three or more incident feature lines
   */

  MeshVertexList const& verts = this->mesh->get_vertices();

  /* Initialize the vertex classes to smooth. */
  this->vclasses.clear();
  this->vclasses.resize(verts.size(), LVC_SMOOTH);

  for (std::size_t i = 0; i < verts.size(); ++i)
  {
    VertexInfo const& vi = this->vinfo[i];
    LoopVClass& vc = this->vclasses[i];

    /* Check if external features are defined. */
    if (this->features->empty())
    {
      /* Only handle borders if there are no features given. */
      if (vi.vclass == VERTEX_CLASS_BORDER)
      {
        if (vi.adj_faces.size() == 3)
          vc = LVC_CREASE_REG;
        else
          vc = LVC_CREASE_NONREG;
      }
    }
    else
    {
      /* Features are given. Get the features for the current vertex. */
      VertexInfo::VertexList vlist;
      this->vinfo->get_adjacent_vertices(i, vlist);
      FeatureVertexEdges e = this->features->get_edges_for_vertex(i, vlist);

      /* Classify the vertex accordingly. */
      if (e.size() == 0)
        vc = LVC_SMOOTH;
      else if (e.size() == 1)
        vc = LVC_DART;
      else if (e.size() >= 3)
        vc = LVC_CORNER;
      else /* e.size() == 2 */
      {
        /* Handle the boundary vertices. */
        if (vi.vclass == VERTEX_CLASS_BORDER)
        {
          /* Check if we have a regular border, and only the boundary
           * is given in the the feature information, tag as regular. */
          if (vi.adj_faces.size() == 3
              && (e[0] == vlist[0] || e[0] == vlist[3])
              && (e[1] == vlist[0] || e[1] == vlist[3]))
            vc = LVC_CREASE_REG;
          else
            vc = LVC_CREASE_NONREG;
        }
        else
        {
          /* Init the vertex as non-regular. Check compare regular creases
           * against the given feature and mark as regular if found. */
          vc = LVC_CREASE_NONREG;
          if (vi.adj_faces.size() == 6)
          {
            for (std::size_t j = 0; j < 6; ++j)
              if (vlist[j] == e[0] && vlist[(j + 3) % 6] == e[1])
              {
                vc = LVC_CREASE_REG;
                break;
              }
          }
        }
      }
    }

    //std::cout << "Classified vertex " << i << " as "
    //    << vclasses[i] << std::endl;
  }
}

/* ---------------------------------------------------------------- */

/* Scheme for smooth edges.    Scheme for the border case.
 *           1
 *         /   \
 *       3 ===== 3                    1 ===== 1
 *         \   /                        \   /
 *           1                            0
 */
void
SubdivLoop::smooth_subdiv (std::size_t v0idx, std::size_t v1idx,
    std::size_t v2idx, std::size_t v3idx)
{
  MeshVertexList& verts = this->mesh->get_vertices();
  Vec3f const& v1 = verts[v1idx];
  Vec3f const& v2 = verts[v2idx];

  Vec3f new_vert;

  if (v3idx == MAX_SIZE_T)
    new_vert = (v1 + v2) / 2.0f;
  else
  {
    Vec3f const& v0 = verts[v0idx];
    Vec3f const& v3 = verts[v3idx];
    new_vert = (v1 * 3.0f + v2 * 3.0f + v0 + v3) / 8.0f;
  }

  verts.push_back(new_vert);
}

/* ---------------------------------------------------------------- */

/* Scheme for regular crease edges.
 *           0
 *         /   \
 *       1 ===== 1
 *         \   /
 *           0
 */
void
SubdivLoop::regular_crease_subdiv (std::size_t v1idx, std::size_t v2idx)
{
  MeshVertexList& verts = this->mesh->get_vertices();
  Vec3f const& v1 = verts[v1idx];
  Vec3f const& v2 = verts[v2idx];

  Vec3f new_vert = (v1 + v2) / 2.0f;
  verts.push_back(new_vert);
}

/* ---------------------------------------------------------------- */

/* Scheme for non-regular crease edges.
 * The regular vertex is weigted with 5, which is v1idx by convention.
 *           0
 *         /   \
 *       5 ===== 3
 *         \   /
 *           0
 */
void
SubdivLoop::nonregular_crease_subdiv (std::size_t v1idx, std::size_t v2idx)
{
  MeshVertexList& verts = this->mesh->get_vertices();
  Vec3f const& v1 = verts[v1idx];
  Vec3f const& v2 = verts[v2idx];

  Vec3f new_vert = (v1 * 5.0f + v2 * 3.0f) / 8.0f;
  verts.push_back(new_vert);
}

REMESHER_NAMESPACE_END
