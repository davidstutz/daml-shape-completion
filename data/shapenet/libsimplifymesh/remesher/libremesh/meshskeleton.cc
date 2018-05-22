#include <set>
#include <iostream>

#include "exception.h"
#include "meshskeleton.h"

REMESHER_NAMESPACE_BEGIN

void
MeshSkeleton::extract (void)
{
  if (this->features.get() == 0)
    throw Exception("MeshSkeleton: No features set");
  if (this->vinfo.get() == 0)
    throw Exception("MeshSkeleton: No vertex info set");

  this->clear();
  this->num_closed = 0;
  this->num_open = 0;

  /*
   * Strategy is as follows:
   * - Find a feature edge E not yet "consumed"
   * - Initially mark backbone as closed
   * - Remember first edge
   * - Iterate...
   *   - First check if second vertex of E is inner-edge
   *   - If not inner-edge, mark backbone as open, set resume flag
   *   - If inner-edge but second vertex equals remembered edge, bail out
   *   - Add the vertex at the end of the backbone list, consume edge
   *   - Advance edge in direction of the second vertex
   * - If resume flag is set, restart with first edge and reverse iterate
   *   - During that pass, push vertices to the front of the backbone
   */

  /* Create a copy of the feature edges for consumption. */
  FeatureEdgesPtr edge_pool(this->features->create_copy());

  /* Start finding a not-yet consumed feature edge. */
  for (std::size_t i = 0; i < this->features->size(); ++i)
  {
    FeatureVertexEdges& edges = this->features[i];
    for (std::size_t j = 0; j < edges.size(); ++j)
      if (edge_pool->is_feature_edge(i, edges[j]))
        this->collect(FeatureEdge(i, edges[j]), edge_pool);
  }
}

/* ---------------------------------------------------------------- */

void
MeshSkeleton::collect (FeatureEdge edge, FeatureEdgesPtr pool)
{
  /* Create a new backbone. */
  this->push_back(MeshBackbone());
  MeshBackbone& b = this->back();
  b.closed = true;
  b.verts.push_back(edge.first);

  FeatureEdge iter = edge;
  while (true)
  {
    /* Insert current (second) vertex. */
    if (b.closed)
      b.verts.push_back(iter.second);
    else
      b.verts.push_front(iter.second);

    /* Consume edge. */
    pool->rm_feature_edge(iter.first, iter.second);

    /* Advance edge. */
    bool has_next = this->advance_edge(iter);

    /* If there's no next edge during reverse mode, it's done. */
    if (!has_next && !b.closed)
      break;

    /* On a corner vertex, mark backbone as open and reverse. */
    if (!has_next)
    {
      b.closed = false;
      iter.first = edge.second;
      iter.second = edge.first;
      b.verts.pop_front();
      continue;
    }

    /* If walked around a closed backbone, bail out. */
    if (iter.second == edge.first)
    {
      pool->rm_feature_edge(iter.first, iter.second);
      break;
    }
  }

  if (b.closed)
    this->num_closed += 1;
  else
    this->num_open += 1;
}

/* ---------------------------------------------------------------- */

bool
MeshSkeleton::advance_edge (FeatureEdge& edge) const
{
  std::size_t i = edge.first;
  std::size_t j = edge.second;

  FeatureVertexEdges const& edges = this->features->get_edges_for_vertex(j);

  /* If it's a corner vertex, indicate that there is no next edge. */
  if (edges.size() != 2)
    return false;

  edge.first = j;
  edge.second = (i == edges[0] ? edges[1] : edges[0]);

  if (this->max_angle != 0)
  {
    MeshVertexList const& verts = this->mesh->get_vertices();
    Vec3f const v[3] = { verts[i], verts[j], verts[edge.second] };
    Vec3f const e[2] = { v[1] - v[0], v[2] - v[1] };
    float cosangle = e[0].scalar(e[1]) / std::sqrt
        (e[0].qlength() * e[1].qlength());
    if (cosangle < std::cos(this->max_angle))
        return false;
  }

  return true;
}

/* ---------------------------------------------------------------- */

std::size_t
MeshSkeleton::get_corner_amount (void) const
{
  /* Collect start and end vertex of each open backbone (no dups!). */
  std::set<std::size_t> corners;
  for (std::size_t i = 0; i < this->size(); ++i)
  {
    if (!this->at(i).closed)
    {
      corners.insert(this->at(i).verts.front());
      corners.insert(this->at(i).verts.back());
    }
  }
  return corners.size();
}

/* ---------------------------------------------------------------- */

void
MeshSkeleton::debug_print (void) const
{
  typedef std::list<std::size_t>::const_iterator VertexIter;
  for (std::size_t i = 0; i < this->size(); ++i)
  {
    MeshBackbone const& b(this->at(i));
    std::cout << (b.closed ? "Closed" : "Open") << " backbone: ";
    for (VertexIter j = b.verts.begin(); j != b.verts.end(); ++j)
      std::cout << *j << " ";
    std::cout << std::endl;
  }
}

REMESHER_NAMESPACE_END
