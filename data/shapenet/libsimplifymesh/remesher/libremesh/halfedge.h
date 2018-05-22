#ifndef REMESHER_HALFEDGE_HEADER
#define REMESHER_HALFEDGE_HEADER

#include "defines.h"

REMESHER_NAMESPACE_BEGIN

struct HEEdge; // Forward declaration of half edges

struct HEVertex
{
  std::size_t id; // ID of the vertex
  HEEdge* edge; // One of the emantating half edges

  HEVertex (void) : id(MAX_SIZE_T), edge(0) {}
  HEVertex (std::size_t id) : id(id), edge(0) {}
  HEVertex (std::size_t id, HEEdge* edge) : id(id), edge(edge) {}
};

struct HEFace
{
  HEEdge* edge; // One of the adjacent half edges

  HEFace (void) : edge(0) {}
  HEFace (HEEdge* edge) : edge(edge) {}
};

struct HEEdge
{
  HEVertex* vert; // Vertex at the end of the half edge
  HEEdge* pair; // Oppositely oriented adjacent half-edge
  HEFace* face; // Face adjacent to the half edge
  HEEdge* next; // Next half edge in ccw orientation around the face

  HEEdge (void) : vert(0), pair(0), face(0), next(0) {}
  HEEdge (HEVertex* v, HEEdge* p, HEFace* f, HEEdge* n)
    : vert(v), pair(p), face(f), next(n) {}
};

/* ---------------------------------------------------------------- */

class HalfEdge
{
  public:
    typedef std::vector<HEVertex*> VertexList;
    typedef std::vector<HEEdge*> EdgeList;
    typedef std::vector<HEFace*> FaceList;

  protected:
    VertexList he_verts;
    EdgeList he_edges;
    FaceList he_faces;

  public:
    HalfEdge (void);
    ~HalfEdge (void);

    void push_back (HEVertex* vert);
    void push_back (HEEdge* edge);
    void push_back (HEFace* face);

    VertexList const& get_vertices (void) const;
    EdgeList const& get_edges (void) const;
    FaceList const& get_faces (void) const;
};

/* ---------------------------------------------------------------- */

inline
HalfEdge::HalfEdge (void)
{
}

inline
HalfEdge::~HalfEdge (void)
{
  for (std::size_t i = 0; i < this->he_verts.size(); ++i)
    delete this->he_verts[i];
  for (std::size_t i = 0; i < this->he_faces.size(); ++i)
    delete this->he_faces[i];
  for (std::size_t i = 0; i < this->he_edges.size(); ++i)
    delete this->he_edges[i];
}

inline void
HalfEdge::push_back (HEVertex* vert)
{
  this->he_verts.push_back(vert);
}

inline void
HalfEdge::push_back (HEEdge* edge)
{
  this->he_edges.push_back(edge);
}

inline void
HalfEdge::push_back (HEFace* face)
{
  this->he_faces.push_back(face);
}

inline HalfEdge::VertexList const&
HalfEdge::get_vertices (void) const
{
  return this->he_verts;
}

inline HalfEdge::EdgeList const&
HalfEdge::get_edges (void) const
{
  return this->he_edges;
}

inline HalfEdge::FaceList const&
HalfEdge::get_faces (void) const
{
  return this->he_faces;
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_HALFEDGE_HEADER */
