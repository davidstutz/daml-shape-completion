#include <fstream> //REMOVE
#include <iostream>
#include <queue>
#include <set>

#include "modelwriter.h"
#include "patch3d.h"

REMESHER_NAMESPACE_BEGIN

void
Patch3d::write_debug_mesh (std::string const& filename)
{
  typedef std::set<std::size_t> MeshFaceSet;
  TriangleMeshPtr save = TriangleMesh::create();
  MeshFaceList& sfl = save->get_faces();
  MeshFaceList const& ofl = this->mesh->get_faces();

  save->get_vertices() = this->mesh->get_vertices();

  for (MeshFaceSet::iterator i = this->faces.begin();
      i != this->faces.end(); ++i)
  {
    sfl.push_back(ofl[*i * 3 + 0]);
    sfl.push_back(ofl[*i * 3 + 1]);
    sfl.push_back(ofl[*i * 3 + 2]);
  }

  ModelWriter::save_off_model(filename, save);
}

/* ================================================================ */

void
Patch3dTriangle::create_patch (std::size_t f1, std::size_t f2, std::size_t f3)
{
  this->add_faces(f1, f2, f3);
  this->add_faces(f2, f1, f3);
  this->add_faces(f3, f1, f2);
}

/* ---------------------------------------------------------------- */

void
Patch3dTriangle::add_faces (std::size_t f1, std::size_t f2, std::size_t f3)
{
  bool seen_f2 = (f1 == f2);
  bool seen_f3 = (f1 == f3);

  MeshFaceList const& faces = this->mesh->get_faces();

  /* A patch is created for face f1. After faces have been
   * collected, the faces are transfered to the real patch. */
  MeshFaceSet patch;
  patch.insert(f1);

  /* Queue for the breath-first search. */
  std::queue<std::size_t> q;
  q.push(f1);

  bool growing_big = false;
  while (!q.empty() && !(seen_f2 && seen_f3))
  {
    //if (q.size() > 500)
    //  std::cout << "Warning: Patch is growing big!" << std::endl;

    /* Handle the oldest triangle in the queue. */
    std::size_t ti = q.front();
    q.pop();

    /* Add all adjacent triangles for this current triangle. This is
     * done by adding faces of adjacent triangles to a temporary set. */
    MeshFaceSet tset;
    for (int i = 0; i < 3; ++i)
    {
      VertexInfo const& vi = this->vinfo[faces[ti * 3 + i]];
      for (std::size_t j = 0; j < vi.adj_faces.size(); ++j)
        tset.insert(vi.adj_faces[j]);
    }

    /* Insert the new faces to the patch. Add each successfully
     * inserted face to the queue for later processing. */
    for (MeshFaceSet::iterator iter = tset.begin(); iter != tset.end(); ++iter)
    {
      if (*iter == f2)
        seen_f2 = true;

      if (*iter == f3)
        seen_f3 = true;

      bool inserted = patch.insert(*iter).second;
      if (inserted)
        q.push(*iter);
    }

    if (patch.size() > 10000 && !growing_big)
    {
      std::cout << "Warngin: Patch is growing big!" << std::endl;
      growing_big = true;
    }
  }

  /* After all faces has been collected, they are transfered to the patch. */
  this->faces.insert(patch.begin(), patch.end());
}

/* ---------------------------------------------------------------- */

bool
Patch3dTriangle::add_face (std::size_t f)
{
  /* TODO check for topology. */
  return this->faces.insert(f).second;
}

/* ================================================================ */

void
Patch3dNRing::create_patch (std::size_t index, std::size_t n_rings)
{
  typedef std::set<std::size_t> VertexSet;
  //typedef std::map<std::size_t, std::size_t> vmap;
  VertexSet vset;
  vset.insert(index); // the 0-ring

  for (std::size_t i = 0; i < n_rings; ++i)
  {
    std::vector<std::size_t> new_verts;
    for (VertexSet::iterator iter = vset.begin(); iter != vset.end(); ++iter)
    {
      VertexInfo::VertexList vlist;
      this->vinfo->get_adjacent_vertices(*iter, vlist);
      for (std::size_t j = 0; j < vlist.size(); ++j)
        new_verts.push_back(vlist[j]);
    }
    vset.insert(new_verts.begin(), new_verts.end());
  }

  MeshFaceList const& mfaces = this->mesh->get_faces();
  std::size_t face_amount = mfaces.size() / 3;

  for (std::size_t i = 0; i < face_amount; ++i)
  {
    std::size_t v[3] = { mfaces[i * 3 + 0],
        mfaces[i * 3 + 1], mfaces[i * 3 + 2] };

    bool has_all = true;
    for (std::size_t j = 0; j < 3; ++j)
      if (vset.find(v[j]) == vset.end())
        has_all = false;
    if (has_all)
      this->faces.insert(i);     
  }

  std::cout << "Created " << n_rings << "-ring patch with " << vset.size()
      << " vertices and " << this->faces.size() << " faces" << std::endl;

  #if 0
  std::cout << "Writing vertices..." << std::flush;
  std::ofstream out("/tmp/mesh_verts_xxx.off");
  out << "OFF" << std::endl;
  out << vset.size() << " 0 0" << std::endl;
  for (VertexSet::iterator i = vset.begin(); i != vset.end(); ++i)
    out << this->mesh->get_vertices()[*i][0] << " "
        << this->mesh->get_vertices()[*i][1] << " "
        << this->mesh->get_vertices()[*i][2] << std::endl;
  out.close();
  std::cout << "done" << std::endl;
  #endif
}

REMESHER_NAMESPACE_END
