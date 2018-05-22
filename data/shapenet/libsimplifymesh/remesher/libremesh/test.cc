#include <iostream>
#include <vector>

#include "interface.h"
#include "vec3.h"
#include "matrix3.h"
#include "plane.h"
#include "helpers.h"
#include "modelloader.h"
#include "modelwriter.h"
#include "simplification.h"
#include "straightline3.h"
#include "straightline2.h"
#include "triangle2.h"
#include "triangle3.h"
#include "micropatch.h"
#include "microdelaunay.h"
#include "vertexref.h"
#include "relocation.h"
#include "patch3d.h"
#include "patch2d.h"
#include "exception.h"
#include "meshcleanup.h"
#include "logging.h"
#include "delaunayflips.h"
#include "subdivloop.h"
#include "meshoptimize.h"
#include "polygon2.h"
#include "thread.h"
#include "oversampling.h"
#include "densitytriangle2.h"
#include "permute.h"
#include "meshslicing.h"
#include "incdelaunay.h"
#include "patch2d_fb.h"
#include "averageplane.h"

using namespace Remesher;

int
main (int argc, char** argv)
{
  #if 1
  std::srand(0); // Init randomizer with constant value

  /* Test average plane calculation. */
  Remesher::AveragePlane ap;
  for (std::size_t i = 0; i < 1000; ++i)
  {
    float pos[3];
    for (int j = 0; j < 3; ++j)
        pos[j] = (float)std::rand() / (float)RAND_MAX;

    pos[2] *= 0.5f;
    ap.add_point(Vec3f(pos));
  }

  Plane3f plane = ap.get_average();
  std::cout << "Normal: " << plane.n << std::endl;


  #endif


  #if 0
  Remesher::TriangleMeshPtr mesh(ModelLoader::load_model("/gris/gris-f/home/sfuhrman/offmodels/"));
  Remesher::Patch2dFBPtr p2d = Remesher::Patch2dFB::create();
  p2d
  #endif

  /* Test NAN checks. */
  #if 0
  float val = std::acos(1.5f);
  REMESHER_NAN_CHECK(val)
  #endif


  /* Test complex vertex loop. */
  #if 0
  if (argc != 2)
  {
    std::cout << "Usage: " << argv[0] << " <mesh>" << std::endl;
    return 1;
  }

  Remesher::Interface iface;
  iface.load_model(argv[1]);

  TriangleMeshPtr mesh(iface.get_reference_mesh());
  VertexInfoListPtr vinfo(Remesher::VertexInfoList::create(mesh));

  for (std::size_t i = 0; i < vinfo->size(); ++i)
  {
    VertexInfo::VertexList vlist;
    vinfo->get_adjacent_vertices(i, vlist);
    if (vinfo->is_complex_vertex_loop(vlist))
    {
      std::cout << "Complex vertex loop for vertex " << i << std::endl;
      Patch3dNRingPtr p3d(Patch3dNRing::create(mesh, vinfo));
      p3d->create_patch(i, 2);
      p3d->write_debug_mesh("/tmp/mesh_" + Helpers::get_string(i) + ".off");
    }
  }


  #endif

  /* Test complex vertex loop. */
  #if 0
  Remesher::VertexInfoListPtr vinfo = Remesher::VertexInfoList::create();
  Remesher::VertexInfo::VertexList list;
  list.push_back(0);
  list.push_back(1);
  list.push_back(2);
  list.push_back(3);
  list.push_back(4);
  list.push_back(5);
  vinfo->is_complex_vertex_loop(list);
  #endif

  /* Test is_mesh_edge from vertex info. */
  #if 0
  if (argc != 4)
  {
    std::cout << "Usage: " << argv[0] << " <mesh> <v1> <v2>" << std::endl;
    return 1;
  }

  std::size_t v1 = Helpers::get_sizet_from_string(argv[2]);
  std::size_t v2 = Helpers::get_sizet_from_string(argv[3]);

  Remesher::Interface iface;
  iface.load_model(argv[1]);
  TriangleMeshPtr mesh = iface.get_reference_mesh();
  VertexInfoListPtr vinfo = VertexInfoList::create(mesh);

  Remesher::VertexInfo const& vi = vinfo[v1];

  std::cout << "Is edge " << v1 << " -> " << v2 << ": "
     << vinfo->is_mesh_edge(v1, v2) << std::endl;
  return 0;
  #endif


  /* Test resampling. */
  #if 0
  Remesher::Interface iface;
  iface.load_model((argc == 2 ? argv[1] : ""));
  iface.exec_feature_extraction();
  iface.exec_resampling();
  #endif


  /* Test incremental delaunay. */
  #if 0
  Vec3f v1(0.0f, 0.0f, 0.0f);
  Vec3f v2(1.0f, 0.0f, 0.0f);
  Vec3f v3(0.5f, 5.0f, 0.0f);
  IncDelaunay id(v1, v2, v3);

  std::size_t insert = 1000;

  if (argc == 2)
  {
    insert = Helpers::get_sizet_from_string(argv[1]);
  }

  try
  {
    for (std::size_t i = 0; i < insert; ++i)
    {
      if (i > 0 && i % 500 == 0)
        std::cout << "Inserting sample " << i << "..." << std::endl;

      float rnd1 = (float)std::rand() / (float)RAND_MAX;
      float rnd2 = (float)std::rand() / (float)RAND_MAX;
      float r1 = std::min(rnd1, rnd2);
      float r2 = std::max(rnd1, rnd2);
      bool inserted = id.insert_sample(Vec2f(r1, r2 - r1));
      insert += (int)(!inserted);
    }
  }
  catch (Exception& e)
  {
    std::cout << "Error while inserting samples: " << e.what() << std::endl;
  }

  id.write_debug_mesh("/tmp/sampling.off");
  #endif


  /* Triangle area test. */
  #if 0
  Triangle2f t(Vec2f(0.0f, 0.0f), Vec2f(1.0f, 0.0f), Vec2f(0.0f, 1.0f));
  std::cout << t.get_area() << std::endl;
  #endif

  /* Test area based remeshing. */
  #if 0
  Vec3f c(0.0f, 0.0f, 0.0f);
  Vec3f v[6];
  v[0] = Vec3f(4.0f, 0.0f, 0.0f);
  v[1] = Vec3f(2.0f, 3.4f, 0.0f);
  v[2] = Vec3f(-2.0f, 3.4f, 0.0f);
  v[3] = Vec3f(-4.0f, 0.0f, 0.0f);
  v[4] = Vec3f(-2.0f, -3.4f, 0.0f);
  v[5] = Vec3f(2.0f, -3.4f, 0.0f);

  Matrix2f m;
  m[0] = 0.0f;
  m[1] = 0.0f;
  m[2] = 0.0f;
  m[3] = 0.0f;
  Vec2f b(0.0f, 0.0f);

  for (std::size_t i = 0; i < 6; ++i)
  {
    std::size_t ip1 = (i + 1) % 6;
    float xip1xi = (v[ip1][0] - v[i][0]);
    float yiyip1 = (v[i][1] - v[ip1][1]);

    m[0] += yiyip1 * yiyip1;
    m[1] += xip1xi * yiyip1;
    m[2] += m[1];
    m[3] += xip1xi * xip1xi;

    float xiyip1xip1yi = v[i][0] * v[ip1][1] - v[ip1][0] * v[i][1];
    b[0] += xiyip1xip1yi * yiyip1;
    b[1] += xiyip1xip1yi * xip1xi;
  }

  float det = m[0] * m[3] - m[1] * m[2];
  Matrix2f invm;
  invm[0] = m[3] / det;
  invm[1] = -m[1] / det;
  invm[2] = -m[2] / det;
  invm[3] = m[0] / det;

  std::cout << "M is: " << m << std::endl;
  std::cout << "det is: " << det << std::endl;
  std::cout << "M^-1 is: " << invm << std::endl;
  std::cout << "Solution is: " << (invm * b) << std::endl;


  #endif


  /* Debug loop subdivision. */
  #if 0
  TriangleMeshPtr mesh = ModelLoader::load_off_model
      ("/data/offmodels/testing/testmesh.off");
  mesh->ensure_normals();
  VertexInfoListPtr vinfo = VertexInfoList::create(mesh);

  FeatureEdgesPtr features = FeatureEdges::create();
  features->set_mesh(mesh);
  features->set_vertex_info(vinfo);
  features->extract_features();
  features->add_feature_edge(0, 1);
  features->add_feature_edge(0, 3);
  features->add_feature_edge(0, 4);

  SubdivLoop sl;
  sl.set_mesh(mesh);
  sl.set_feature_edges(features);
  sl.set_vertex_info(vinfo);
  sl.start_subdiv();


  #endif

  /* Test PTT distance. */
  #if 0
  PTTDist ptt;
  ptt.set_triangle(Vec3f(-1.0f, 0.0f, 0.0f), Vec3f(1.0f, 0.0f, 0.0f),
      Vec3f(0.0f, 1.0f, 0.0f));
  std::cout << "Distance: " << ptt.get_squared_dist(Vec3f(0.0f, -1.0f, 0.0f))
      << std::endl;

  #endif


  /* Debugging patch creation. */
  #if 0
  TriangleMeshPtr mesh = ModelLoader::load_off_model
      ("../results/crank_remeshed_2.off");
  mesh->scale_and_center();
  VertexInfoList vinfo(mesh);

  Patch3dTrianglePtr p3d = Patch3dTriangle::create(mesh, vinfo);
  p3d->create_patch(31750, 47229, 47223);

  p3d->write_debug_mesh("../debug/defect_delaunay.off");
  #endif


  /* Test mesh slicing. */
  #if 0
  TriangleMeshPtr mesh = ModelLoader::load_off_model
      ("../gtkremesher/offmodels/cad_models/crank.off");
  mesh->scale_and_center();
  VertexInfoList vinfo(mesh);

  GridSlicer ms;
  ms.set_mesh(mesh);
  //ms.slice(Plane3f(Vec3f(0.0f, 1.0f, 0.0f), 0.0f));
  ms.grid_slice(100);

  ModelWriter::save_off_model("../debug/sliced1.off", mesh);
  #endif


  /* Test permutations. */
  #if 0
  std::vector<double> v;
  v.push_back(0.5);
  v.push_back(1.5);
  v.push_back(2.5);
  v.push_back(3.5);
  v.push_back(4.5);
  v.push_back(5.5);

  std::vector<std::size_t> p;
  p.push_back(0); // 0
  p.push_back(4); // 1
  p.push_back(5); // 2
  p.push_back(2); // 3
  p.push_back(1); // 4
  p.push_back(3); // 5

  Permute<double, std::size_t> pa;
  pa.permute(v, p);

  for (std::size_t i = 0; i < v.size(); ++i)
    std::cout << "Value " << i << ": " << v[i] << std::endl;

  #endif

  /* Testing circumcircle in 3D. */
  #if 0
  Triangle3f t(Vec3f(3, 3, 0), Vec3f(4, 5, 0), Vec3f(3, 7, 0));
  Vec3f center;
  float qradius;
  t.get_circumcircle(qradius, center);
  std::cout << "Q-Radius: " << qradius << ", Center: " << center << std::endl;
  #endif

  /* Testing non-uniform centroid calc. */
  #if 0
  Vec2f v1(0.0f, 0.0f);
  Vec2f v2(8.0f, 8.0f);
  Vec2f v3(0.0f, 16.0f);

  Vec3f vx(v1[0], v2[0], v3[0]);
  Vec3f vy(v1[1], v2[1], v3[1]);

  Vec3f d(4.0f, 2.0f, 4.0f);

  Matrix3f m;
  m[0] = d[0];        m[1] = d[0] / 2.0f; m[2] = d[0] / 2.0f;
  m[3] = d[1] / 2.0f; m[4] = d[1];        m[5] = d[1] / 2.0f;
  m[6] = d[2] / 2.0f; m[7] = d[2] / 2.0f; m[8] = d[2];

  Vec3f f(1.0f / 6.0f, 1.0f / 6.0f, 1.0f / 6.0f);

  float area = v2[0] * v3[1] + v3[0] * v1[1] + v1[0] * v2[1]
      - v1[1] * v2[0] - v2[1] * v3[0] - v3[1] * v1[0];
  area = area / 2.0f;
  float mass = area * (d[0] + d[1] + d[2]) / 3.0f;

  Vec2f c(m.mult(vx).scalar(f), m.mult(vy).scalar(f));
  c *= area / mass;

  std::cout << "Solution: " << c << std::endl;
  std::cout << "Area: " << area << ", Mass: " << mass
      << ", Factor: " << area / mass <<  std::endl;


  DensityTriangle2f td(v1, v2, v3, d);
  std::cout << td.get_centroid() << std::endl;

  #endif

  /* Testing vertex info, edge info */
  #if 0
  TriangleMeshPtr mesh = ModelLoader::load_off_model
      ("../gtkremesher/offmodels/testing/testmesh.off");
  VertexInfoList vinfo(mesh);

  std::size_t v1idx = 5;
  std::size_t v2idx = 4;

  std::size_t face1 = MAX_SIZE_T;
  std::size_t face2 = MAX_SIZE_T;
  std::size_t v0idx = MAX_SIZE_T;
  std::size_t v3idx = MAX_SIZE_T;
  vinfo.get_edge_info(v1idx, v2idx, face1, face2, v0idx, v3idx);

  std::cout << "Face1: " << face1 << ", Face2: " << face2
      << ", v0: " << v0idx << ", v3: " << v3idx << std::endl;


  #endif


  /* Testing clipping of polygons. */
  #if 0
  Polygon2 p;
  p.push_back(Vec2f(0.0f, 0.0f));
  p.push_back(Vec2f(2.0f, 2.0f));
  p.push_back(Vec2f(0.0f, 4.0f));
  p.push_back(Vec2f(-2.0f, 2.0f));
  p.clip_with_line(Vec2f(1, 4), Vec2f(-2, 1));
  for (std::size_t i = 0; i < p.size(); ++i)
    std::cout << p[i] << ", ";
  std::cout << std::endl;

  #endif


  /* Testing feature edges. */
  #if 0
  TriangleMeshPtr mesh = ModelLoader::load_off_model
      ("../gtkremesher/offmodels/remeshes/box_r.off");
  mesh->ensure_normals(true, false);
  VertexInfoList vinfo(mesh);
  FeatureEdges features;
  features.set_mesh(mesh);
  features.set_vertex_info(vinfo);
  features.extract_features();
  FeatureEdges::VertexEdges ve = features.get_edges_for_vertex(94);
  for (std::size_t i = 0; i < ve.size(); ++i)
    std::cout << "Feature edge: " << ve[i] << std::endl;
  #endif


  /* Testing area/centroid calculation of a polygon. */
  #if 0
  std::vector<Vec2f> poly;
  float area;
  Vec2f centroid;

  /* A square. */
  poly.push_back(Vec2f(0, 0));
  poly.push_back(Vec2f(2, 0));
  poly.push_back(Vec2f(2, 2));
  poly.push_back(Vec2f(0, 2));
  calc_area_and_centroid(poly, area, centroid);
  std::cout << "Area: " << area << ", Centroid: " << centroid << std::endl;

  /* A more complex polygon. */
  poly.clear();
  poly.push_back(Vec2f(0, 0));
  poly.push_back(Vec2f(2, 0));
  poly.push_back(Vec2f(4, 1.5f));
  poly.push_back(Vec2f(2, 3));
  poly.push_back(Vec2f(0, 3));
  poly.push_back(Vec2f(0, 2));
  poly.push_back(Vec2f(0, 1));
  calc_area_and_centroid(poly, area, centroid);
  std::cout << "Area: " << area << ", Centroid: " << centroid << std::endl;
  #endif


  /* String normalize test. */
  #if 0
  std::string str = "cool test string";
  Helpers::normalize_string(str);
  std::cout << str << std::endl;
  str = "cool    test         string";
  Helpers::normalize_string(str);
  std::cout << str << std::endl;
  str = "cool\t\ttest\t   \t\t string";
  Helpers::normalize_string(str);
  std::cout << str << std::endl;
  #endif


  /* Oversampling testing. */
  #if 0
  TriangleMeshPtr mesh = ModelLoader::load_off_model
      ("../gtkremesher/offmodels/animals/duck.off");
  mesh->ensure_normals(true, true);

  OversamplingConf conf;
  conf.vertex_limit = 1000;

  Oversampling o;
  o.set_mesh(mesh, true);
  o.set_config(conf);
  try
  {
    o.start_oversampling();
  }
  catch (Exception& e)
  {
    std::cout << "Error: " << e << std::endl;
  }

  ModelWriter::save_off_model("../debug/oversampled1.off", o.get_mesh());
  #endif


  /* Straight line testings. */
  #if 0
  {
    StraightLine2f line(Vec2f(0.0f, 0.0f), Vec2f(3.0f, 3.0f));
    std::cout << line.edge_equation(Vec2f(1.0f, 1.0f)) << std::endl;
  }
  #endif
  #if 0
  {
    StraightLine2f line(Vec2f(0.0f, 0.0f), Vec2f(0.0f, 3.0f));
    std::cout << line.point_qdist(Vec2f(1.0f, 0.0f)) << std::endl;
    std::cout << line.point_qdist(Vec2f(1.0f, 1.0f)) << std::endl;
    std::cout << line.point_qdist(Vec2f(1.0f, -1.0f)) << std::endl;
  }

  {
    StraightLine2f line(Vec2f(0.0f, 0.0f), Vec2f(3.0f, 3.0f));
    std::cout << line.point_qdist(Vec2f(3.0f, 0.0f)) << std::endl;
  }

  {
    StraightLine2f line(Vec2f(0.0f, 0.0f), Vec2f(-1.0f, 2.0f));
    std::cout << line.point_qdist(Vec2f(0.50f, -1.0f)) << std::endl;
  }
  #endif


  #if 0
  TriangleMeshPtr mesh = ModelLoader::load_model
      ("../gtkremesher/offmodels/cad_models/fandisk.m");
  VertexInfoListPtr vinfo = VertexInfoList::create(mesh);

  MeshOptimize o;
  o.set_mesh(mesh);
  o.set_vertex_info(vinfo);
  //o.randomize();
  //vinfo.calc_for_mesh(mesh);
  //o.set_vertex_info(vinfo);
  o.optimize(24);

  ModelWriter::save_off_model("../debug/optimized_fandisk2.off", mesh);
  #endif



  #if 0
  LOG(LOG_INFO) << "Some info message and " << 213 << " numbers";
  LOG(LOG_DEBUG) << "Some debugging message and " << 213 << " numbers";
  LOG(LOG_DEBUG3) << "Some geek debugging and " << 213 << " numbers";
  #endif


  #if 0
  TriangleMeshPtr mesh = ModelLoader::load_off_model("../gtkremesher/offmodels/testmesh3.off");
  MeshCleanup mc(mesh);
  mc.cleanup_duplicated(0.15f);
  ModelWriter::save_off_model("../debug/cleaned.off", mesh);

  #endif


  #if 0
  TriangleMeshPtr mesh = ModelLoader::load_off_model("/data/offmodels/animals/bunny.off");
  VertexInfoList vinfo(mesh);
  Patch3dTrianglePtr patch3d = Patch3dTriangle::create(mesh, vinfo);
  patch3d->create_patch(10045, 10373, 10050);
  Patch2dPtr patch = Patch2d::create((Patch3dPtr)patch3d);

  Patch2d::VertexMapping const& vmap = patch->get_vertex_mapping();
  Patch2d::FaceMapping const& fmap = patch->get_face_mapping();

  std::cout << "Vertex mapping: " << std::endl;
  for (std::size_t i = 0; i < vmap.size(); ++i)
    std::cout << i << " -> " << vmap[i] << std::endl;
  std::cout << std::endl;

  std::cout << "Face mapping: " << std::endl;
  for (std::size_t i = 0; i < fmap.size(); ++i)
    std::cout << i << " -> " << fmap[i] << std::endl;
  std::cout << std::endl;

  // Does not work, vertices are reordered during parametrization
  std::cout << patch->lookup_vertex(6456) << std::endl;

  return 0;
  #endif


  #if 0
  TriangleMeshPtr mesh = ModelLoader::load_off_model("../gtkremesher/offmodels/animals/cow.off");
  VertexInfoList vinfo(mesh);

  DelaunayFlips d3d;
  d3d.set_mesh(mesh);
  d3d.set_vertex_info(vinfo);
  try
  { d3d.restore_delaunay(); }
  catch (Exception& e)
  {
    std::cout << "Excetion: " << e << std::endl;
  }

  ModelWriter::save_off_model("../debug/edgeflip.off", mesh);
  #endif


  #if 0
  TriangleMeshPtr mesh = ModelLoader::load_off_model("../gtkremesher/offmodels/bunny_7k.off");
  VertexInfoList vinfo(mesh);

  VertexRef vr1, vr2, vr3;
  vr1.face = 682;
  vr2.face = 840;
  vr3.face = 716;
  vr1.bary = Vec2f(0.33f, 0.33f);
  vr2.bary = Vec2f(0.33f, 0.33f);
  vr3.bary = Vec2f(0.33f, 0.33f);

  RelocationPtr reloc = Relocation::create(mesh, vinfo);
  VertexRef nr = reloc->relocate(vr1, vr2, vr3, Vec2f(0.33, 0.33));
  Triangle3f tri(mesh, nr.face);
  Vec3f pnt = tri.get_point_for_bary(nr.bary);

  std::cout << "New face: " << nr.face << ", bary: " << nr.bary
      << ", new point: " << pnt << std::endl;

  Triangle3f t1(mesh, vr1.face);
  Triangle3f t2(mesh, vr2.face);
  Triangle3f t3(mesh, vr3.face);
  Vec3f v1 = t1.get_point_for_bary(vr1.bary);
  Vec3f v2 = t2.get_point_for_bary(vr2.bary);
  Vec3f v3 = t3.get_point_for_bary(vr3.bary);

  //Patch3dTrianglePtr p = Patch3dTriangle::create(mesh, vinfo);
  //p->create_patch(5755, 5714, 5955); // Bunny invalid
  //p->create_patch(529, 136, 598); // Bunny7k medium patch
  //p->create_patch(753, 715, 844); // Bunny7k small patch
  //p->create_patch(38670, 2064, 31159); // Bunny big patch
  //p->create_patch(126103, 149336, 1398); // Elephant huge patch
  //p->create_patch(582, 582, 873);

  //Patch2dPtr p2 = Patch2d::create(p);
  #endif


  //MicroDelaunay mp;
  #if 0 // Delaunay
  mp.set_center_vertex(Vec3f(0, 0, 0));
  mp.append_adjacent_vertex(Vec3f(1, -5, 0));
  mp.append_adjacent_vertex(Vec3f(2, 0, 0));
  mp.append_adjacent_vertex(Vec3f(1, 2, 0));
  mp.append_adjacent_vertex(Vec3f(-1, 2, 0));
  mp.append_adjacent_vertex(Vec3f(-2, 0, 0));
  mp.append_adjacent_vertex(Vec3f(-1, -5, 0));
  #endif
  #if 0 // Delaunay
  mp.set_center_vertex(Vec3f(0, 0, 0));
  mp.append_adjacent_vertex(Vec3f(2, -2, 0));
  mp.append_adjacent_vertex(Vec3f(2, 2, 0));
  mp.append_adjacent_vertex(Vec3f(-2, 2, 0));
  mp.append_adjacent_vertex(Vec3f(-2, -2, 0));
  #endif
  #if 0 // Non-delaunay
  mp.set_center_vertex(Vec3f(0, 0, 0));
  mp.append_adjacent_vertex(Vec3f(2, -2, 0));
  mp.append_adjacent_vertex(Vec3f(2, 2, 0));
  mp.append_adjacent_vertex(Vec3f(-2, 2, 0));
  mp.append_adjacent_vertex(Vec3f(-2, -2, 0));
  mp.append_adjacent_vertex(Vec3f(0, -3.9999999, 0));
  #endif

  #if 0
  mp.calculate_patch();
  ModelWriter::save_off_model("patch1.off", mp.get_debug_mesh());

  std::cout << "Is delaunay: " << mp.is_delaunay() << std::endl;
  if (mp.is_delaunay())
  {
    Vec2f centroid = mp.get_voronoi_centroid();
    mp.debug_bary_coords(centroid);
    std::cout << "Centroid: " << centroid << std::endl;

    MicroPatch::AdjacentVertices2D& adj2d = mp.get_adjacent_vertices();
    for (unsigned int i = 0; i < adj2d.size(); ++i)
      adj2d[i] -= centroid;
  }
  ModelWriter::save_off_model("patch2.off", mp.get_debug_mesh());
  #endif

  #if 0
  Triangle2f t;
  t.set(Vec2f(-1.0f, 0.0f), Vec2f(1.0f, 0.0f), Vec2f(0.0f, 1.0f));
  std::cout << t.get_bary_coords(Vec2f(0.75f, -0.0f)) << std::endl;
  #endif

  #if 0
  TriangleMeshPtr mesh = ModelLoader::load_off_model
      ("../gtkremesher/offmodels/testmesh.off");
  SimplificationConf conf;
  conf.vertex_limit = 1;
  Simplification s;
  s.set_mesh(mesh, true);
  s.set_config(conf);
  s.start_algorithm();
  #endif

  #if 0
  TriangleMeshPtr mesh = ModelLoader::load_off_model("dragon_head.off");

  float min[3] = { 999999, 99999, 99999 };
  float max[3] = { -999999, -999999, -999999 };
  MeshVertexList const& vertices = mesh->get_vertices();
  for (unsigned int i = 0; i < vertices.size(); ++i)
  {
    if (vertices[i][0] > max[0]) max[0] = vertices[i][0];
    if (vertices[i][1] > max[1]) max[1] = vertices[i][1];
    if (vertices[i][2] > max[2]) max[2] = vertices[i][2];
    if (vertices[i][0] < min[0]) min[0] = vertices[i][0];
    if (vertices[i][1] < min[1]) min[1] = vertices[i][1];
    if (vertices[i][2] < min[2]) min[2] = vertices[i][2];
  }

  std::vector<Vec3f> new_vertices;
  float threshold = min[1] + 1.5f / 3.0f * (max[1] - min[1]);
  for (unsigned int i = 0; i < vertices.size(); ++i)
  {
    if (vertices[i][1] > threshold)
      new_vertices.push_back(vertices[i]);
  }

  mesh->get_vertices() = new_vertices;
  mesh->get_faces().clear();
  mesh->get_face_normals().clear();
  ModelWriter::save_off_model("dragon_head.off", mesh);
  #endif

  return 0;
}
