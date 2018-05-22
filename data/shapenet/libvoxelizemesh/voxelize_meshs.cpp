#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cfloat>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <boost/filesystem.hpp>
#include <json/json.h>
#include <H5Cpp.h>

// OpenMP
#include <omp.h>

#include "box_triangle/aabb_triangle_overlap.h"
#include "box_ray/vector3.h"
#include "box_ray/ray.h"
#include "box_ray/box.h"
#include "triangle_point/vec.h"
#include "triangle_point/makelevelset3.h"
#include "triangle_ray/raytri.h"

/** \brief Simple struct representing arbitrary boxes. */
struct Box {
  /** \brief Minimum coordinates. */
  Eigen::Vector3f min;
  /** \brief Maximum coordinates. */
  Eigen::Vector3f max;

  /** \brief Constructor.
   */
  Box() {
    this->min = Eigen::Vector3f::Zero();
    this->max = Eigen::Vector3f::Zero();
  }

  /** \brief Constructor; construct a box from minimum and maximum coordinates.
   * \param[in] min minimum coordinates
   * \param[in] max maximum coordinates
   */
  Box(const Eigen::Vector3f &min, const Eigen::Vector3f &max) {
    this->min = min;
    this->max = max;
  }

  /** \brief Check if a point is inside the box.
   * \param[in] point point to check
   * \return contained
   */
  bool contains(const Eigen::Vector3f &point) {
    bool contains = true;

    for (int d = 0; d < 3; d++) {
      if (point(d) < this->min(d) || point(d) > this->max(d)) {
        contains = false;
      }
    }

    return contains;
  }

  /** \brief Scale the box.
   * \param[in] scale scale per axis
   */
  void scale(const Eigen::Vector3f &scale) {
    for (int d = 0; d < 3; d++) {
      min(d) *= scale(d);
      max(d) *= scale(d);
    }
  }

  /** \brief Scale the box.
   * \param[in] scale scale over all axes
   */
  void scale(float scale) {
    this->scale(Eigen::Vector3f(scale, scale, scale));
  }

  /** \brief Scale the box.
   * \param[in] translation translation per axis
   */
  void translate(const Eigen::Vector3f &translation) {
    for (int d = 0; d < 3; d++) {
      min(d) += translation(d);
      max(d) += translation(d);
    }
  }
};

/** \brief Test box ray intersection.
 * \param[in] box box
 * \param[in] ray ray to intersect with
 * \return intersects
 */
bool box_ray_intersection(const Box &box, const Eigen::Vector3f &ray) {
  for (int d = 0; d < 3; d++) {
    assert(box.min(d) < box.max(d));
  }

  // convert to used data structure
  bri::Vector3 _o(0, 0, 0);
  bri::Vector3 _d(ray(0), ray(1), ray(2));
  bri::Ray _ray(_o, _d);
  bri::Vector3 _min(box.min(0), box.min(1), box.min(2));
  bri::Vector3 _max(box.max(0), box.max(1), box.max(2));
  bri::Box _box(_min, _max);

  // interval for valid its is [0,1] as the ray will
  // be the full-length vector to an observed point
  return _box.intersect(_ray, 0, 0.99);
}

/** \brief Test box ray intersection.
 * \param[in] box box
 * \param[in] origin origin point
 * \param[in] dest destination point
 * \return intersects
 */
bool box_ray_intersection(const Box &box, const Eigen::Vector3f &origin, const Eigen::Vector3f &dest) {
  for (int d = 0; d < 3; d++) {
    assert(box.min(d) < box.max(d));
  }

  // convert to used data structure
  bri::Vector3 _o(origin(0), origin(1), origin(2));
  bri::Vector3 _d(dest(0), dest(1), dest(2));
  bri::Ray _ray(_o, _d);
  bri::Vector3 _min(box.min(0), box.min(1), box.min(2));
  bri::Vector3 _max(box.max(0), box.max(1), box.max(2));
  bri::Box _box(_min, _max);

  // interval for valid its is [0,1] as the ray will
  // be the full-length vector to an observed point
  return _box.intersect(_ray, 0, 0.99);
}

/** \brief Compute triangle box intersection.
 * \param[in] voxel voxel to test intersection for
 * \param[in] v1 first vertex
 * \param[in] v2 second vertex
 * \param[in] v3 third vertex
 * \return intersects
 */
bool triangle_box_intersection(const Box &voxel, const Eigen::Vector3f &v1, const Eigen::Vector3f &v2, const Eigen::Vector3f &v3) {
  float half_size[3] = {
    (voxel.max(0) - voxel.min(0))/2.,
    (voxel.max(1) - voxel.min(1))/2.,
    (voxel.max(2) - voxel.min(2))/2.
  };

  float center[3] = {
    voxel.max(0) - half_size[0],
    voxel.max(1) - half_size[1],
    voxel.max(2) - half_size[2]
  };

  float vertices[3][3] = {{v1(0), v1(1), v1(2)}, {v2(0), v2(1), v2(2)}, {v3(0), v3(1), v3(2)}};
  return triBoxOverlap(center, half_size, vertices);
}

/** \brief Compute triangle point distance and corresponding closest point.
 * \param[in] point point
 * \param[in] v1 first vertex
 * \param[in] v2 second vertex
 * \param[in] v3 third vertex
 * \param[out] ray corresponding closest point
 * \return distance
 */
float triangle_point_distance(const Eigen::Vector3f &point, const Eigen::Vector3f &v1, const Eigen::Vector3f &v2, const Eigen::Vector3f &v3,
    Eigen::Vector3f &closest_point) {

  Vec3f x0(point.data());
  Vec3f x1(v1.data());
  Vec3f x2(v2.data());
  Vec3f x3(v3.data());

  Vec3f r(0);
  float distance = point_triangle_distance_field(x0, x1, x2, x3, r);

  for (int d = 0; d < 3; d++) {
    closest_point(d) = r[d];
  }

  return distance;
}

/** \brief Test triangle ray intersection.
 * \param[in] origin origin of ray
 * \param[in] dest destination of ray
 * \param[in] v1 first vertex
 * \param[in] v2 second vertex
 * \param[in] v3 third vertex
 * \return intersects
 */
bool triangle_ray_intersection(const Eigen::Vector3f &origin, const Eigen::Vector3f &dest,
    const Eigen::Vector3f &v1, const Eigen::Vector3f &v2, const Eigen::Vector3f &v3, float &t) {

  double _origin[3] = {origin(0), origin(1), origin(2)};
  double _dir[3] = {dest(0) - origin(0), dest(1) - origin(1), dest(2) - origin(2)};
  double _v1[3] = {v1(0), v1(1), v1(2)};
  double _v2[3] = {v2(0), v2(1), v2(2)};
  double _v3[3] = {v3(0), v3(1), v3(2)};

  // t is the distance, u and v are barycentric coordinates
  // http://fileadmin.cs.lth.se/cs/personal/tomas_akenine-moller/code/raytri_tam.pdf
  double _t, u, v;
  int success = intersect_triangle(_origin, _dir, _v1, _v2, _v3, &_t, &u, &v);
  t = _t;

  if (success) {
    return true;
  }

  return false;
}

/** \brief Simple struct encapsulating a 3D bounding box. */
struct BoundingBox {
  /** \brief Size of the bounding box. */
  Eigen::Vector3f size;
  /** \brief Translation (i.e. center) of bounding box. */
  Eigen::Vector3f translation;
  /** \brief Rotation of bounding box (in radians per axis). */
  Eigen::Vector3f rotation;
  /** \brief Meta information for writing. */
  std::string meta;

  /** \brief Constructor.
   */
  BoundingBox() {
    this->size = Eigen::Vector3f::Zero();
    this->translation = Eigen::Vector3f::Zero();
    this->rotation = Eigen::Vector3f::Zero();
    this->meta = "";
  }

  /** \brief Assignment operator.
   * \param[in] bounding_box bounding box to assign
   * \return this
   */
  BoundingBox& operator=(const BoundingBox &bounding_box) {
    this->size = bounding_box.size;
    this->translation = bounding_box.translation;
    this->rotation = bounding_box.rotation;
    this->meta = bounding_box.meta;
    return *this;
  }

  /** \brief Convert to regular box.
   * \param[out] box box to convert to
   */
  void to_box(Box &box) const {
    for (int d = 0 ; d < 3; d++) {
      box.min(d) = this->translation(d) - this->size(d)/2;
      box.max(d) = this->translation(d) + this->size(d)/2;
    }
  }

  /** \brief Given the angle in radians, construct a rotation matrix around the y-axis.
   * \param[in] radians angle in radians
   * \param[out] rotation rotation matrix
   */
  static void rotation_matrix_y(const float radians, Eigen::Matrix3f &rotation) {
    rotation = Eigen::Matrix3f::Zero();

    rotation(0, 0) = std::cos(radians); rotation(0, 2) = std::sin(radians);
    rotation(1, 1) = 1;
    rotation(2, 0) = -std::sin(radians); rotation(2, 2) = std::cos(radians);
  }

  /** \brief Check if a point lies inside the bounding box.
   * \param[in] point point to check
   * \return contained in bounding box
   */
  bool contains(const Eigen::Vector3f &point) const {
    Eigen::Vector3f transformed = point - this->translation;
    Eigen::Matrix3f rotation;
    BoundingBox::rotation_matrix_y(this->rotation(1), rotation);
    transformed = rotation*transformed;

    bool contained = true;
    for (int d = 0; d < 3; d++) {
      if (transformed(d) > this->size(d)/2 || transformed(d) < -this->size(d)/2) {
        contained = false;
      }
    }

    return contained;
  }

  /** \brief Get largest scale of bounding box.
   * \return scale
   */
  float scale() const {
    float scale = 0;

    for (int d = 0; d < 3; d++) {
      if (this->size(d) > scale) {
        scale = this->size(d);
      }
    }

    assert(scale > 0);
    return scale;
  }
};

/** \brief Just encapsulating vertices and faces. */
class Mesh {
public:
  /** \brief Empty constructor. */
  Mesh() {

  }

  /** \brief Add a vertex.
   * \param[in] vertex vertex to add
   */
  void add_vertex(Eigen::Vector3f& vertex) {
    this->vertices.push_back(vertex);
  }

  /** \brief Get the number of vertices.
   * \return number of vertices
   */
  int num_vertices() {
    return static_cast<int>(this->vertices.size());
  }

  /** \brief Add a face.
   * \param[in] face face to add
   */
  void add_face(Eigen::Vector3i& face) {
    this->faces.push_back(face);
  }

  /** \brief Get the number of faces.
   * \return number of faces
   */
  int num_faces() {
    return static_cast<int>(this->faces.size());
  }

  /** \brief Translate the mesh.
   * \param[in] translation translation vector
   */
  void translate(const Eigen::Vector3f& translation) {
    for (int v = 0; v < this->num_vertices(); ++v) {
      for (int i = 0; i < 3; ++i) {
        this->vertices[v](i) += translation(i);
      }
    }
  }

  /** \brief Scale the mesh.
   * \param[in] scale scale vector
   */
  void scale(const Eigen::Vector3f& scale) {
    for (int v = 0; v < this->num_vertices(); ++v) {
      for (int i = 0; i < 3; ++i) {
        this->vertices[v](i) *= scale(i);
      }
    }
  }

  /** \brief Voxelize the given mesh into a dense volume.
   * \param[in] mesh mesh to voxelize
   * \param[in] n batch index in dense
   * \param[in] dense dense pre-initialized volume
   */
  void voxelize_UNSCALED(int n, Eigen::Tensor<float, 4, Eigen::RowMajor>& sdf) {

    int height = sdf.dimension(1);
    int width = sdf.dimension(2);
    int depth = sdf.dimension(3);

    for (int h = 0; h < height; ++h) {
      for (int w = 0; w < width; ++w) {
        for (int d = 0; d < depth; ++d) {
          sdf(n, h, w, d) = FLT_MAX;

          // [Data] 99% percentile dimensions:
          // [Data]   4.750000 1.870000 1.970000
          // [Data]   2.540102 1        1.053475
          // [Data]   56       22       24

          float padding_factor = 1.2;
          float width_factor = (padding_factor*2.545454545)/2.545454545; // = padding_factor
          float height_factor = padding_factor/2.545454545;
          float depth_factor = (padding_factor*1.090909091)/2.545454545;

          // the box corresponding to this voxel
          Eigen::Vector3f min(width_factor*static_cast<float>(w)/width, height_factor*static_cast<float>(h)/height,
            depth_factor*static_cast<float>(d)/depth);
          Eigen::Vector3f max(width_factor*static_cast<float>(w + 1)/width, height_factor*static_cast<float>(h + 1)/height,
            depth_factor*static_cast<float>(d + 1)/depth);

          Box voxel_box(min, max);
          voxel_box.translate(Eigen::Vector3f(-width_factor/2, -height_factor/2, -depth_factor/2));
          Eigen::Vector3f center((voxel_box.max(0) + voxel_box.min(0))/2., (voxel_box.max(1) + voxel_box.min(1))/2., (voxel_box.max(2) + voxel_box.min(2))/2.);

          // count number of intersections.
          int num_intersect = 0;
          for (unsigned int f = 0; f < this->faces.size(); ++f) {

            Eigen::Vector3f v1 = this->vertices[this->faces[f](0)];
            Eigen::Vector3f v2 = this->vertices[this->faces[f](1)];
            Eigen::Vector3f v3 = this->vertices[this->faces[f](2)];

            Eigen::Vector3f closest_point;
            triangle_point_distance(center, v1, v2, v3, closest_point);
            float distance = (center - closest_point).norm();

            if (distance < sdf(n, h, w, d)) {
              sdf(n, h, w, d) = distance;
            }

            bool intersect = triangle_ray_intersection(center, Eigen::Vector3f(0, 0, 0), v1, v2, v3, distance);

            if (intersect && distance >= 0) {
              num_intersect++;
            }
          }

          if (num_intersect%2 == 1) {
            sdf(n, h, w, d) *= -1;
          }
        }
      }
    }
  }

  /** \brief Voxelize the given mesh into a dense volume.
   * \param[in] mesh mesh to voxelize
   * \param[in] n batch index in dense
   * \param[in] dense dense pre-initialized volume
   */
  void voxelize_SCALED(int n, Eigen::Tensor<float, 4, Eigen::RowMajor>& sdf) {

    int height = sdf.dimension(1);
    int width = sdf.dimension(2);
    int depth = sdf.dimension(3);

    #pragma omp parallel
    {
      #pragma omp for
      for (int i = 0; i < height*width*depth; i++) {
        int d = i%depth;
        int w = (i/depth)%width;
        int h = (i/depth)/width;
        sdf(n, h, w, d) = FLT_MAX;

        // the box corresponding to this voxel
        Eigen::Vector3f min(static_cast<float>(w), static_cast<float>(h), static_cast<float>(d));
        Eigen::Vector3f max(static_cast<float>(w + 1), static_cast<float>(h + 1), static_cast<float>(d + 1));

        Box voxel_box(min, max);
        //voxel_box.translate(Eigen::Vector3f(-static_cast<float>(width)/2, -static_cast<float>(height)/2, -static_cast<float>(depth)/2));
        Eigen::Vector3f center((voxel_box.max(0) + voxel_box.min(0))/2., (voxel_box.max(1) + voxel_box.min(1))/2., (voxel_box.max(2) + voxel_box.min(2))/2.);

        // count number of intersections.
        int num_intersect = 0;
        for (unsigned int f = 0; f < this->faces.size(); ++f) {

          Eigen::Vector3f v1 = this->vertices[this->faces[f](0)];
          Eigen::Vector3f v2 = this->vertices[this->faces[f](1)];
          Eigen::Vector3f v3 = this->vertices[this->faces[f](2)];

          Eigen::Vector3f closest_point;
          triangle_point_distance(center, v1, v2, v3, closest_point);
          float distance = (center - closest_point).norm();

          if (distance < sdf(n, h, w, d)) {
            sdf(n, h, w, d) = distance;
          }

          bool intersect = triangle_ray_intersection(center, Eigen::Vector3f(0, 0, 0), v1, v2, v3, distance);

          if (intersect && distance >= 0) {
            num_intersect++;
          }
        }

        if (num_intersect%2 == 1) {
          sdf(n, h, w, d) *= -1;
        }
      }
    }
  }

  /** \brief Sample points from the mesh
   * \param[in] mesh mesh to sample from
   * \param[in] n batch index in points
   * \param[in] points pre-initialized tensor holding points
   */
  void sample(int n, Eigen::Tensor<float, 3, Eigen::RowMajor>& points) {

    // The number of points to sample.
    const int N = points.dimension(1);

    // Stores the areas of faces.
    std::vector<float> areas(this->faces.size());
    float sum = 0;

    // Build a probability distribution over faces.
    for (int f = 0; f < this->faces.size(); f++) {
      Eigen::Vector3f a = this->vertices[this->faces[f][0]];
      Eigen::Vector3f b = this->vertices[this->faces[f][1]];
      Eigen::Vector3f c = this->vertices[this->faces[f][2]];

      // Angle between a->b and a->c.
      Eigen::Vector3f ab = b - a;
      Eigen::Vector3f ac = c - a;
      float cos_angle = ab.dot(ac)/(ab.norm()*ac.norm());
      float angle = std::acos(cos_angle);

      // Compute triangle area.
      float area = std::max(0., 0.5*ab.norm()*ac.norm()*std::sin(angle));
      //std::cout << area << " " << std::pow(area, 1./4.) << " " << angle << " " << ab.norm() << " " << ac.norm() << " " << std::sin(angle) << std::endl;

      // Accumulate.
      //area = std::sqrt(area);
      areas[f] = area;
      sum += area;
      //areas.push_back(1);
      //sum += 1;
    }

    //std::cout << sum << std::endl;
    assert(sum > 1e-8);

    for (int f = 0; f < this->faces.size(); f++) {
      //std::cout << areas[f] << " ";
      areas[f] /= sum;
      //std::cout << areas[f] << std::endl;
    }

    std::vector<float> cum_areas(areas.size());
    cum_areas[0] = areas[0];

    for (int f = 1; f < this->faces.size(); f++) {
      cum_areas[f] = areas[f] + cum_areas[f - 1];
    }

    for (int i = 0; i < N; i++) {
      float r = static_cast<float>(std::rand())/static_cast<float>(RAND_MAX);
      int face = 0;

      while (r > cum_areas[face + 1] && face < this->faces.size() - 1) {
        face++;
      }

      assert(face >= 0 && face < this->faces.size());
      //int face = std::rand()%this->faces.size();

      float r1 = 0;
      float r2 = 0;
      do {
        r1 = static_cast<float>(std::rand())/static_cast<float>(RAND_MAX);
        r2 = static_cast<float>(std::rand())/static_cast<float>(RAND_MAX);
      }
      while (r1 + r2 > 1.f);

      int s = std::rand()%3;
      //std::cout << face << " " << areas[face] << std::endl;

      Eigen::Vector3f a = this->vertices[this->faces[face](s)];
      Eigen::Vector3f b = this->vertices[this->faces[face]((s + 1)%3)];
      Eigen::Vector3f c = this->vertices[this->faces[face]((s + 2)%3)];

      Eigen::Vector3f ab = b - a;
      Eigen::Vector3f ac = c - a;

      Eigen::Vector3f point = a + r1*ab + r2*ac;
      points(n, i, 0) = point(0);
      points(n, i, 1) = point(1);
      points(n, i, 2) = point(2);
    }
  }

  /** \brief Write mesh to OFF file.
   * \param[in] filepath path to OFF file to write
   * \return success
   */
  bool to_off(const std::string filepath) {
    std::ofstream* out = new std::ofstream(filepath, std::ofstream::out);
    if (!static_cast<bool>(out)) {
      return false;
    }

    (*out) << "OFF" << std::endl;
    (*out) << this->vertices.size() << " " << this->faces.size() << " 0" << std::endl;

    for (unsigned int v = 0; v < this->vertices.size(); v++) {
      (*out) << this->vertices[v](0) << " " << this->vertices[v](1) << " " << this->vertices[v](2) << std::endl;
    }

    for (unsigned int f = 0; f < this->faces.size(); f++) {
      (*out) << "3 " << this->faces[f](0) << " " << this->faces[f](1) << " " << this->faces[f](2) << std::endl;
    }

    out->close();
    delete out;

    return true;
  }

private:

  /** \brief Vertices as (x,y,z)-vectors. */
  std::vector<Eigen::Vector3f> vertices;

  /** \brief Faces as list of vertex indices. */
  std::vector<Eigen::Vector3i> faces;
};

/** \brief Reading an off file and returning the vertices x, y, z coordinates and the
 * face indices.
 * \param[in] filepath path to the OFF file
 * \param[out] mesh read mesh with vertices and faces
 * \return success
 */
bool read_off(const std::string filepath, Mesh& mesh) {

  std::ifstream* file = new std::ifstream(filepath.c_str());
  std::string line;
  std::stringstream ss;
  int line_nb = 0;

  std::getline(*file, line);
  ++line_nb;

  if (line != "off" && line != "OFF") {
    std::cout << "[Error] Invalid header: \"" << line << "\", " << filepath << std::endl;
    return false;
  }

  size_t n_edges;
  std::getline(*file, line);
  ++line_nb;

  int n_vertices;
  int n_faces;
  ss << line;
  ss >> n_vertices;
  ss >> n_faces;
  ss >> n_edges;

  for (size_t v = 0; v < n_vertices; ++v) {
    std::getline(*file, line);
    ++line_nb;

    ss.clear();
    ss.str("");

    Eigen::Vector3f vertex;
    ss << line;
    ss >> vertex(0);
    ss >> vertex(1);
    ss >> vertex(2);

    mesh.add_vertex(vertex);
  }

  size_t n;
  for (size_t f = 0; f < n_faces; ++f) {
    std::getline(*file, line);
    ++line_nb;

    ss.clear();
    ss.str("");

    size_t n;
    ss << line;
    ss >> n;

    if(n != 3) {
      std::cout << "[Error] Not a triangle (" << n << " points) at " << (line_nb - 1) << std::endl;
      return false;
    }

    Eigen::Vector3i face;
    ss >> face(0);
    ss >> face(1);
    ss >> face(2);

    mesh.add_face(face);
  }

  if (n_vertices != mesh.num_vertices()) {
    std::cout << "[Error] Number of vertices in header differs from actual number of vertices." << std::endl;
    return false;
  }

  if (n_faces != mesh.num_faces()) {
    std::cout << "[Error] Number of faces in header differs from actual number of faces." << std::endl;
    return false;
  }

  file->close();
  delete file;

  return true;
}

/** \brief Write the given set of volumes to h5 file.
 * \param[in] filepath h5 file to write
 * \param[in] n number of volumes
 * \param[in] height height of volumes
 * \param[in] width width of volumes
 * \param[in] depth depth of volumes
 * \param[in] dense volume data
 */
bool write_hdf5(const std::string filepath, Eigen::Tensor<float, 4, Eigen::RowMajor>& dense) {

  try {

    /*
     * Turn off the auto-printing when failure occurs so that we can
     * handle the errors appropriately
     */
    H5::Exception::dontPrint();

    /*
     * Create a new file using H5F_ACC_TRUNC access,
     * default file creation properties, and default file
     * access properties.
     */
    H5::H5File file(filepath, H5F_ACC_TRUNC);

    /*
     * Define the size of the array and create the data space for fixed
     * size dataset.
     */
    hsize_t rank = 4;
    hsize_t dimsf[rank];
    dimsf[0] = dense.dimension(0);
    dimsf[1] = dense.dimension(1);
    dimsf[2] = dense.dimension(2);
    dimsf[3] = dense.dimension(3);
    H5::DataSpace dataspace(rank, dimsf);

    /*
     * Define datatype for the data in the file.
     * We will store little endian INT numbers.
     */
    H5::IntType datatype(H5::PredType::NATIVE_FLOAT);
    datatype.setOrder(H5T_ORDER_LE);

    /*
     * Create a new dataset within the file using defined dataspace and
     * datatype and default dataset creation properties.
     */
    H5::DataSet dataset = file.createDataSet("tensor", datatype, dataspace);

    /*
     * Write the data to the dataset using default memory space, file
     * space, and transfer properties.
     */
    float* data = static_cast<float*>(dense.data());
    dataset.write(data, H5::PredType::NATIVE_FLOAT);
  }  // end of try block

  // catch failure caused by the H5File operations
  catch(H5::FileIException error) {
    error.printError();
    return false;
  }

  // catch failure caused by the DataSet operations
  catch(H5::DataSetIException error) {
    error.printError();
    return false;
  }

  // catch failure caused by the DataSpace operations
  catch(H5::DataSpaceIException error) {
    error.printError();
    return false;
  }

  // catch failure caused by the DataSpace operations
  catch(H5::DataTypeIException error) {
    error.printError();
    return false;
  }

  return true;
}

/** \brief Write the given set of volumes to h5 file.
 * \param[in] filepath h5 file to write
 * \param[in] n number of volumes
 * \param[in] height height of volumes
 * \param[in] width width of volumes
 * \param[in] depth depth of volumes
 * \param[in] dense volume data
 */
bool write_hdf5(const std::string filepath, Eigen::Tensor<float, 3, Eigen::RowMajor>& dense) {

  try {

    /*
     * Turn off the auto-printing when failure occurs so that we can
     * handle the errors appropriately
     */
    H5::Exception::dontPrint();

    /*
     * Create a new file using H5F_ACC_TRUNC access,
     * default file creation properties, and default file
     * access properties.
     */
    H5::H5File file(filepath, H5F_ACC_TRUNC);

    /*
     * Define the size of the array and create the data space for fixed
     * size dataset.
     */
    hsize_t rank = 3;
    hsize_t dimsf[rank];
    dimsf[0] = dense.dimension(0);
    dimsf[1] = dense.dimension(1);
    dimsf[2] = dense.dimension(2);
    H5::DataSpace dataspace(rank, dimsf);

    /*
     * Define datatype for the data in the file.
     * We will store little endian INT numbers.
     */
    H5::IntType datatype(H5::PredType::NATIVE_FLOAT);
    datatype.setOrder(H5T_ORDER_LE);

    /*
     * Create a new dataset within the file using defined dataspace and
     * datatype and default dataset creation properties.
     */
    H5::DataSet dataset = file.createDataSet("tensor", datatype, dataspace);

    /*
     * Write the data to the dataset using default memory space, file
     * space, and transfer properties.
     */
    float* data = static_cast<float*>(dense.data());
    dataset.write(data, H5::PredType::NATIVE_FLOAT);
  }  // end of try block

  // catch failure caused by the H5File operations
  catch(H5::FileIException error) {
    error.printError();
    return false;
  }

  // catch failure caused by the DataSet operations
  catch(H5::DataSetIException error) {
    error.printError();
    return false;
  }

  // catch failure caused by the DataSpace operations
  catch(H5::DataSpaceIException error) {
    error.printError();
    return false;
  }

  // catch failure caused by the DataSpace operations
  catch(H5::DataTypeIException error) {
    error.printError();
    return false;
  }

  return true;
}

/** \brief Write the given set of volumes to h5 file.
 * \param[in] filepath h5 file to write
 * \param[in] n number of volumes
 * \param[in] height height of volumes
 * \param[in] width width of volumes
 * \param[in] depth depth of volumes
 * \param[in] dense volume data
 */
bool write_hdf5(const std::string filepath, Eigen::Tensor<int, 4, Eigen::RowMajor>& dense) {

  try {

    /*
     * Turn off the auto-printing when failure occurs so that we can
     * handle the errors appropriately
     */
    H5::Exception::dontPrint();

    /*
     * Create a new file using H5F_ACC_TRUNC access,
     * default file creation properties, and default file
     * access properties.
     */
    H5::H5File file(filepath, H5F_ACC_TRUNC);

    /*
     * Define the size of the array and create the data space for fixed
     * size dataset.
     */
    hsize_t rank = 4;
    hsize_t dimsf[rank];
    dimsf[0] = dense.dimension(0);
    dimsf[1] = dense.dimension(1);
    dimsf[2] = dense.dimension(2);
    dimsf[3] = dense.dimension(3);
    H5::DataSpace dataspace(rank, dimsf);

    /*
     * Define datatype for the data in the file.
     * We will store little endian INT numbers.
     */
    H5::IntType datatype(H5::PredType::NATIVE_INT);
    datatype.setOrder(H5T_ORDER_LE);

    /*
     * Create a new dataset within the file using defined dataspace and
     * datatype and default dataset creation properties.
     */
    H5::DataSet dataset = file.createDataSet("tensor", datatype, dataspace);

    /*
     * Write the data to the dataset using default memory space, file
     * space, and transfer properties.
     */
    int* data = static_cast<int*>(dense.data());
    dataset.write(data, H5::PredType::NATIVE_INT);
  }  // end of try block

  // catch failure caused by the H5File operations
  catch(H5::FileIException error) {
    error.printError();
    return false;
  }

  // catch failure caused by the DataSet operations
  catch(H5::DataSetIException error) {
    error.printError();
    return false;
  }

  // catch failure caused by the DataSpace operations
  catch(H5::DataSpaceIException error) {
    error.printError();
    return false;
  }

  // catch failure caused by the DataSpace operations
  catch(H5::DataTypeIException error) {
    error.printError();
    return false;
  }

  return true;
}

/** \brief Read all files in a directory matching the given extension.
 * \param[in] directory path to directory
 * \param[out] files read file paths
 * \param[in] extension extension to filter for
 */
void read_directory(const std::string directory, std::map<int, std::string>& files, const std::string extension = ".off") {

  boost::filesystem::path dir(directory);
  boost::filesystem::directory_iterator end;

  files.clear();
  for (boost::filesystem::directory_iterator it(dir); it != end; ++it) {
    if (it->path().extension().string() == extension) {
      int number = std::stoi(it->path().filename().string());
      files.insert(std::pair<int, std::string>(number, it->path().string()));
    }
  }
}

/** \brief Small helper for safely retrieving values from Json root element. */
class JSON {
public:

  /** \brief Constructor.
   * \param[in] root JSON root element
   */
  JSON(Json::Value root) {
    this->root = root;
  }

  /** \brief Get the value of the given key as int.
   * \param[in] key key to retrieve
   * \return value
   */
  int get_int(std::string key) {
    if (!this->root.isMember(key)) {
      std::cout << "[Error] key " << key << " not found" << std::endl;
      exit(1);
    }

    return this->root[key].asInt();
  }

  /** \brief Get the value of the given key as float.
   * \param[in] key key to retrieve
   * \return value
   */
  float get_float(std::string key) {
    if (!this->root.isMember(key)) {
      std::cout << "[Error] key " << key << " not found" << std::endl;
      exit(1);
    }

    return this->root[key].asFloat();
  }

  /** \brief Get the value of the given key as string.
   * \param[in] key key to retrieve
   * \return value
   */
  std::string get_string(std::string key) {
    if (!this->root.isMember(key)) {
      std::cout << "[Error] key " << key << " not found" << std::endl;
      exit(1);
    }

    return this->root[key].asString();
  }

private:

  /** \brief Root element. */
  Json::Value root;
};

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cout << "[Error] Usage: voxelize_cuboids config.json" << std::endl;
    exit(1);
  }

  Json::Value root;
  Json::Reader reader;

  std::string config_file = argv[1];

  if (!boost::filesystem::is_regular_file(boost::filesystem::path(config_file))) {
    std::cout << "[Error] Config file not found: " << config_file << std::endl;
    exit(1);
  }

  std::ifstream json_file(config_file, std::ifstream::in | std::ifstream::binary);
  reader.parse(json_file, root, false);
  JSON json(root);

  // Unsafe!
  int height = json.get_int("height");
  int width = json.get_int("width");
  int depth = json.get_int("depth");

  int image_height = json.get_int("image_height");
  int image_width = json.get_int("image_width");

  std::string suffix = json.get_string("suffix");
  int multiplier = json.get_int("multiplier");

  std::string sdf_file = json.get_string("sdf_file");
  if (sdf_file.empty()) {
    std::cout << "[Error] Read invalid output file" << std::endl;
    exit(1);
  }

  sdf_file = sdf_file + "_" + std::to_string(multiplier) + "_"
      + std::to_string(image_height) + "x" + std::to_string(image_width) + "_"
      + std::to_string(height) + "x" + std::to_string(width) + "x" + std::to_string(depth)
      + suffix + ".h5";

  std::string output_file = json.get_string("output_file");
  if (output_file.empty()) {
    std::cout << "[Error] Read invalid output file" << std::endl;
    exit(1);
  }

  output_file = output_file + "_" + std::to_string(multiplier) + "_"
      + std::to_string(image_height) + "x" + std::to_string(image_width) + "_"
      + std::to_string(height) + "x" + std::to_string(width) + "x" + std::to_string(depth)
      + suffix + ".h5";

  std::string off_directory = json.get_string("off_dir");
  if (off_directory.empty()) {
    std::cout << "[Error] Read invalid off directory" << std::endl;
    exit(1);
  }

  off_directory = off_directory + "_" + std::to_string(multiplier) + "_"
      + std::to_string(image_height) + "x" + std::to_string(image_width) + "_"
      + std::to_string(height) + "x" + std::to_string(width) + "x" + std::to_string(depth) + suffix;

  std::string off_gt_directory = json.get_string("off_gt_dir");
  if (off_gt_directory.empty()) {
    std::cout << "[Error] Read invalid off gt directory" << std::endl;
    exit(1);
  }

  off_gt_directory = off_gt_directory + "_" + std::to_string(multiplier) + "_"
      + std::to_string(image_height) + "x" + std::to_string(image_width) + "_"
      + std::to_string(height) + "x" + std::to_string(width) + "x" + std::to_string(depth) + suffix;

  if (!boost::filesystem::is_directory(boost::filesystem::path(off_gt_directory))) {
    boost::filesystem::create_directories(boost::filesystem::path(off_gt_directory));
  }

  std::string point_file = json.get_string("point_file");
  if (point_file.empty()) {
    std::cout << "[Error] Read invalid point file" << std::endl;
    exit(1);
  }

  point_file = point_file + "_" + std::to_string(multiplier) + "_"
      + std::to_string(image_height) + "x" + std::to_string(image_width) + "_"
      + std::to_string(height) + "x" + std::to_string(width) + "x" + std::to_string(depth)
      + suffix + ".h5";

  if (!boost::filesystem::is_directory(boost::filesystem::path(off_directory))) {
    std::cout << "[Error] OFF directory " << off_directory << " not found" << std::endl;
    exit(1);
  }

  std::map<int, std::string> files;
  read_directory(off_directory, files);

  int N = files.size();
  Eigen::Tensor<float, 4, Eigen::RowMajor> sdf(N, height, width, depth);
  sdf.setZero();

  int N_points = json.get_int("n_points");
  Eigen::Tensor<float, 3, Eigen::RowMajor> points(N, N_points, 3);
  points.setZero();

  std::vector<int> indices;
  for (std::map<int, std::string>::iterator it = files.begin(); it != files.end(); it++) {
    indices.push_back(it->first);
  }

  int n = 0;
  #pragma omp parallel
  {
    #pragma omp for
    for (unsigned int i = 0; i< indices.size(); i++) {
      int n = indices[i];

      Mesh mesh;
      bool success = read_off(files[n], mesh);

      if (!success) {
        std::cout << "[Error] could not read " << files[n] << std::endl;
        exit(1);
      }

      float padding_factor = 1 + json.get_float("padding");
      float scale_factor = static_cast<float>(width)/padding_factor;
      mesh.scale(Eigen::Vector3f(scale_factor, scale_factor, scale_factor));
      mesh.translate(Eigen::Vector3f(static_cast<float>(width)/2, static_cast<float>(height)/2, static_cast<float>(depth)/2));

      std::cout << "[Data] sampling mesh " << n << "/" << files.size() << " " << files[n] << std::endl;
      mesh.sample(i, points);

      std::string off_file = off_gt_directory + "/" + std::to_string(i) + ".off";
      std::cout << "[Data] writing " << off_file << std::endl;
      mesh.to_off(off_file);

      if (!success) {
        std::cout << "[Error] error reading " << files[n] << std::endl;
        exit(1);
      }

      std::cout << "[Data] voxelizing mesh " << n << "/" << files.size() << " " << files[n] << std::endl;
      mesh.voxelize_SCALED(i, sdf);

      n++;
    }
  }

  bool success = write_hdf5(sdf_file, sdf);

  if (success) {
    std::cout <<"[Data] wrote " << sdf_file << std::endl;
  }
  else {
    std::cout << "[Error] error writing " << sdf_file << std::endl;
  }

  success = write_hdf5(point_file, points);

  if (success) {
    std::cout <<"[Data] wrote " << point_file << std::endl;
  }
  else {
    std::cout << "[Error] error writing " << point_file << std::endl;
  }

  Eigen::Tensor<int, 4, Eigen::RowMajor> occ(N, height, width, depth);
  occ.setZero();

  for (int n = 0; n < N; n++) {
    for (int h = 0; h < height; h++) {
      for (int w = 0; w < width; w++) {
        for (int d = 0; d < depth; d++) {
          if (sdf(n, h, w, d) <= 0) {
            occ(n, h, w, d) = 1;
          }
        }
      }
    }
  }

  success = write_hdf5(output_file, occ);

  if (success) {
    std::cout <<"[Data] wrote " << output_file << std::endl;
  }
  else {
    std::cout << "[Error] error writing " << output_file << std::endl;
  }

  exit(0);
}