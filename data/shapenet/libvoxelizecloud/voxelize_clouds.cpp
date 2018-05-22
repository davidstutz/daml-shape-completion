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

#include "box_ray/vector3.h"
#include "box_ray/ray.h"
#include "box_ray/box.h"

/** \brief Class representing a point cloud in 3D. */
class PointCloud {
public:
  /** \brief Constructor. */
  PointCloud() {

  }

  /** \brief Constructor.
   * \param[in] points points
   */
  PointCloud(const std::vector<Eigen::Vector3f> &points) {
    this->points = points;
    this->colors = std::vector<Eigen::Vector3i>(this->points.size(), Eigen::Vector3i(0, 0, 0));
  }

  /** \brief Copy constructor.
   * \param[in] point_cloud point cloud to copy
   */
  PointCloud(const PointCloud &point_cloud) {
    this->points.clear();
    this->colors.clear();

    for (unsigned int i = 0; i < point_cloud.points.size(); i++) {
      this->points.push_back(point_cloud.points[i]);
      this->colors.push_back(point_cloud.colors[i]);
    }
  }

  /** \brief Destructor. */
  ~PointCloud() {

  }

  /** \brief Read point cloud from binary file.
   * \param[in] filepath path to binary file
   * \param[out] point_cloud read point cloud
   */
  static bool from_bin(const std::string &filepath, PointCloud &point_cloud) {

    // start with clear point cloud
    point_cloud.points.clear();
    point_cloud.colors.clear();

    // allocate 4 MB buffer (only ~130*4*4 KB are needed)
    int buffer_size = 1000000;
    float* data = new float[buffer_size];

    // important to start with clean memory to avoid reading points
    // from previous files!
    for (int i = 0; i < buffer_size; i++) {
      data[i] = 0;
    }

    // pointers
    float* px = data + 0;
    float* py = data + 1;
    float* pz = data + 2;

    // load point cloud
    std::ifstream* in = new std::ifstream(filepath, std::ifstream::in | std::ifstream::binary);

    if (!static_cast<bool>(*in)) {
      return false;
    }

    // https://stackoverflow.com/questions/19614581/reading-floating-numbers-from-bin-file-continuosly-and-outputting-in-console-win
    in->seekg(0, std::ios::end);
    int n_points = in->tellg()/sizeof(float)/4; // 4 floats per point
    //std::cout << "[Data] " << n_points << " points" << std::endl;
    in->seekg(0, std::ios::beg);

    in->read(reinterpret_cast<char*>(data), n_points*sizeof(float));

    for (int i = 0; i < n_points; i++) {
      Eigen::Vector3f point;
      // bounding boxes are in camera coordinate system
      // x = right, y = down, z = forward
      // but velodyne point cloud
      // is in velodyne system
      // x = forward, y = left, z = up
      point(0) = -*py;
      point(1) = *pz; // note the flip!
      point(2) = *px;

      if (std::abs(point(0)) > 1e-8 || std::abs(point(1)) > 1e-8 || std::abs(point(2)) > 1e-8) {
        if (std::abs(point(0)) < 1e8 && std::abs(point(1)) < 1e8 && std::abs(point(2)) < 1e8) {
          point_cloud.add_point(point);
        }
      }

      px += 4;
      py += 4;
      pz += 4;
    }

    in->close();
    delete in;
    delete[] data;

    return true;
  }

  /** \brief Add a point to the point cloud.
   * \param[in] point point to add
   */
  void add_point(const Eigen::Vector3f &point) {
    this->points.push_back(point);
    this->colors.push_back(Eigen::Vector3i::Zero());
  }

  /** \brief Add a colored point to the point cloud.
   * \param[in] point point to add
   * \param[in] color color of the point
   */
  void add_point(const Eigen::Vector3f &point, const Eigen::Vector3i &color) {
    this->points.push_back(point);
    this->colors.push_back(color);
  }

  /** \brief Get number of points.
   * \return number of points
   */
  unsigned int num_points() const {
    return this->points.size();
  }

  /** \brief Write point cloud to txt file.
   * \param[in] filepath path to file
   * \return success
   */
  bool to_txt(const std::string &filepath) {
    std::ofstream* out = new std::ofstream(filepath, std::ofstream::out);
    if (!static_cast<bool>(*out)) {
      return false;
    }

    (*out) << this->points.size() << std::endl;
    for (unsigned int i = 0; i < this->points.size(); i++) {
     (*out) << this->points[i](0) << " " << this->points[i](1) << " " << this->points[i](2) << std::endl;
    }

    out->close();
    delete out;

    return true;
  }

  /** \brief Write point cloud to txt file.
   * \param[in] filepath path to file
   * \return success
   */
  bool to_ply(const std::string &filepath) {
    std::ofstream* out = new std::ofstream(filepath, std::ofstream::out);
    if (!static_cast<bool>(*out)) {
      return false;
    }

    if (this->points.size() != this->colors.size()) {
      return false;
    }

    (*out) << "ply" << std::endl;
    (*out) << "format ascii 1.0" << std::endl;
    (*out) << "element vertex " << this->points.size() << std::endl;
    (*out) << "property float32 x" << std::endl;
    (*out) << "property float32 y" << std::endl;
    (*out) << "property float32 z" << std::endl;
    (*out) << "property uchar red" << std::endl;
    (*out) << "property uchar green" << std::endl;
    (*out) << "property uchar blue" << std::endl;
    (*out) << "end_header" << std::endl;

    for (unsigned int i = 0; i < this->points.size(); i++) {
      //out << this->points[i](0) << " " << this->points[i](1) << " " << this->points[i](2) << std::endl;
      (*out) << this->points[i](0) << " " << this->points[i](1) << " " << this->points[i](2) << " "
        << this->colors[i](0) << " " << this->colors[i](1) << " " << this->colors[i](2) << std::endl;
    }

    out->close();
    delete out;

    return true;
  }

  /** \brief Scale the point cloud.
   * \param[in] scale scales for each axis
   */
  void scale(const Eigen::Vector3f &scale) {
    for (unsigned int i = 0; i < this->points.size(); i++) {
      for (int d = 0; d < 3; d++) {
        this->points[i](d) *= scale(d);
      }
    }
  }

  /** \brief Scale the point cloud.
   * \param[in] scale overlal scale for all axes
   */
  void scale(float scale) {
    this->scale(Eigen::Vector3f(scale, scale, scale));
  }

  /** \brief Rotate the point cloud around the origin.
   * \param[in] rotation rotation matrix
   */
  void rotate(const Eigen::Matrix3f &rotation) {
    for (unsigned int i = 0; i < this->points.size(); i++) {
      this->points[i] = rotation*this->points[i];
    }
  }

  /** \brief Translate the point cloud.
   * \brief translation translation vector
   */
  void translate(const Eigen::Vector3f &translation) {
    for (unsigned int i = 0; i < this->points.size(); i++) {
      this->points[i] += translation;
    }
  }

private:
  /** \brief The points of the point cloud. */
  std::vector<Eigen::Vector3f> points;
  /** \brief Colors of the points. */
  std::vector<Eigen::Vector3i> colors;

};

/** \brief Read a Hdf5 file into an Eigen tensor.
 * \param[in] filepath path to file
 * \param[out] dense Eigen tensor
 * \return success
 */
template<int RANK>
bool read_hdf5(const std::string filepath, Eigen::Tensor<float, RANK, Eigen::RowMajor>& dense) {

  try {
    H5::H5File file(filepath, H5F_ACC_RDONLY);
    H5::DataSet dataset = file.openDataSet("tensor");

    /*
     * Get filespace for rank and dimension
     */
    H5::DataSpace filespace = dataset.getSpace();

    /*
     * Get number of dimensions in the file dataspace
     */
    size_t rank = filespace.getSimpleExtentNdims();

    if (rank != RANK) {
      std::cout << "[Error] invalid rank read: " << rank << std::endl;
      exit(1);
    }

    /*
     * Get and print the dimension sizes of the file dataspace
     */
    hsize_t dimsf[rank];
    filespace.getSimpleExtentDims(dimsf);

    //std::cout << "Read " << filepath << ": ";
    //for (int i = 0; i < RANK; ++i) {
    //  std::cout << dimsf[i] << " ";
    //}
    //std::cout << std::endl;

    /*
     * Define the memory space to read dataset.
     */
    H5::DataSpace mspace(rank, dimsf);

    //size_t buffer_size = 1;
    //for (int i = 0; i < RANK; ++i) {
    //  buffer_size *= dimsf[i];
    //}

    float* buffer = static_cast<float*>(dense.data());
    dataset.read(buffer, H5::PredType::NATIVE_FLOAT, mspace, filespace);

    //for (int i = 0; i < buffer_size; ++i) {
    //  std::cout << buffer[i] << std::endl;
    //}
  }

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
}

/** \brief Write the given set of volumes to h5 file.
 * \param[in] filepath h5 file to write
 * \param[in] dense volume data
 */
template<int RANK>
bool write_hdf5(const std::string filepath, Eigen::Tensor<float, RANK, Eigen::RowMajor>& dense) {

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
    hsize_t dimsf[RANK];
    for (int i = 0; i < RANK; ++i) {
      dimsf[i] = dense.dimension(i);
    }

    H5::DataSpace dataspace(RANK, dimsf);

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

/** \brief Get the rotation matrix by reading the y-angle from angles.
 * \param[in] angles tensor holding the y-angles
 * \param[in] n which angle to take
 * \param[out] R rotation matrix
 */
void get_rotation_matrix(const Eigen::Tensor<float, 2, Eigen::RowMajor> &angles, const int n, Eigen::Matrix3f &R) {
  R = Eigen::Matrix3f::Zero();
  float radians = angles(n, 0);///180.f*M_PI;
  R(0, 0) = std::cos(radians); R(0, 2) = std::sin(radians);
  R(1, 1) = 1;
  R(2, 0) = - std::sin(radians); R(2, 2) = std::cos(radians);
}

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

  /** \brief Test box ray intersection.
   * \param[in] ray ray to intersect with
   * \return intersects
   */
  bool intersect(const Eigen::Vector3f &ray) {
    for (int d = 0; d < 3; d++) {
      assert(this->min(d) < this->max(d));
    }

    // convert to used data structure
    box_ray_intersection::Vector3 _o(0, 0, 0);
    box_ray_intersection::Vector3 _d(ray(0), ray(1), ray(2));
    box_ray_intersection::Ray _ray(_o, _d);
    box_ray_intersection::Vector3 _min(this->min(0), this->min(1), this->min(2));
    box_ray_intersection::Vector3 _max(this->max(0), this->max(1), this->max(2));
    box_ray_intersection::Box _box(_min, _max);

    // interval for valid its is [0,1] as the ray will
    // be the full-length vector to an observed point
    return _box.intersect(_ray, 0, 0.99);
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

/** \brief Voxelization to emulate real KITTI setting.
 * \param[in] depth_maps depth maps to calculate point clouds
 * \param[in] angles angles used for rendering
 * \param[in] max_depth_value maximum depth value in depthmaps
 * \param[in] K camera parameters from rendering
 * \param[in] mesh_center center of mesh used for rendering
 * \param[out] input volumes containing voxelized point clouds
 * \param[out] space volumes containing voxelized space
 */
void voxelize_realistic_UNSCALED(const Eigen::Tensor<float, 3, Eigen::RowMajor> &depth_maps,
    const Eigen::Tensor<float, 2, Eigen::RowMajor> &angles, float max_depth_value,
    const Eigen::Matrix3f& K, const Eigen::Vector3f& mesh_center,
    Eigen::Tensor<float, 4, Eigen::RowMajor>& input, Eigen::Tensor<float, 4, Eigen::RowMajor>& input_sdf,
    Eigen::Tensor<float, 4, Eigen::RowMajor>& space,
    const std::string txt_directory, float padding_factor, float lambda = 1, float ignore = 0) {

  const int N = input.dimension(0);
  const int height = input.dimension(1);
  const int width = input.dimension(2);
  const int depth = input.dimension(3);

  int image_height = depth_maps.dimension(1);
  int image_width = depth_maps.dimension(2);

  Eigen::Matrix3f Kinv = K.inverse();

  for (int n = 0; n < N; n++) {
    std::cout << "[Data] realisticly voxelizing " << n << std::endl;

    Eigen::Matrix3f R;
    get_rotation_matrix(angles, n, R);
    Eigen::Matrix3f Rinv = R.transpose();

    std::vector<Eigen::Vector3f> scaled_points;
    std::vector<Eigen::Vector3f> points;
    std::vector<float> depths;

    // get all points corresponding to the depth map
    for (int h = 0; h < image_height; h++) {
      for (int w = 0; w < image_width; w++) {
        if (depth_maps(n, h, w) < max_depth_value) {
          Eigen::Vector3f pixel;
          // Assume points in the middle of pixels.
          pixel(0) = w + 0.5f;
          pixel(1) = h + 0.5f;
          pixel(2) = 1;

          Eigen::Vector3f point = Kinv*pixel;
          point(2) = depth_maps(n, h, w);

          float r = std::rand() / (float) RAND_MAX;
          if (r < ignore) {
            point(2) = max_depth_value;
          }
          else if (lambda > 0) {
            r = std::rand() / (float) RAND_MAX;
            point(2) -= std::log(r)/lambda;
          }

          point(0) *= point(2); // !
          point(1) *= point(2); // !

          point -= mesh_center;

          // Now point cloud is centered at (0,0,0) and
          // can be rotated.
          point = Rinv*point;

          // Now put point cloud back to mesh center.
          //point -= Rinv*(-mesh_center);

          // Scale point.
          Eigen::Vector3f scaled_point = point;
          scaled_point(0) *= static_cast<float>(width)/padding_factor;
          scaled_point(1) *= static_cast<float>(width)/padding_factor;
          scaled_point(2) *= static_cast<float>(width)/padding_factor;

          // Now put point cloud to voxel grid center.
          scaled_point += Eigen::Vector3f(static_cast<float>(width)/2, static_cast<float>(height)/2, static_cast<float>(depth)/2);

          // Now put point cloud back to mesh center.
          point -= Rinv*(-mesh_center);
          points.push_back(point);
          scaled_points.push_back(scaled_point);
          depths.push_back(depth_maps(n, h, w));

          //point -= mesh_center;
          //point += Eigen::Vector3f(0.5f, 0.5f, 0.5f);

          //point(0) *= width;
          //point(1) *= height;
          //point(2) *= depth;

          //int hh = static_cast<int>(point(1));
          //int ww = static_cast<int>(point(0));
          //int dd = static_cast<int>(point(2));

          //if (hh >= 0 && hh < height && ww >= 0 && ww < width && dd >= 0 && dd < depth) {
          //  input(n, hh, ww, dd) = 1;
          //}
        }
      }
    }

    PointCloud point_cloud(scaled_points);

    std::string txt_file = txt_directory + "/" + std::to_string(n) + ".txt";
    std::cout << "[Data] writing " << txt_file << std::endl;
    bool success = point_cloud.to_txt(txt_file);

    if (!success) {
      std::cout << "[Error] could not write " << txt_file << std::endl;
      exit(1);
    }

    for (int h = 0; h < height; h++) {
      for (int w = 0; w < width; w++) {
        for (int d = 0; d < depth; d++) {
          input_sdf(n, h, w, d) = FLT_MAX;

          // [Data] 99% percentile dimensions:
          // [Data]   4.750000 1.870000 1.970000
          // [Data]   2.540102 1        1.053475
          // [Data]   56       22       24

          float width_factor = (padding_factor*2.545454545)/2.545454545;
          float height_factor = padding_factor/2.545454545;
          float depth_factor = (padding_factor*1.090909091)/2.545454545;

          // the box corresponding to this voxel
          Eigen::Vector3f min(width_factor*static_cast<float>(w)/width, height_factor*static_cast<float>(h)/height,
            depth_factor*static_cast<float>(d)/depth);
          Eigen::Vector3f max(width_factor*static_cast<float>(w + 1)/width, height_factor*static_cast<float>(h + 1)/height,
            depth_factor*static_cast<float>(d + 1)/depth);
          Box voxel_box(min, max);

          voxel_box.translate(Eigen::Vector3f(-width_factor/2, -height_factor/2, -depth_factor/2));
          voxel_box.translate(-Rinv*(-mesh_center));

          Eigen::Vector3f center(w + 0.5f, h + 0.5f, d + 0.5f);

          // no rotation as the intersection does not work
          // and the point cloud is assumed to be rotated

          // go over all "rays"
          for (unsigned int i = 0; i < points.size(); i++) {
            float distance = (center - scaled_points[i]).norm();
            if (distance < input_sdf(n, h, w, d)) {
              input_sdf(n, h, w, d) = distance;
            }

            if (voxel_box.contains(points[i])) {
              input(n, h, w, d) = 1;
            }
            else if (voxel_box.intersect(points[i]) && points[i](2) < max_depth_value) {
              space(n, h, w, d) = 1;
            }
          }
        }
      }
    }
  }
}

/** \brief Voxelization to emulate real KITTI setting.
 * \param[in] depth_maps depth maps to calculate point clouds
 * \param[in] angles angles used for rendering
 * \param[in] max_depth_value maximum depth value in depthmaps
 * \param[in] K camera parameters from rendering
 * \param[in] mesh_center center of mesh used for rendering
 * \param[out] input volumes containing voxelized point clouds
 * \param[out] space volumes containing voxelized space
 */
void voxelize_realistic_SCALED(const Eigen::Tensor<float, 3, Eigen::RowMajor> &depth_maps,
    const Eigen::Tensor<float, 2, Eigen::RowMajor> &angles, float max_depth_value,
    const Eigen::Matrix3f& K, const Eigen::Vector3f& mesh_center,
    Eigen::Tensor<float, 4, Eigen::RowMajor>& input, Eigen::Tensor<float, 4, Eigen::RowMajor>& input_sdf,
    Eigen::Tensor<float, 4, Eigen::RowMajor>& space,
    const std::string txt_directory, float padding_factor, float lambda = 1, float ignore = 0) {

  const int N = input.dimension(0);
  const int height = input.dimension(1);
  const int width = input.dimension(2);
  const int depth = input.dimension(3);

  int image_height = depth_maps.dimension(1);
  int image_width = depth_maps.dimension(2);

  Eigen::Matrix3f Kinv = K.inverse();

  for (int n = 0; n < N; n++) {
    std::cout << "[Data] realisticly voxelizing " << n << std::endl;

    Eigen::Matrix3f R;
    get_rotation_matrix(angles, n, R);
    Eigen::Matrix3f Rinv = R.transpose();

    std::vector<Eigen::Vector3f> points;
    std::vector<float> depths;

    // get all points corresponding to the depth map
    for (int h = 0; h < image_height; h++) {
      for (int w = 0; w < image_width; w++) {
        if (depth_maps(n, h, w) < max_depth_value) {
          Eigen::Vector3f pixel;
          // Assume points in the middle of pixels.
          pixel(0) = w + 0.5f;
          pixel(1) = h + 0.5f;
          pixel(2) = 1;

          Eigen::Vector3f point = Kinv*pixel;
          point(2) = depth_maps(n, h, w);

          float r = std::rand() / (float) RAND_MAX;
          if (r < ignore) {
            point(2) = max_depth_value;
          }
          else if (lambda > 0) {
            r = std::rand() / (float) RAND_MAX;
            point(2) -= std::log(r)/lambda;
          }

          point(0) *= point(2); // !
          point(1) *= point(2); // !

          point -= mesh_center;

          // Now point cloud is centered at (0,0,0) and
          // can be rotated.
          point = Rinv*point;

          // Now put point cloud back to mesh center.
          //point -= Rinv*(-mesh_center);

          // Scale point.
          point(0) *= static_cast<float>(width)/padding_factor;
          point(1) *= static_cast<float>(width)/padding_factor;
          point(2) *= static_cast<float>(width)/padding_factor;

          // Now put point cloud to voxel grid center.
          point += Eigen::Vector3f(static_cast<float>(width)/2, static_cast<float>(height)/2, static_cast<float>(depth)/2);

          points.push_back(point);
          depths.push_back(depth_maps(n, h, w));

          //point -= mesh_center;
          //point += Eigen::Vector3f(0.5f, 0.5f, 0.5f);

          //point(0) *= width;
          //point(1) *= height;
          //point(2) *= depth;

          //int hh = static_cast<int>(point(1));
          //int ww = static_cast<int>(point(0));
          //int dd = static_cast<int>(point(2));

          //if (hh >= 0 && hh < height && ww >= 0 && ww < width && dd >= 0 && dd < depth) {
          //  input(n, hh, ww, dd) = 1;
          //}
        }
      }
    }

    PointCloud point_cloud(points);

    std::string txt_file = txt_directory + "/" + std::to_string(n) + ".txt";
    std::cout << "[Data] writing " << txt_file << std::endl;
    bool success = point_cloud.to_txt(txt_file);

    if (!success) {
      std::cout << "[Error] could not write " << txt_file << std::endl;
      exit(1);
    }

    for (int h = 0; h < height; h++) {
      for (int w = 0; w < width; w++) {
        for (int d = 0; d < depth; d++) {
          input_sdf(n, h, w, d) = FLT_MAX;

          // the box corresponding to this voxel
          Eigen::Vector3f min(static_cast<float>(w), static_cast<float>(h), static_cast<float>(d));
          Eigen::Vector3f max(static_cast<float>(w + 1), static_cast<float>(h + 1), static_cast<float>(d + 1));
          Eigen::Vector3f center(min(0) + (max(0) - min(0))/2, min(1) + (max(1) - min(1))/2, min(2) + (max(2) - min(2))/2);
          Box voxel_box(min, max);

          //voxel_box.translate(Eigen::Vector3f(-static_cast<float>(width)/2, -static_cast<float>(height)/2, -static_cast<float>(depth)/2));
          //voxel_box.translate(-Rinv*(-mesh_center));

          // no rotation as the intersection does not work
          // and the point cloud is assumed to be rotated

          // go over all "rays"
          for (unsigned int i = 0; i < points.size(); i++) {
            float distance = (center - points[i]).norm();
            if (distance < input_sdf(n, h, w, d)) {
              input_sdf(n, h, w, d) = distance;
            }

            if (voxel_box.contains(points[i])) {
              input(n, h, w, d) = 1;
            }
            else if (voxel_box.intersect(points[i]) && depths[i] < max_depth_value) {
              space(n, h, w, d) = 1;
            }
          }
        }
      }
    }
  }
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

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cout << "[Error] Usage: voxelize_cuboids config.json" << std::endl;
    exit(1);
  }

  Json::Value root;
  Json::Reader reader;

  std::ifstream json_file(argv[1], std::ifstream::in | std::ifstream::binary);
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

  std::string depth_file = json.get_string("depth_file") + "_" + std::to_string(multiplier) + "_"
      + std::to_string(image_height) + "x" + std::to_string(image_width) + "_"
      + std::to_string(height) + "x" + std::to_string(width) + "x" + std::to_string(depth) + suffix + ".h5";

  if (!boost::filesystem::exists(boost::filesystem::path(depth_file))) {
    std::cout << "[Error] depth file " << depth_file << " not found" << std::endl;
    exit(1);
  }

  std::string angles_file = json.get_string("render_orientation_file") + "_" + std::to_string(multiplier) + "_"
      + std::to_string(image_height) + "x" + std::to_string(image_width) + "_"
      + std::to_string(height) + "x" + std::to_string(width) + "x" + std::to_string(depth) + suffix + ".h5";

  if (!boost::filesystem::exists(boost::filesystem::path(angles_file))) {
    std::cout << "[Error] angles file " << angles_file << " not found" << std::endl;
    exit(1);
  }

  std::string off_directory = json.get_string("off_dir");
  if (off_directory.empty()) {
    std::cout << "[Error] Read invalid off directory" << std::endl;
    exit(1);
  }

  off_directory = off_directory + "_" + std::to_string(multiplier) + "_"
      + std::to_string(image_height) + "x" + std::to_string(image_width) + "_"
      + std::to_string(height) + "x" + std::to_string(width) + "x" + std::to_string(depth) + suffix;

  std::string txt_directory = json.get_string("txt_gt_dir");
  if (txt_directory.empty()) {
    std::cout << "[Error] Read invalid ply gt directory" << std::endl;
    exit(1);
  }

  txt_directory = txt_directory + "_" + std::to_string(multiplier) + "_"
      + std::to_string(image_height) + "x" + std::to_string(image_width) + "_"
      + std::to_string(height) + "x" + std::to_string(width) + "x" + std::to_string(depth) + suffix;

  if (!boost::filesystem::is_directory(boost::filesystem::path(txt_directory))) {
    boost::filesystem::create_directories(boost::filesystem::path(txt_directory));
  }

  std::map<int, std::string> files;
  read_directory(off_directory, files);

  int N = files.size();
  Eigen::Tensor<float, 3, Eigen::RowMajor> depth_maps(N, image_height, image_width);
  depth_maps.setZero();
  read_hdf5<3>(depth_file, depth_maps);

  Eigen::Tensor<float, 2, Eigen::RowMajor> angles(N, 1);
  angles.setZero();
  read_hdf5<2>(angles_file, angles);

  float max_depth_value = json.get_float("max_depth_value");
  Eigen::Matrix3f K = Eigen::Matrix3f::Zero();
  K(0, 0) = json.get_float("focal_length_x");
  K(0, 2) = json.get_float("principal_point_x");
  K(1, 1) = json.get_float("focal_length_y");
  K(1, 2) = json.get_float("principal_point_y");
  K(2, 2) = 1;

  Eigen::Vector3f mesh_center;
  mesh_center(0) = json.get_float("mesh_center_x");
  mesh_center(1) = json.get_float("mesh_center_y");
  mesh_center(2) = json.get_float("mesh_center_z");

  Eigen::Tensor<float, 4, Eigen::RowMajor> real_input(N, height, width, depth);
  real_input.setZero();
  Eigen::Tensor<float, 4, Eigen::RowMajor> real_input_sdf(N, height, width, depth);
  real_input_sdf.setZero();
  Eigen::Tensor<float, 4, Eigen::RowMajor> real_space(N, height, width, depth);
  real_space.setZero();

  float padding_factor = 1 + json.get_float("padding");
  float lambda = json.get_float("lambda");
  float ignore = json.get_float("ignore");
  voxelize_realistic_UNSCALED(depth_maps, angles, max_depth_value, K, mesh_center, real_input, real_input_sdf, real_space, txt_directory,
    padding_factor, lambda, ignore);

  std::string real_input_file = json.get_string("input_file") + "_" + std::to_string(multiplier) + "_"
      + std::to_string(image_height) + "x" + std::to_string(image_width) + "_"
      + std::to_string(height) + "x" + std::to_string(width) + "x" + std::to_string(depth) + suffix + ".h5";
  bool success = write_hdf5(real_input_file, real_input);

  if (success) {
    std::cout <<"[Data] wrote " << real_input_file << std::endl;
  }
  else {
    std::cout << "[Error] error writing " << real_input_file << std::endl;
  }

  std::string real_input_sdf_file = json.get_string("input_sdf_file") + "_" + std::to_string(multiplier) + "_"
      + std::to_string(image_height) + "x" + std::to_string(image_width) + "_"
      + std::to_string(height) + "x" + std::to_string(width) + "x" + std::to_string(depth) + suffix + ".h5";
  success = write_hdf5(real_input_sdf_file, real_input_sdf);

  if (success) {
    std::cout <<"[Data] wrote " << real_input_sdf_file << std::endl;
  }
  else {
    std::cout << "[Error] error writing " << real_input_sdf_file << std::endl;
  }

  std::string real_space_file = json.get_string("space_file") + "_" + std::to_string(multiplier) + "_"
      + std::to_string(image_height) + "x" + std::to_string(image_width) + "_"
      + std::to_string(height) + "x" + std::to_string(width) + "x" + std::to_string(depth) + suffix + ".h5";
  success = write_hdf5(real_space_file, real_space);

  if (success) {
    std::cout <<"[Data] wrote " << real_space_file << std::endl;
  }
  else {
    std::cout << "[Error] error writing " << real_space_file << std::endl;
  }

  exit(0);
}
