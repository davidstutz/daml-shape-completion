#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cstdarg>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <boost/filesystem.hpp>
#include <json/json.h>
#include <H5Cpp.h>

#include "box/vector3.h"
#include "box/ray.h"
#include "box/box.h"
#include "icp/icpPointToPoint.h"

/** \brief Format a string as done by sprintf.
 * \see https://stackoverflow.com/questions/2342162/stdstring-formatting-like-sprintf
 */
std::string format_string(const std::string fmt_str, ...) {
  int final_n, n = ((int)fmt_str.size()) * 2; /* Reserve two times as much as the length of the fmt_str */
  std::string str;
  std::unique_ptr<char[]> formatted;
  va_list ap;

  while(1) {
    formatted.reset(new char[n]); /* Wrap the plain char array into the unique_ptr */
    strcpy(&formatted[0], fmt_str.c_str());
    va_start(ap, fmt_str);
    final_n = vsnprintf(&formatted[0], n, fmt_str.c_str(), ap);
    va_end(ap);
    if (final_n < 0 || final_n >= n) {
      n += abs(final_n - n + 1);
    }
    else {
      break;
    }
  }

  return std::string(formatted.get());
}

/** \brief Format a vector to a string.
 * \param[in] vector vector to format
 * \return string
 */
template<typename T>
std::string format_vector(Eigen::Matrix<T, 3, 1> vector) {
  std::string str;
  for (int d = 0; d < 3; d++) {
    str += std::to_string(vector(d)) + " ";
  }

  return str;
}

/** \brief Read all files in a directory matching the given extension.
 * \param[in] directory path to directory
 * \param[out] files read file paths
 * \param[in] extension extension to filter for
 */
void read_directory(const boost::filesystem::path directory, std::map<int, boost::filesystem::path>& files, const std::string extension = ".txt") {
  files.clear();
  boost::filesystem::directory_iterator end;

  for (boost::filesystem::directory_iterator it(directory); it != end; ++it) {
    if (it->path().extension().string() == extension) {
      if (!boost::filesystem::is_empty(it->path()) && !it->path().empty() && it->path().filename().string() != "") {
        int number = std::stoi(it->path().filename().string());
        files.insert(std::pair<int, boost::filesystem::path>(number, it->path()));
      }
    }
  }
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
    return _box.intersect(_ray, 0, 1);
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

/** \brief Read the number of bounding boxes, i.e. a single integer from the given file.
 * \param[in] filepath path to file
 * \return number of bounding boxes
 */
int read_num_bounding_boxes(const std::string filepath) {
  std::ifstream file(filepath.c_str());
  std::string line;
  std::stringstream ss;

  std::getline(file, line);
  ss << line;

  int n_bounding_boxes;
  ss >> n_bounding_boxes;

  return n_bounding_boxes;
}

/** \brief Read file containing split indices.
 * \param[in] filepath path to file
 * \param[out] indices read indices
 */
bool read_split_file(const std::string filepath, std::vector<int> &indices) {
  std::ifstream file(filepath.c_str());

  std::string line;
  while (std::getline(file, line))
  {
      std::istringstream iss(line);
      int index = -1;
      iss >> index;
      indices.push_back(index);
  }

  return true;
}

/** \brief Read the given TXT file and extract the corresponding bounding boxes.
 * \param[in] filepath path to file to read
 * \param[out] bounding_boxes read bounding boxes
 * \return success
 */
bool read_bounding_boxes(const std::string filepath, std::vector<BoundingBox> &bounding_boxes) {
  std::ifstream file(filepath.c_str());
  std::string line;
  std::stringstream ss;

  std::getline(file, line);
  ss << line;

  int n_bounding_boxes;
  ss >> n_bounding_boxes;

  if (n_bounding_boxes < 0) {
    return false;
  }

  for (int i = 0; i < n_bounding_boxes; i++) {
    std::getline(file, line);

    ss.clear();
    ss.str("");
    ss << line;

    BoundingBox bounding_box;

    ss >> bounding_box.size(0);
    ss >> bounding_box.size(1);
    ss >> bounding_box.size(2);

    ss >> bounding_box.translation(0);
    ss >> bounding_box.translation(1);
    ss >> bounding_box.translation(2);

    ss >> bounding_box.rotation(0);
    ss >> bounding_box.rotation(1);
    ss >> bounding_box.rotation(2);

    bounding_boxes.push_back(bounding_box);
  }

  return true;
}

/** \brief Write bounding boxes.
 * \param[in] filepath path to file to write to
 * \param[in] bounding_boxes bounding boxes
 * \return success
 */
bool write_bounding_boxes(const std::string filepath, const std::vector<BoundingBox> &bounding_boxes) {
  std::ofstream* file = new std::ofstream(filepath.c_str());

  if (!static_cast<bool>(*file)) {
    return false;
  }

  (*file) << bounding_boxes.size() << std::endl;

  for (unsigned int i = 0; i < bounding_boxes.size(); i++) {
    BoundingBox bounding_box = bounding_boxes[i];

    (*file) << bounding_box.size(0) << " ";
    (*file) << bounding_box.size(1) << " ";
    (*file) << bounding_box.size(2) << " ";

    (*file) << bounding_box.translation(0) << " ";
    (*file) << bounding_box.translation(1) << " ";
    (*file) << bounding_box.translation(2) << " ";

    (*file) << bounding_box.rotation(0) << " ";
    (*file) << bounding_box.rotation(1) << " ";
    (*file) << bounding_box.rotation(2) << " ";

    (*file) << bounding_box.meta << std::endl;
  }

  file->close();
  delete file;

  return true;
}

/** \brief Class representing a point cloud in 3D. */
class PointCloud {
public:
  /** \brief Constructor. */
  PointCloud() {

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

    // allocate 4 MB buffer (only ~130*4*4 KB are needed)
    int32_t num = 1000000;
    float *data = (float*)malloc(num*sizeof(float));

    // pointers
    float *px = data+0;
    float *py = data+1;
    float *pz = data+2;
    float *pr = data+3;

    // load point cloud
    FILE *stream;
    stream = fopen (filepath.c_str(), "rb");
    num = fread(data,sizeof(float),num,stream)/4;
    for (int32_t i=0; i<num; i++) {
      point_cloud.add_point(Eigen::Vector3f(-*py,*pz,*px));
      //std::cout << *px << " " << *px << " " << *pz << std::endl;
      px+=4; py+=4; pz+=4; pr+=4;
    }

    fclose(stream);

    //exit(1);
    return true;
  }

  /** \brief Given the angle in radians, construct a rotation matrix around the x-axis.
   * \param[in] radians angle in radians
   * \param[out] rotation rotation matrix
   */
  static void rotation_matrix_x(const float radians, Eigen::Matrix3f &rotation) {
    rotation = Eigen::Matrix3f::Zero();

    rotation(0, 0) = 1;
    rotation(1, 1) = std::cos(radians); rotation(1, 2) = -std::sin(radians);
    rotation(2, 1) = std::sin(radians); rotation(2, 2) = std::cos(radians);
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

  /** \brief Given the angle in radians, construct a rotation matrix around the z-axis.
   * \param[in] radians angle in radians
   * \param[out] rotation rotation matrix
   */
  static void rotation_matrix_z(const float radians, Eigen::Matrix3f &rotation) {
    rotation = Eigen::Matrix3f::Zero();

    rotation(0, 0) = std::cos(radians); rotation(0, 1) = -std::sin(radians);
    rotation(1, 0) = std::sin(radians); rotation(1, 1) = std::cos(radians);
    rotation(2, 2) = 1;
  }

  /** \brief Computes the rotation matrix corresponding to the given ray.
   * \param[in] ray ray defining the direction to rotate to
   * \param[out] rotation final rotation
   */
  static void rotation_matrix(const Eigen::Vector3f ray, Eigen::Matrix3f &rotation) {
    Eigen::Matrix3f rotation_x;
    Eigen::Matrix3f rotation_y;
    Eigen::Matrix3f rotation_z;

    Eigen::Vector3f axis_x = Eigen::Vector3f(1, 0, 0);
    Eigen::Vector3f axis_y = Eigen::Vector3f(0, 1, 0);
    Eigen::Vector3f axis_z = Eigen::Vector3f(0, 0, 1);

    Eigen::Vector3f ray_x = ray; ray_x(0) = 0; ray_x /= ray_x.norm();
    Eigen::Vector3f ray_y = ray; ray_y(1) = 0; ray_y /= ray_y.norm();
    Eigen::Vector3f ray_z = ray; ray_z(2) = 0; ray_y /= ray_z.norm();

    float radians_x = std::acos(axis_x.dot(ray_x));
    PointCloud::rotation_matrix_x(radians_x, rotation_x);

    float radians_y = std::acos(axis_y.dot(ray_y));
    PointCloud::rotation_matrix_y(radians_y, rotation_y);

    float radians_z = std::acos(axis_z.dot(ray_z));
    PointCloud::rotation_matrix_z(radians_z, rotation_z);
    std::cout << "[Data] radians " << radians_x << " " << radians_y << " " << radians_z << std::endl;

    rotation = Eigen::Matrix3f::Zero();
    rotation = rotation_z*rotation_y*rotation_x;
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

  /** \brief Add points from a point cloud.
   * \param[in] point_cloud point cloud whose points to add
   */
  void add_points(const PointCloud &point_cloud) {
    for (unsigned int i = 0; i < point_cloud.num_points(); i++) {
      this->add_point(point_cloud.points[i], point_cloud.colors[i]);
    }
  }

  /** \brief Merge/add points from another point cloud.
   * \param[in] point_cloud point_cloud to take points from
   */
  void merge(const PointCloud &point_cloud) {
    for (unsigned int i = 0; i < point_cloud.points.size(); i++) {
      this->add_point(point_cloud.points[i], point_cloud.colors[i]);
    }
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
      //(*out) << this->points[i](0) << " " << this->points[i](1) << " " << this->points[i](2) << std::endl;
      (*out) << this->points[i](0) << " " << this->points[i](1) << " " << this->points[i](2) << " "
        << this->colors[i](0) << " " << this->colors[i](1) << " " << this->colors[i](2) << std::endl;
    }

    out->close();
    delete out;

    return true;
  }

  /** \brief Write point cloud to txt file.
   * \param[in] bounding_box bounding box to write
   * \param[in] filepath path to file
   * \return success
   */
  void extract(const BoundingBox bounding_box, PointCloud &point_cloud) {
    for (unsigned int i = 0; i < this->points.size(); i++) {
      if (bounding_box.contains(this->points[i])) {
        point_cloud.add_point(this->points[i]);
      }
    }
  }

  /** \brief Get the extents of the point cloud in all three axes.
   * \param[out] min minimum coordinate values per axis
   * \param[out] max maximum coordinate values per axis
   */
  void extents(Eigen::Vector3f &min, Eigen::Vector3f &max) {
    for (int d = 0; d < 3; d++) {
      min(d) = FLT_MAX;
      max(d) = FLT_MIN;
    }

    for (unsigned int i = 0; i < this->points.size(); i++) {
      for (int d = 0; d < 3; d++) {
        if (this->points[i](d) < min(d)) {
          min(d) = this->points[i](d);
        }
        if (this->points[i](d) > max(d)) {
          max(d) = this->points[i](d);
        }
      }
    }
  }

  /** \brief Mark the origin and its axes.
   */
  void mark() {
      this->add_point(Eigen::Vector3f(0, 0, 0), Eigen::Vector3i(255, 0, 0));
      this->add_point(Eigen::Vector3f(1, 0, 0), Eigen::Vector3i(255, 0, 0));
      this->add_point(Eigen::Vector3f(0, 1, 0), Eigen::Vector3i(0, 255, 0));
      this->add_point(Eigen::Vector3f(0, 0, 1), Eigen::Vector3i(0, 0, 255));
  }

  /** \brief Mark the corners of the bounding box in the point cloud.
   * \param[in] bounding_box bounding box to mark
   */
  void mark(const BoundingBox &bounding_box) {
    Eigen::Matrix3f rotation;
    PointCloud::rotation_matrix_y(bounding_box.rotation(1), rotation);
    rotation.transposeInPlace();

    for (unsigned int i = 0; i < this->points.size(); i++) {
      if (bounding_box.contains(this->points[i])) {
        this->colors[i] = Eigen::Vector3i(0, 0, 255);
      }
    }

    {
      Eigen::Vector3f point(-bounding_box.size(0)/2, -bounding_box.size(1)/2, -bounding_box.size(2)/2);
      point = rotation*point + bounding_box.translation;
      this->add_point(point, Eigen::Vector3i(0, 255, 0));
    }
    {
      Eigen::Vector3f point(-bounding_box.size(0)/2, -bounding_box.size(1)/2, bounding_box.size(2)/2);
      point = rotation*point + bounding_box.translation;
      this->add_point(point, Eigen::Vector3i(0, 255, 0));
    }
    {
      Eigen::Vector3f point(-bounding_box.size(0)/2, bounding_box.size(1)/2, -bounding_box.size(2)/2);
      point = rotation*point + bounding_box.translation;
      this->add_point(point, Eigen::Vector3i(0, 255, 0));
    }
    {
      Eigen::Vector3f point(-bounding_box.size(0)/2, bounding_box.size(1)/2, bounding_box.size(2)/2);
      point = rotation*point + bounding_box.translation;
      this->add_point(point, Eigen::Vector3i(0, 255, 0));
    }
    {
      Eigen::Vector3f point(bounding_box.size(0)/2, -bounding_box.size(1)/2, -bounding_box.size(2)/2);
      point = rotation*point + bounding_box.translation;
      this->add_point(point, Eigen::Vector3i(0, 255, 0));
    }
    {
      Eigen::Vector3f point(bounding_box.size(0)/2, -bounding_box.size(1)/2, bounding_box.size(2)/2);
      point = rotation*point + bounding_box.translation;
      this->add_point(point, Eigen::Vector3i(0, 255, 0));
    }
    {
      Eigen::Vector3f point(bounding_box.size(0)/2, bounding_box.size(1)/2, -bounding_box.size(2)/2);
      point = rotation*point + bounding_box.translation;
      this->add_point(point, Eigen::Vector3i(0, 255, 0));
    }
    {
      Eigen::Vector3f point(bounding_box.size(0)/2, bounding_box.size(1)/2, bounding_box.size(2)/2);
      point = rotation*point + bounding_box.translation;
      this->add_point(point, Eigen::Vector3i(0, 255, 0));
    }
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

  /** \brief Voxelize the point cloud into the given dense array.
   * Assumes the point cloud to be rotated such that the bounding
   * box is axis-aligned.
   * \param[in] bounding_box bounding box to voxelize
   * \param[in] dense dense tensor to voxelize to
   * \param[in] n index of first dimension in dense
   */
  void voxelize_DEPRECATED(const BoundingBox &bounding_box, Eigen::Tensor<float, 4, Eigen::RowMajor> &dense, const int n) {
    const int height = dense.dimension(1);
    const int width = dense.dimension(2);
    const int depth = dense.dimension(3);

    Eigen::Matrix3f rotation;
    BoundingBox::rotation_matrix_y(bounding_box.rotation(1), rotation);

    for (unsigned int i = 0; i < this->points.size(); i++) {
      //Eigen::Vector3f point = this->points[i] - bounding_box.translation;
      //point = rotation*point;

      // After rotating, axis 1 will always be the longer axis
      // of the bounding box.
      //point /= bounding_box.size(0);
      //point += Eigen::Vector3f(0.5, bounding_box.size(1)/bounding_box.size(0)/2, bounding_box.size(2)/bounding_box.size(0)/2);

      //if (point(0) <= 1 && point(0) >= 0
      //    && point(1) <= bounding_box.size(1)/bounding_box.size(0) && point(1) >= 0
      //    && point(2) <= bounding_box.size(2)/bounding_box.size(0) && point(2) >= 0) {

      //  int h = static_cast<int>(point(1)*height); // y
      //  int w = static_cast<int>(point(0)*width); // x
      //  int d = static_cast<int>(point(2)*depth); // z

        // Mark!
      //  this->colors[i] = Eigen::Vector3i(255, 0, 0);

      //  if (h >= 0 && h < height && w >= 0 && w < width && d >= 0 && d < depth) {
      //    dense(n, h, w, d) = 1;
      //  }
      //}
      if (bounding_box.contains(this->points[i])) {

        Eigen::Vector3f point = this->points[i] - bounding_box.translation;
        point = rotation*point;

        point /= bounding_box.size(0);
        point += Eigen::Vector3f(0.5, bounding_box.size(1)/bounding_box.size(0)/2, bounding_box.size(2)/bounding_box.size(0)/2);

        int h = static_cast<int>(point(1)*height); // y
        int w = static_cast<int>(point(0)*width); // x
        int d = static_cast<int>(point(2)*depth); // z

        // Mark!
        this->colors[i] = Eigen::Vector3i(255, 0, 0);

        if (h >= 0 && h < height && w >= 0 && w < width && d >= 0 && d < depth) {
          dense(n, h, w, d) = 1;
        }
      }
    }
  }

  /** \brief Voxelizes free space; needs the corresponding depth image and
   * assumes the point cloud to be aligned (i.e. rotated) with the bounding box sides.
   * \param[in] camera camera used for rendering
   * \param[in] input tensor to hold voxelized point clouds
   * \param[in] full_space tensor to hold voxelized free space
   * \param[in] part_space tensor to hold voxelized (partial) free space
   * \param[in] index of first dimension in dense and depth
   */
  void voxelize_UNSCALED(const BoundingBox &bounding_box,
      Eigen::Tensor<float, 4, Eigen::RowMajor> &input,
      Eigen::Tensor<float, 4, Eigen::RowMajor> &full_space,
      Eigen::Tensor<float, 4, Eigen::RowMajor> &part_space, const int n) {

    const int height = input.dimension(1);
    const int width = input.dimension(2);
    const int depth = input.dimension(3);

    Box box;
    bounding_box.to_box(box);
    std::vector<int> considered_rays;
    std::vector<int> considered_points;

    for (unsigned int i = 0; i < this->points.size(); i++) {
      if (box.intersect(this->points[i])) {
        considered_rays.push_back(i);
      }
      if (box.contains(this->points[i])) {
        considered_points.push_back(i);
        this->colors[i] = Eigen::Vector3i(0, 255, 0);
      }
    }

    for (int h = 0; h < height; h++) {
      for (int w = 0; w < width; w++) {
        for (int d = 0; d < depth; d++) {

          // [Data] average dimensions:
          // [Data]   3.892065 1.618767 1.529863
          // [Data] max dimensions:
          // [Data]   6.670000 1.990000 2.480000
          // [Data]   3.375    1        1.25
          // [Data]   54       16       20

          float width_factor = 3.375;
          float height_factor = 1;
          float depth_factor = 1.25;

          // [Data] 99% percentile dimensions:
          // [Data]   4.750000 1.870000 1.970000
          // [Data]   2.540102 1        1.053475
          // [Data]   56       22       24

          width_factor = 2.545454545;
          height_factor = 1;
          depth_factor = 1.090909091;

          // the box corresponding to this voxel
          Eigen::Vector3f min(width_factor*static_cast<float>(w)/width, height_factor*static_cast<float>(h)/height,
            depth_factor*static_cast<float>(d)/depth);
          Eigen::Vector3f max(width_factor*static_cast<float>(w + 1)/width, height_factor*static_cast<float>(h + 1)/height,
            depth_factor*static_cast<float>(d + 1)/depth);
          Box voxel_box(min, max);

          voxel_box.translate(Eigen::Vector3f(-width_factor/2, -height_factor/2, -depth_factor/2));
          voxel_box.scale(bounding_box.size(1));
          voxel_box.translate(bounding_box.translation);

          // no rotation as the intersection does not work
          // and the point cloud is assumed to be rotated

          // go over all "rays"
          for (unsigned int i = 0; i < considered_rays.size(); i++) {
            if (voxel_box.intersect(this->points[considered_rays[i]])) {
              full_space(n, h, w, d) = 1;
            }
          }

          for (unsigned int i = 0; i < considered_points.size(); i++) {
            if (voxel_box.contains(this->points[considered_points[i]])) {
              //this->colors[considered_points[i]] = Eigen::Vector3i(0, 0, 255);
              input(n, h, w, d) = 1;
            }
            if (voxel_box.intersect(this->points[considered_points[i]])) {
              part_space(n, h, w, d) = 1;
            }
          }
        }
      }
    }
  }

  /** \brief Voxelizes free space; needs the corresponding depth image and
   * assumes the point cloud to be aligned (i.e. rotated) with the bounding box sides.
   * \param[in] camera camera used for rendering
   * \param[in] input tensor to hold voxelized point clouds
   * \param[in] full_space tensor to hold voxelized free space
   * \param[in] part_space tensor to hold voxelized (partial) free space
   * \param[in] index of first dimension in dense and depth
   */
  void voxelize_SCALED(const BoundingBox &bounding_box,
    const BoundingBox &bounding_box_small,
      Eigen::Tensor<float, 4, Eigen::RowMajor> &input,
      Eigen::Tensor<float, 4, Eigen::RowMajor> &input_sdf,
      Eigen::Tensor<float, 4, Eigen::RowMajor> &full_space,
      Eigen::Tensor<float, 4, Eigen::RowMajor> &part_space, const int n) {

    const int height = input.dimension(1);
    const int width = input.dimension(2);
    const int depth = input.dimension(3);

    Box box;
    bounding_box_small.to_box(box);
    std::vector<int> considered_rays;
    std::vector<int> considered_points;

    for (unsigned int i = 0; i < this->points.size(); i++) {
      if (box.intersect(this->points[i])) {
        considered_rays.push_back(i);
      }
      if (box.contains(this->points[i])) {
        considered_points.push_back(i);
        this->colors[i] = Eigen::Vector3i(0, 255, 0);
      }
    }

    //std::cout << "[Data] SCALED 4" << std::endl;
    //std::cout << "[Data] box " << box.min(0) << " - " << box.max(0)
    //  << " " << box.min(1) << " - " << box.max(1) << " " << box.min(2) << " - " << box.max(2) << std::endl;
    //std::cout << "[Data] voxel grid "
    //  << -static_cast<float>(width)/2 + bounding_box.translation(0) << " - " << static_cast<float>(width)/2 + bounding_box.translation(0) << " "
    //  << -static_cast<float>(height)/2 + bounding_box.translation(1) << " - " << static_cast<float>(height)/2 + bounding_box.translation(1) << " "
    //  << -static_cast<float>(depth)/2 + bounding_box.translation(2) << " - " << static_cast<float>(depth)/2 + bounding_box.translation(2)
    //  << std::endl;

    int input_count = 0;
    for (int h = 0; h < height; h++) {
      for (int w = 0; w < width; w++) {
        for (int d = 0; d < depth; d++) {
          input_sdf(n, h, w, d) = FLT_MAX;

          // the box corresponding to this voxel
          Eigen::Vector3f min(static_cast<float>(w), static_cast<float>(h), static_cast<float>(d));
          Eigen::Vector3f max(static_cast<float>(w + 1), static_cast<float>(h + 1), static_cast<float>(d + 1));
          Eigen::Vector3f center(w + 0.5f, h + 0.5f, d + 0.5f);
          Box voxel_box(min, max);

          center += Eigen::Vector3f(-static_cast<float>(width)/2, -static_cast<float>(height)/2, -static_cast<float>(depth)/2);
          voxel_box.translate(Eigen::Vector3f(-static_cast<float>(width)/2, -static_cast<float>(height)/2, -static_cast<float>(depth)/2));
          center += bounding_box.translation;
          //voxel_box.scale(bounding_box.size(1));
          voxel_box.translate(bounding_box.translation);

          // no rotation as the intersection does not work
          // and the point cloud is assumed to be rotated

          // go over all "rays"
          for (unsigned int i = 0; i < considered_rays.size(); i++) {
            if (voxel_box.intersect(this->points[considered_rays[i]])) {
              full_space(n, h, w, d) = 1;
            }
          }

          for (unsigned int i = 0; i < considered_points.size(); i++) {
            float distance = (center - this->points[considered_points[i]]).norm();
            //std::cout << "distance: " << distance << std::endl;

            if (distance < input_sdf(n, h, w, d)) {
              input_sdf(n, h, w, d) = distance;
            }

            if (voxel_box.contains(this->points[considered_points[i]])) {
              //this->colors[considered_points[i]] = Eigen::Vector3i(0, 0, 255);
              input(n, h, w, d) = 1;
              input_count++;
            }
            if (voxel_box.intersect(this->points[considered_points[i]])) {
              part_space(n, h, w, d) = 1;
            }
          }
        }
      }
    }

    //std::cout << "[Data] " << input_count << " observed voxels" << std::endl;
  }

  /** \brief Voxelizes free space; needs the corresponding depth image and
   * assumes the point cloud to be aligned (i.e. rotated) with the bounding box sides.
   * \param[in] camera camera used for rendering
   * \param[in] input tensor to hold voxelized point clouds
   * \param[in] full_space tensor to hold voxelized free space
   * \param[in] part_space tensor to hold voxelized (partial) free space
   * \param[in] index of first dimension in dense and depth
   */
  void voxelize_SCALED(const BoundingBox &bounding_box,
    const BoundingBox &bounding_box_small,
      Eigen::Tensor<float, 5, Eigen::RowMajor> &input,
      Eigen::Tensor<float, 5, Eigen::RowMajor> &full_space,
      Eigen::Tensor<float, 5, Eigen::RowMajor> &part_space,
      const int n, const int m) {

    const int height = input.dimension(2);
    const int width = input.dimension(3);
    const int depth = input.dimension(4);

    Box box;
    bounding_box_small.to_box(box);
    std::vector<int> considered_rays;
    std::vector<int> considered_points;

    for (unsigned int i = 0; i < this->points.size(); i++) {
      if (box.intersect(this->points[i])) {
        considered_rays.push_back(i);
      }
      if (box.contains(this->points[i])) {
        considered_points.push_back(i);
        this->colors[i] = Eigen::Vector3i(0, 255, 0);
      }
    }

    //std::cout << "[Data] SCALED 5" << std::endl;
    //std::cout << "[Data] box " << box.min(0) << " - " << box.max(0)
    //  << " " << box.min(1) << " - " << box.max(1) << " " << box.min(2) << " - " << box.max(2) << std::endl;
    //std::cout << "[Data] voxel grid "
    //  << -static_cast<float>(width)/2 + bounding_box.translation(0) << " - " << static_cast<float>(width)/2 + bounding_box.translation(0) << " "
    //  << -static_cast<float>(height)/2 + bounding_box.translation(1) << " - " << static_cast<float>(height)/2 + bounding_box.translation(1) << " "
    //  << -static_cast<float>(depth)/2 + bounding_box.translation(2) << " - " << static_cast<float>(depth)/2 + bounding_box.translation(2)
    //  << std::endl;

    int input_count = 0;
    for (int h = 0; h < height; h++) {
      for (int w = 0; w < width; w++) {
        for (int d = 0; d < depth; d++) {
          // the box corresponding to this voxel
          Eigen::Vector3f min(static_cast<float>(w), static_cast<float>(h), static_cast<float>(d));
          Eigen::Vector3f max(static_cast<float>(w + 1), static_cast<float>(h + 1), static_cast<float>(d + 1));
          Eigen::Vector3f center(w + 0.5f, h + 0.5f, d + 0.5f);
          Box voxel_box(min, max);

          center += Eigen::Vector3f(-static_cast<float>(width)/2, -static_cast<float>(height)/2, -static_cast<float>(depth)/2);
          voxel_box.translate(Eigen::Vector3f(-static_cast<float>(width)/2, -static_cast<float>(height)/2, -static_cast<float>(depth)/2));
          center += bounding_box.translation;
          //voxel_box.scale(bounding_box.size(1));
          voxel_box.translate(bounding_box.translation);

          // no rotation as the intersection does not work
          // and the point cloud is assumed to be rotated

          // go over all "rays"
          for (unsigned int i = 0; i < considered_rays.size(); i++) {
            if (voxel_box.intersect(this->points[considered_rays[i]])) {
              full_space(n, m, h, w, d) = 1;
            }
          }

          for (unsigned int i = 0; i < considered_points.size(); i++) {
            if (voxel_box.contains(this->points[considered_points[i]])) {
              //this->colors[considered_points[i]] = Eigen::Vector3i(0, 0, 255);
              input(n, m, h, w, d) = 1;
              input_count++;
            }
            if (voxel_box.intersect(this->points[considered_points[i]])) {
              part_space(n, m, h, w, d) = 1;
            }
          }
        }
      }
    }

    //std::cout << "[Data] " << input_count << " observed voxels" << std::endl;
  }

  /** \brief Voxelizes free space; needs the corresponding depth image and
   * assumes the point cloud to be aligned (i.e. rotated) with the bounding box sides.
   * \param[in] camera camera used for rendering
   * \param[in] input tensor to hold voxelized point clouds
   * \param[in] full_space tensor to hold voxelized free space
   * \param[in] part_space tensor to hold voxelized (partial) free space
   * \param[in] index of first dimension in dense and depth
   */
  void voxelize_SCALED_NOSPACE_2(const BoundingBox &bounding_box,
    const BoundingBox &bounding_box_small,
      Eigen::Tensor<float, 5, Eigen::RowMajor> &input,
      const int n, const int m) {

    const int height = input.dimension(2);
    const int width = input.dimension(3);
    const int depth = input.dimension(4);

    Box box;
    bounding_box_small.to_box(box);
    std::vector<int> considered_rays;
    std::vector<int> considered_points;

    for (unsigned int i = 0; i < this->points.size(); i++) {
      if (box.intersect(this->points[i])) {
        considered_rays.push_back(i);
      }
      if (box.contains(this->points[i])) {
        considered_points.push_back(i);
        this->colors[i] = Eigen::Vector3i(0, 255, 0);
      }
    }

    //std::cout << "[Data] SCALED 5" << std::endl;
    //std::cout << "[Data] box " << box.min(0) << " - " << box.max(0)
    //  << " " << box.min(1) << " - " << box.max(1) << " " << box.min(2) << " - " << box.max(2) << std::endl;
    //std::cout << "[Data] voxel grid "
    //  << -static_cast<float>(width)/2 + bounding_box.translation(0) << " - " << static_cast<float>(width)/2 + bounding_box.translation(0) << " "
    //  << -static_cast<float>(height)/2 + bounding_box.translation(1) << " - " << static_cast<float>(height)/2 + bounding_box.translation(1) << " "
    //  << -static_cast<float>(depth)/2 + bounding_box.translation(2) << " - " << static_cast<float>(depth)/2 + bounding_box.translation(2)
    //  << std::endl;

    int input_count = 0;
    for (int h = 0; h < height; h++) {
      for (int w = 0; w < width; w++) {
        for (int d = 0; d < depth; d++) {
          // the box corresponding to this voxel
          Eigen::Vector3f min(static_cast<float>(w), static_cast<float>(h), static_cast<float>(d));
          Eigen::Vector3f max(static_cast<float>(w + 1), static_cast<float>(h + 1), static_cast<float>(d + 1));
          Eigen::Vector3f center(w + 0.5f, h + 0.5f, d + 0.5f);
          Box voxel_box(min, max);

          center += Eigen::Vector3f(-static_cast<float>(width)/2, -static_cast<float>(height)/2, -static_cast<float>(depth)/2);
          voxel_box.translate(Eigen::Vector3f(-static_cast<float>(width)/2, -static_cast<float>(height)/2, -static_cast<float>(depth)/2));
          center += bounding_box.translation;
          //voxel_box.scale(bounding_box.size(1));
          voxel_box.translate(bounding_box.translation);

          // no rotation as the intersection does not work
          // and the point cloud is assumed to be rotated

          // go over all "rays"
          //for (unsigned int i = 0; i < considered_rays.size(); i++) {
          //  if (voxel_box.intersect(this->points[considered_rays[i]])) {
          //    full_space(n, m, h, w, d) = 1;
          //  }
          //}

          for (unsigned int i = 0; i < considered_points.size(); i++) {
            if (voxel_box.contains(this->points[considered_points[i]])) {
              //this->colors[considered_points[i]] = Eigen::Vector3i(0, 0, 255);
              input(n, m, h, w, d) = 1;
              input_count++;
            }
            //if (voxel_box.intersect(this->points[considered_points[i]])) {
            //  part_space(n, m, h, w, d) = 1;
            //}
          }
        }
      }
    }

    //std::cout << "[Data] " << input_count << " observed voxels" << std::endl;
  }

  /** \brief Voxelizes free space; needs the corresponding depth image and
   * assumes the point cloud to be aligned (i.e. rotated) with the bounding box sides.
   * \param[in] camera camera used for rendering
   * \param[in] input tensor to hold voxelized point clouds
   * \param[in] index of first dimension in dense and depth
   */
  void voxelize_SCALED_NOSPACE(Eigen::Tensor<float, 4, Eigen::RowMajor> &input, const int n) {

    const int height = input.dimension(1);
    const int width = input.dimension(2);
    const int depth = input.dimension(3);

    int input_count = 0;
    for (int h = 0; h < height; h++) {
      for (int w = 0; w < width; w++) {
        for (int d = 0; d < depth; d++) {

          // the box corresponding to this voxel
          Eigen::Vector3f min(static_cast<float>(w), static_cast<float>(h), static_cast<float>(d));
          Eigen::Vector3f max(static_cast<float>(w + 1), static_cast<float>(h + 1), static_cast<float>(d + 1));
          Eigen::Vector3f center(w + 0.5f, h + 0.5f, d + 0.5f);
          Box voxel_box(min, max);

          //center += Eigen::Vector3f(-static_cast<float>(width)/2, -static_cast<float>(height)/2, -static_cast<float>(depth)/2);
          //voxel_box.translate(Eigen::Vector3f(-static_cast<float>(width)/2, -static_cast<float>(height)/2, -static_cast<float>(depth)/2));
          //center += bounding_box.translation;
          //voxel_box.translate(bounding_box.translation);

          // no rotation as the intersection does not work
          // and the point cloud is assumed to be rotated

          for (unsigned int i = 0; i < this->points.size(); i++) {
            if (voxel_box.contains(this->points[i])) {
              //this->colors[considered_points[i]] = Eigen::Vector3i(0, 0, 255);
              input(n, h, w, d) = 1;
              input_count++;
            }
          }
        }
      }
    }

    //std::cout << "[Data] " << input_count << " observed voxels" << std::endl;
  }

  /** \brief Get a point.
   * \param[in] n
   * \return point n
   */
  Eigen::Vector3f get(int n) const {
    assert(n < this->num_points() && n >= 0);
    return this->points[n];
  }

private:
  /** \brief The points of the point cloud. */
  std::vector<Eigen::Vector3f> points;
  /** \brief Colors of the points. */
  std::vector<Eigen::Vector3i> colors;

};

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

  /** \brief Check if key exists.
   * \param[in] key
   * \return is set
   */
  bool is_set(std::string key) {
    return this->root.isMember(key);
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

  /** \brief Get the value of the given key as bool.
   * \param[in] key key to retrieve
   * \return value
   */
  bool get_bool(std::string key) {
    if (!this->root.isMember(key)) {
      std::cout << "[Error] key " << key << " not found" << std::endl;
      exit(1);
    }

    return this->root[key].asBool();
  }

  /** \brief Get the value as a vector of integers.
   * \param[in] key key to retrieve
   * \param[out] vector vector to write to
   */
  void get_int_vector(std::string key, std::vector<int> &vector) {
    if (!this->root.isMember(key)) {
      std::cout << "[Error] key " << key << " not found" << std::endl;
      exit(1);
    }

    const Json::Value elements = root[key];
    for (int i = 0; i < elements.size(); ++i) {
      vector.push_back(elements[i].asInt());
    }
  }

private:

  /** \brief Root element. */
  Json::Value root;
};

/** \brief Perform ICP using a point-to-point distance.
 * \param[in] point_cloud_from
 * \param[in] point_cloud_to
 * \param[out] rotation
 * \param[out] translation
 */
void point_to_point_icp(const PointCloud &point_cloud_from, const PointCloud &point_cloud_to,
    Eigen::Matrix3f &rotation, Eigen::Vector3f &translation, float &residual, int &inliers) {

  double* M = new double[3*point_cloud_to.num_points()];
  for (unsigned int i = 0; i < point_cloud_to.num_points(); i++) {
    for (int d = 0; d < 3; d++) {
      M[i*3 + d] = point_cloud_to.get(i)(d);
    }
  }

  double* T = new double[3*point_cloud_from.num_points()];
  for (unsigned int i = 0; i < point_cloud_from.num_points(); i++) {
    for (int d = 0; d < 3; d++) {
      T[i*3 + d] = point_cloud_from.get(i)(d);
    }
  }

  unsigned int I = std::min(point_cloud_from.num_points(), point_cloud_to.num_points());
  translation(0) = 0; translation(1) = 0; translation(2) = 0;

  for (unsigned int i = 0; i < I; i++) {
    translation += point_cloud_to.get(i) - point_cloud_from.get(i);
  }

  translation /= I;

  Matrix R = Matrix::eye(3);
  Matrix t(3, 1);
  t.val[0][0] = translation(0);
  t.val[1][0] = translation(1);
  t.val[2][0] = translation(2);

  IcpPointToPoint icp(M, point_cloud_to.num_points(), 3);
  icp.setMaxIterations(1000);
  icp.setMinDeltaParam(0.0001);

  residual = icp.fit(T, point_cloud_from.num_points(), R, t, 54);
  inliers = icp.getInlierCount();

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      rotation(i, j) = R.val[i][j];
    }

    translation(i) = t.val[i][0];
  }

  delete[] M;
  delete[] T;
}

void extract_ALTERNATIVE(const PointCloud &point_cloud, const BoundingBox &bounding_box, PointCloud &bounding_box_point_cloud) {
  for (unsigned int l = 0; l < point_cloud.num_points(); l++) {
    bool contained = true;
    for (int d = 0; d < 3; d++) {
      if (point_cloud.get(l)(d) < bounding_box.translation(d) - bounding_box.size(d)/2) {
        contained = false;
      }
      if (point_cloud.get(l)(d) > bounding_box.translation(d) + bounding_box.size(d)/2) {
        contained = false;
      }
    }

    if (contained) {
      bounding_box_point_cloud.add_point(point_cloud.get(l));
    }
  }
}

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cout << "[Error] Usage: voxelize_kitti config.json" << std::endl;
    exit(1);
  }

  // Read the JSON configuration file.
  ////////////////////////////////////

  Json::Value root;
  Json::Reader reader;

  boost::filesystem::path json_file(argv[1]);
  std::string set = json_file.stem().string();
  std::cout << "[Data] processing " << set << std::endl;

  std::ifstream json_stream(json_file.string(), std::ifstream::in | std::ifstream::binary);
  reader.parse(json_stream, root, false);
  JSON json(root);

  int height = json.get_int("height");
  int width = json.get_int("width");
  int depth = json.get_int("depth");
  int multiplier = json.get_int("multiplier");

  // assumes padding fixed (not random) and the same in all dimensions
  float padding = json.get_float("padding");

  int min_points = json.get_int("min_points");

  bool gt = json.get_bool("gt");
  int gt_range = 0;
  int gt_skip = 1;

  if (gt) {
    gt_range = json.get_int("gt_range");
    gt_skip = json.get_int("gt_skip");
  }

  boost::filesystem::path bin_directory(json.get_string("velodyne_directory"));
  if (!boost::filesystem::is_directory(bin_directory)) {
    std::cout << "[Error] Directory " << bin_directory.string() << " not found" << std::endl;
    exit(1);
  }

  boost::filesystem::path bounding_box_directory = boost::filesystem::path(json.get_string("bounding_box_directory")
    + "_" + set + "_" + std::to_string(multiplier) + "_" + std::to_string(height) + "x" + std::to_string(width) + "x" + std::to_string(depth));
  if (!boost::filesystem::is_directory(bounding_box_directory)) {
    std::cout << "[Error] Directory " << bounding_box_directory.string() << " not found" << std::endl;
    exit(1);
  }

  boost::filesystem::path num_bounding_box_file = boost::filesystem::path(json.get_string("num_bounding_box_file")
    + "_" + set + "_" + std::to_string(multiplier) + "_" + std::to_string(height) + "x" + std::to_string(width) + "x" + std::to_string(depth) + ".txt");
  if (!boost::filesystem::is_regular_file(num_bounding_box_file)) {
    std::cout << "[Error] file " << num_bounding_box_file.string() << " not found" << std::endl;
    exit(1);
  }

  boost::filesystem::path velodyne_ply_directory = boost::filesystem::path(json.get_string("velodyne_ply_directory")
    + "_" + set + "_" + std::to_string(multiplier) + "_" + std::to_string(height) + "x" + std::to_string(width) + "x" + std::to_string(depth));
  if (!boost::filesystem::is_directory(velodyne_ply_directory)) {
    boost::filesystem::create_directories(velodyne_ply_directory);
    std::cout << "[Data] created " << velodyne_ply_directory.string() << std::endl;
  }

  boost::filesystem::path bounding_box_txt_directory = boost::filesystem::path(json.get_string("bounding_box_txt_directory")
    + "_" + set + "_" + std::to_string(multiplier) + "_" + std::to_string(height) + "x" + std::to_string(width) + "x" + std::to_string(depth));
  if (!boost::filesystem::is_directory(bounding_box_txt_directory)) {
    boost::filesystem::create_directories(bounding_box_txt_directory);
    std::cout << "[Data] created " << bounding_box_txt_directory.string() << std::endl;
  }

  boost::filesystem::path bounding_box_gt_txt_directory = boost::filesystem::path(json.get_string("bounding_box_gt_txt_directory")
    + "_" + set + "_" + std::to_string(multiplier) + "_" + std::to_string(height) + "x" + std::to_string(width) + "x" + std::to_string(depth));
  if (gt && !boost::filesystem::is_directory(bounding_box_gt_txt_directory)) {
    std::cout << "[Error] Directory " << bounding_box_gt_txt_directory.string() << " not found" << std::endl;
    exit(1);
  }

  boost::filesystem::path velodyne_gt_bin_directory = boost::filesystem::path(json.get_string("velodyne_gt_bin_directory")
    + "_" + set + "_" + std::to_string(multiplier) + "_" + std::to_string(height) + "x" + std::to_string(width) + "x" + std::to_string(depth));
  if (gt && !boost::filesystem::is_directory(velodyne_gt_bin_directory)) {
    std::cout << "[Error] Directory " << velodyne_gt_bin_directory.string() << " not found" << std::endl;
    exit(1);
  }

  boost::filesystem::path velodyne_gt_txt_directory = boost::filesystem::path(json.get_string("velodyne_gt_txt_directory")
    + "_" + set + "_" + std::to_string(multiplier) + "_" + std::to_string(height) + "x" + std::to_string(width) + "x" + std::to_string(depth));
  if (gt && !boost::filesystem::is_directory(velodyne_gt_txt_directory)) {
    boost::filesystem::create_directories(velodyne_gt_txt_directory);
  }

  boost::filesystem::path velodyne_individual_gt_txt_directory = boost::filesystem::path(json.get_string("velodyne_individual_gt_txt_directory")
    + "_" + set + "_" + std::to_string(multiplier) + "_" + std::to_string(height) + "x" + std::to_string(width) + "x" + std::to_string(depth));
  if (gt && !boost::filesystem::is_directory(velodyne_individual_gt_txt_directory)) {
    boost::filesystem::create_directories(velodyne_individual_gt_txt_directory);
  }

  boost::filesystem::path bounding_boxes_used_directory;
  if (json.is_set("bounding_boxes_used_directory")) {
    bounding_boxes_used_directory = boost::filesystem::path(json.get_string("bounding_boxes_used_directory")
      + "_" + set + "_" + std::to_string(multiplier) + "_" + std::to_string(height) + "x" + std::to_string(width) + "x" + std::to_string(depth));

    if (!boost::filesystem::is_directory(bounding_boxes_used_directory)) {
      std::cout << "[Error] Directory " << bounding_boxes_used_directory.string() << " not found" << std::endl;
      exit(1);
    }
  }

  // Get the number of bounding boxes and setup output tensors.
  /////////////////////////////////////////////////////////////

  int num_bounding_boxes = read_num_bounding_boxes(num_bounding_box_file.string());
  Eigen::Tensor<float, 4, Eigen::RowMajor> input(num_bounding_boxes, height, width, depth);
  input.setZero();
  std::cout << "[Data] initialized input" << std::endl;

  Eigen::Tensor<float, 4, Eigen::RowMajor> input_sdf(num_bounding_boxes, height, width, depth);
  input_sdf.setZero();
  std::cout << "[Data] initialized input_sdf" << std::endl;

  Eigen::Tensor<float, 4, Eigen::RowMajor> full_space(num_bounding_boxes, height, width, depth);
  full_space.setZero();
  std::cout << "[Data] initialized full_space" << std::endl;

  Eigen::Tensor<float, 4, Eigen::RowMajor> part_space(num_bounding_boxes, height, width, depth);
  std::vector<BoundingBox> all_bounding_boxes(num_bounding_boxes);

  //std::cout << "[Data] " << 2*gt_range/gt_skip + 1 << " channels for gt" << std::endl;
  //Eigen::Tensor<float, 5, Eigen::RowMajor> input_gt(num_bounding_boxes, 2*gt_range/gt_skip + 1, height, width, depth);
  //input_gt.setZero();
  //std::cout << "[Data] initialized input_gt" << std::endl;

  //Eigen::Tensor<float, 5, Eigen::RowMajor> input_sdf_gt(num_bounding_boxes, 2*gt_range/gt_skip + 1, height, width, depth);
  //input_sdf_gt.setZero();

  Eigen::Tensor<float, 4, Eigen::RowMajor> input_combined_gt(num_bounding_boxes, height, width, depth);
  input_combined_gt.setZero();
  std::cout << "[Data] initialized input_combined_gt" << std::endl;

  //Eigen::Tensor<float, 4, Eigen::RowMajor> input_sdf_combined_gt(num_bounding_boxes, height, width, depth);
  //input_sdf_combined_gt.setZero();

  //Eigen::Tensor<float, 5, Eigen::RowMajor> full_space_gt(num_bounding_boxes, 2*gt_range/gt_skip + 1, height, width, depth);
  //full_space_gt.setZero();
  //Eigen::Tensor<float, 5, Eigen::RowMajor> part_space_gt(num_bounding_boxes, 2*gt_range/gt_skip + 1, height, width, depth);
  //part_space_gt.setZero();

  std::vector<int> indices;
  read_split_file(json.get_string("split_file"), indices);

  int n = 0;
  int not_matched_gt = 0;

  for (unsigned int i = 0; i < indices.size(); i++) {

    PointCloud base_point_cloud;
    assert(base_point_cloud.num_points() == 0);
    std::vector<BoundingBox> bounding_boxes;

    // Read the point cloud.
    ////////////////////////
    {
      std::string bin_file = format_string(bin_directory.string() + "/%06d.bin", indices[i]);
      bool success = PointCloud::from_bin(bin_file, base_point_cloud);

      if (!success) {
        std::cout << "[Error] Could not read bin " << bin_file
          << " (" << indices[i] << ")" << std::endl;
        exit(1);
        //continue;
      }

      std::cout << "[Data] Read " << base_point_cloud.num_points()
          << " points from " << bin_file
          << " (" << indices[i] << ")" << std::endl;

      Eigen::Vector3f min;
      Eigen::Vector3f max;
      base_point_cloud.extents(min, max);

      std::cout << "   extents " << min(0) << "/" << max(0) << " " << min(1)
          << "/" << max(1) << " " << min(2) << "/" << max(2) << std::endl;
    }

    // Read the bounding boxes.
    ///////////////////////////
    {
      std::string bounding_box_file = format_string(bounding_box_directory.string() + "/%06d.txt", indices[i]);
      bool success = read_bounding_boxes(bounding_box_file, bounding_boxes);

      if (!success) {
        std::cout << "[Error] Could not read bounding box file " << bounding_box_file
        << " (" << indices[i] << ")" << std::endl;
        exit(1);
        //continue;
      }

      std::cout << "[Data] Read " << bounding_boxes.size()
          << " bounding boxes from " << bounding_box_file
          << " (" << indices[i] << ")" << std::endl;
    }

    // Mark bounding boxes and write point cloud.
    /////////////////////////////////////////////
//    {
//      PointCloud marked_point_cloud(base_point_cloud);
//      marked_point_cloud.mark();
//      for (unsigned int j = 0; j < bounding_boxes.size(); j++) {
//        marked_point_cloud.mark(bounding_boxes[j]);
//      }
//
//      std::string ply_file = format_string(velodyne_ply_directory.string() + "/%06d.ply", indices[i]);
//      bool success = marked_point_cloud.to_ply(ply_file);
//
//      if (!success) {
//        std::cout << "[Error] Could not write point cloud ply " << ply_file
//          << " (" << indices[i] << ")" << std::endl;
//        exit(1);
//        //continue;
//      }
//
//      std::cout << "[Data] Wrote " << marked_point_cloud.num_points()
//          << " points to " << ply_file
//          << " (" << indices[i] << ")" << std::endl;
//    }

    // Voxelize the bounding boxes and compute free space.
    //////////////////////////////////////////////////////
    for (unsigned int j = 0; j < bounding_boxes.size(); j++) {
      BoundingBox original_bounding_box = bounding_boxes[j];

      // copy to manipulate
      BoundingBox bounding_box = original_bounding_box;
      Eigen::Vector3f bounding_box_size = bounding_box.size;
      Eigen::Vector3f bounding_box_translation = bounding_box.translation;

      // add bounding box
      bounding_box.meta = format_string("%06d", indices[i]);

      std::cout << "[Data] processing bounding box " << n
        << " (" << indices[i] << ")" << std::endl;
      std::cout << "   translation: " << format_vector(bounding_box.translation) << std::endl;
      std::cout << "   rotation: " << format_vector(bounding_box.rotation) << std::endl;
      std::cout << "   size: " << format_vector(bounding_box.size) << std::endl;

      // Get a point cloud we can rotate.
      PointCloud point_cloud(base_point_cloud);

      // Rotate it according to bounding box.
      point_cloud.translate(-bounding_box.translation);

      Eigen::Matrix3f rotation;
      PointCloud::rotation_matrix_y(bounding_box.rotation(1), rotation);

      point_cloud.rotate(rotation);
      bounding_box.translation = rotation*bounding_box.translation;
      bounding_box.rotation = Eigen::Vector3f::Zero(); // reset!
      point_cloud.translate(bounding_box.translation);

      //float padding_factor = 1.2;
      float scale_factor = static_cast<float>(width)/bounding_box.size(0);
      std::cout << "[Data] scaling factor " << scale_factor << std::endl;

      point_cloud.scale(scale_factor);
      bounding_box.translation *= scale_factor;
      bounding_box.size *= scale_factor;

      BoundingBox bounding_box_small = bounding_box;
      bounding_box_small.size(0) /= (1 + 2*padding);
      bounding_box_small.size(1) /= (1 + 2*padding);
      bounding_box_small.size(2) /= (1 + 2*padding);

      std::cout << "[Data] small: " << bounding_box_small.size(0) << " " << bounding_box_small.size(2) << " " << bounding_box_small.size(1)
          << " original: " << bounding_box.size(0) << " " << bounding_box.size(2) << " " << bounding_box.size(1)
          << " (" << indices[i] << ")" << std::endl;

      PointCloud bounding_box_point_cloud;
      extract_ALTERNATIVE(point_cloud, bounding_box_small, bounding_box_point_cloud);
      //point_cloud.extract(bounding_box_small, bounding_box_point_cloud);

      int num_points = bounding_box_point_cloud.num_points();
      if (num_points >= min_points) {

        // !
        all_bounding_boxes[n] = original_bounding_box;
        all_bounding_boxes[n].meta = format_string("%06d", indices[i]);

        // Voxelize free space and point cloud.
        point_cloud.voxelize_SCALED(bounding_box, bounding_box_small, input, input_sdf, full_space, part_space, n);
        std::cout << "[Data] space and input for bounding box " << n
          << " (" << indices[i] << ")" << std::endl;

        bounding_box_point_cloud.translate(-bounding_box.translation);
        bounding_box_point_cloud.translate(Eigen::Vector3f(static_cast<float>(width)/2, static_cast<float>(height)/2, static_cast<float>(depth)/2));

        std::string txt_file = format_string(bounding_box_txt_directory.string() + "/%d.txt", n);
        bool success = bounding_box_point_cloud.to_txt(txt_file);

        if (!success) {
          std::cout << "[Error] could not write bounding box txt " << txt_file << std::endl;
          exit(1);
          //continue;
        }

        // Merge with preceding and following point clouds as "ground truth".
        ////////////////////////////////////////////////////////////////////

        if (gt) {
          boost::filesystem::path bb_file = bounding_box_gt_txt_directory / boost::filesystem::path(format_string("%06d_%d.txt", indices[i], 0));

          bool gt_possible = true;
          if (boost::filesystem::is_regular_file(bb_file)) {

            std::vector<BoundingBox> bounding_boxes_ref;
            read_bounding_boxes(bb_file.string(), bounding_boxes_ref);

            int used_j = -1;
            // TODO
            if (!bounding_boxes_used_directory.empty()) {
              boost::filesystem::path used_file = bounding_boxes_used_directory / boost::filesystem::path(format_string("%06d.txt", indices[i]));

              if (!boost::filesystem::is_regular_file(used_file)) {
                std::cout << "[Error] file " << used_file << " not found"
                  << " (" << indices[i] << ")" << std::endl;
                exit(1);
              }

              std::ifstream file(used_file.string().c_str());
              std::vector<int> used;
              int u = -1;

              while (file >> u) {
                assert(u == 1 || u == 0);
                used.push_back(u);
              }

              assert(bounding_boxes_ref.size() == used.size());

              used_j = 0;
              int kk = -1;
              for (unsigned int k = 0; k < used.size(); k++) {
                if (used[i] == 1) {
                  kk++;

                  if (kk == j) {
                    used_j = k;
                  }
                }
              }
            }
            else {
              float min_distance = FLT_MAX;
              for (unsigned int k = 0; k < bounding_boxes.size(); k++) {
                float distance = (bounding_boxes_ref[k].size - bounding_boxes[j].size).norm();

                if (distance < min_distance) {
                  min_distance = distance;
                  used_j = k;
                }
              }
            }

            assert(used_j >= 0);

            BoundingBox bounding_box_ref = bounding_boxes_ref[used_j];
            Eigen::Vector3f bounding_box_size_ref = bounding_box_ref.size;

            boost::filesystem::path bin_file = velodyne_gt_bin_directory / boost::filesystem::path(format_string("%06d_%d.bin", indices[i], 0));

            if (!boost::filesystem::is_regular_file(bin_file)) {
              std::cout << "[Error] " << bin_file.string() << " not found" << std::endl;
              exit(1);
            }

            PointCloud point_cloud_ref;
            bool success = PointCloud::from_bin(bin_file.string(), point_cloud_ref);

            if (!success) {
              std::cout << "[Error] could not read" << bin_file.string() << std::endl;
              exit(1);
            }

            std::cout << "[Data] ATTENTION:" << std::endl;
            std::cout << "   " << bounding_box_translation(0) << " " << bounding_box_translation(1) << " " << bounding_box_translation(2) << std::endl;
            std::cout << "   " << bounding_box_ref.translation(0) << " " << bounding_box_ref.translation(1) << " " << bounding_box_ref.translation(2) << std::endl;
            std::cout << "   " << bounding_box_size(0) << " " << bounding_box_size(1) << " " << bounding_box_size(2) << std::endl;
            std::cout << "   " << bounding_box_ref.size(0) << " " << bounding_box_ref.size(1) << " " << bounding_box_ref.size(2) << std::endl;

            PointCloud bounding_box_point_cloud_ref;
            //extract_bounding_box(point_cloud_ref, bounding_box_ref, width, height, depth, padding, bounding_box_point_cloud);

            point_cloud_ref.translate(-bounding_box_ref.translation);

            // Rotate according to bounding box.
            Eigen::Matrix3f rotation;
            PointCloud::rotation_matrix_y(bounding_box_ref.rotation(1), rotation);
            point_cloud_ref.rotate(rotation);

            bounding_box_ref.translation = rotation*bounding_box_ref.translation;
            bounding_box_ref.rotation = Eigen::Vector3f::Zero(); // reset!
            point_cloud_ref.translate(bounding_box_ref.translation);

            //float padding_factor = 1.2;
            float scale_factor = static_cast<float>(width)/bounding_box_ref.size(0);
            point_cloud_ref.scale(scale_factor);
            bounding_box_ref.translation *= scale_factor;
            bounding_box_ref.size *= scale_factor;

            BoundingBox bounding_box_ref_small = bounding_box_ref;
            bounding_box_ref_small.size(0) /= (1 + 2*padding);
            bounding_box_ref_small.size(1) /= (1 + 2*padding);
            bounding_box_ref_small.size(2) /= (1 + 2*padding);

            int m = 0;
            // Save the bounding box cloud.
            //point_cloud_ref.extract(bounding_box_ref_small, bounding_box_point_cloud_ref);
            extract_ALTERNATIVE(point_cloud_ref, bounding_box_ref_small, bounding_box_point_cloud_ref);

            bounding_box_point_cloud_ref.translate(-bounding_box_ref.translation);
            bounding_box_point_cloud_ref.translate(Eigen::Vector3f(static_cast<float>(width)/2, static_cast<float>(height)/2, static_cast<float>(depth)/2));

            float residual = 1000;
            int inliers = 0;
            Eigen::Matrix3f icp_rotation;
            Eigen::Vector3f icp_translation;
            point_to_point_icp(bounding_box_point_cloud_ref, bounding_box_point_cloud, icp_rotation, icp_translation, residual, inliers);

            std::cout << "[Data] translation" << std::endl;
            std::cout << "   " << icp_translation(0) << " " << icp_translation(1) << " " << icp_translation(2) << std::endl;
            std::cout << "[Data] rotation" << std::endl;
            std::cout << "   " << icp_rotation(0, 0) << " " << icp_rotation(0, 1) << " " << icp_rotation(0, 2) << std::endl;
            std::cout << "   " << icp_rotation(1, 0) << " " << icp_rotation(1, 1) << " " << icp_rotation(1, 2) << std::endl;
            std::cout << "   " << icp_rotation(2, 0) << " " << icp_rotation(2, 1) << " " << icp_rotation(2, 2) << std::endl;

            if (residual < 0.001 && inliers >= bounding_box_point_cloud.num_points()/2) {
              std::cout << "[Data] GT residual " << residual << "/" << inliers << "/" << bounding_box_point_cloud_ref.num_points() << " for ground truth of " << n
                << " (" << indices[i] << ")" << std::endl;

              bounding_box_point_cloud_ref.rotate(icp_rotation);
              bounding_box_point_cloud_ref.translate(icp_translation);

              point_cloud_ref.translate(-bounding_box_ref.translation);
              point_cloud_ref.translate(Eigen::Vector3f(static_cast<float>(width)/2, static_cast<float>(height)/2, static_cast<float>(depth)/2));

              point_cloud_ref.rotate(icp_rotation);
              point_cloud_ref.translate(icp_translation);

              point_cloud_ref.translate(-Eigen::Vector3f(static_cast<float>(width)/2, static_cast<float>(height)/2, static_cast<float>(depth)/2));
              point_cloud_ref.translate(bounding_box_ref.translation);

              //point_cloud_ref.voxelize_SCALED_NOSPACE_2(bounding_box_ref, bounding_box_ref_small, input_gt, n, m);

              std::string txt_file = format_string(velodyne_individual_gt_txt_directory.string() + "/%d_%d_%d.txt", n, 0, m);
              success = bounding_box_point_cloud_ref.to_txt(txt_file);

              if (!success) {
                std::cout << "[Error] could not write " << txt_file << std::endl;
                exit(1);
                //continue;
              }

              m++;
              for (int k = -gt_range; k <= gt_range; k += gt_skip) {
                if (k == 0) {
                  continue;
                }

                boost::filesystem::path bb_file_k = bounding_box_gt_txt_directory / boost::filesystem::path(format_string("%06d_%d.txt", indices[i], k));
                boost::filesystem::path bin_file_k = velodyne_gt_bin_directory / boost::filesystem::path(format_string("%06d_%d.bin", indices[i], k));

                if (!boost::filesystem::is_regular_file(bb_file_k) || !boost::filesystem::is_regular_file(bin_file_k)) {
                  std::cout << "[Data] no ground truth " << k << " for bounding box " << n
                    << " (" << indices[i] << ")" << std::endl;
                  continue;
                }

                std::vector<BoundingBox> bounding_boxes_k;
                read_bounding_boxes(bb_file_k.string(), bounding_boxes_k);

                int j_k = -1;
                for (int jj = 0; jj < bounding_boxes_k.size(); jj++) {
                  if (bounding_box_size_ref.isApprox(bounding_boxes_k[jj].size)) {
                    j_k = jj;
                  }
                }

                if (j_k < 0) {
                  std::cout << "[Data] no bounding box in frame " << k << " for bounding box " << n
                    << " (" << indices[i] << ")" << std::endl;
                  continue;
                }

                BoundingBox bounding_box_k = bounding_boxes_k[j_k];

                PointCloud point_cloud_k;
                bool success = PointCloud::from_bin(bin_file_k.string(), point_cloud_k);

                if (!success) {
                  std::cout << "[Error] could not read " << bin_file_k.string() << std::endl;
                  exit(1);
                }

                PointCloud bounding_box_point_cloud_k;
                //extract_bounding_box(point_cloud_k, bounding_box_k, width, height, depth, padding, bounding_box_point_cloud_k);

                point_cloud_k.translate(-bounding_box_k.translation);

                // Rotate according to bounding box.
                Eigen::Matrix3f rotation;
                PointCloud::rotation_matrix_y(bounding_box_k.rotation(1), rotation);
                point_cloud_k.rotate(rotation);

                bounding_box_k.translation = rotation*bounding_box_k.translation;
                bounding_box_k.rotation = Eigen::Vector3f::Zero(); // reset!
                point_cloud_k.translate(bounding_box_k.translation);

                //float padding_factor = 1.2;
                float scale_factor = static_cast<float>(width)/bounding_box_k.size(0);
                point_cloud_k.scale(scale_factor);
                bounding_box_k.translation *= scale_factor;
                bounding_box_k.size *= scale_factor;

                BoundingBox bounding_box_small_k = bounding_box_k;
                bounding_box_small_k.size(0) /= (1 + 2*padding);
                bounding_box_small_k.size(1) /= (1 + 2*padding);
                bounding_box_small_k.size(2) /= (1 + 2*padding);

                // Save the bounding box cloud.
                //point_cloud_k.extract(bounding_box_small_k, bounding_box_point_cloud_k);
                extract_ALTERNATIVE(point_cloud_k, bounding_box_small_k, bounding_box_point_cloud_k);

                bounding_box_point_cloud_k.translate(-bounding_box_k.translation);
                bounding_box_point_cloud_k.translate(Eigen::Vector3f(static_cast<float>(width)/2, static_cast<float>(height)/2, static_cast<float>(depth)/2));

                bounding_box_point_cloud_k.rotate(icp_rotation);
                bounding_box_point_cloud_k.translate(icp_translation);

                point_cloud_k.translate(-bounding_box_ref.translation);
                point_cloud_k.translate(Eigen::Vector3f(static_cast<float>(width)/2, static_cast<float>(height)/2, static_cast<float>(depth)/2));

                point_cloud_k.rotate(icp_rotation);
                point_cloud_k.translate(icp_translation);

                point_cloud_k.translate(-Eigen::Vector3f(static_cast<float>(width)/2, static_cast<float>(height)/2, static_cast<float>(depth)/2));
                point_cloud_k.translate(bounding_box_ref.translation);

                //point_cloud_k.voxelize_SCALED_NOSPACE_2(bounding_box_k, bounding_box_small_k, input_gt, n, m);

                std::string txt_file = format_string(velodyne_individual_gt_txt_directory.string() + "/%d_%d_%d.txt", n, k, m);
                success = bounding_box_point_cloud_k.to_txt(txt_file);

                if (!success) {
                  std::cout << "[Error] could not write " << txt_file << std::endl;
                  exit(1);
                  //continue;
                }

                //Eigen::Matrix3f rotation;
                //Eigen::Vector3f translation;
                //point_to_plane_icp(bounding_box_point_cloud_k, bounding_box_point_cloud, rotation, translation);

                //bounding_box_point_cloud_k.rotate(rotation);
                //bounding_box_point_cloud_k.translate(translation);

                bounding_box_point_cloud_ref.merge(bounding_box_point_cloud_k);

                m++; // !
              }

              std::cout << "[Data] ground truth for bounding box " << n
                    << " (" << indices[i] << ")" << std::endl;
              txt_file = velodyne_gt_txt_directory.string() + "/" + format_string("%d.txt", n);
              success = bounding_box_point_cloud_ref.to_txt(txt_file);

              if (!success) {
                std::cout << "[Error] could not write " << txt_file << std::endl;
                exit(1);
                //continue;
              }

              bounding_box_point_cloud_ref.voxelize_SCALED_NOSPACE(input_combined_gt, n);
            }
            else {
              gt_possible = false;
              not_matched_gt++;
              std::cout << "[Data] NO-GT residual " << residual << "/" << inliers << "/" << bounding_box_point_cloud_ref.num_points() << " for ground truth of " << n
                << " (" << indices[i] << ")" << std::endl;
            }
          }
          else {
            gt_possible = false;
            std::cout << "[Data] NO-GT" << std::endl;
          }

          if (!gt_possible) {
            std::cout << "[Data] FAKE ground truth for bounding box " << n
                    << " (" << indices[i] << ")" << std::endl;

            //point_cloud.voxelize_SCALED_NOSPACE_2(bounding_box, bounding_box_small, input_gt, n, 0);

            std::string txt_file = velodyne_gt_txt_directory.string() + "/" + format_string("%d.txt", n);
            success = bounding_box_point_cloud.to_txt(txt_file);

            if (!success) {
              std::cout << "[Error] could not write " << txt_file << std::endl;
              exit(1);
              //continue;
            }

            txt_file = format_string(velodyne_individual_gt_txt_directory.string() + "/%d_%d_%d.txt", n, 0, 0);
            success = bounding_box_point_cloud.to_txt(txt_file);

            if (!success) {
              std::cout << "[Error] could not write " << txt_file << std::endl;
              exit(1);
              //continue;
            }

            bounding_box_point_cloud.voxelize_SCALED_NOSPACE(input_combined_gt, n);
          }
        }

        // Increment the current bounding box index.
        ////////////////////////////////////////////
        n++;

        if (n > num_bounding_boxes) {
          std::cout << "[Error] invalid number of bounding boxes: " << n << " processed; " << num_bounding_boxes << " expected" << std::endl;
          exit(1);
        }
      }
    }
  }

  // Write the voxelizations to HDF 5 files.
  //////////////////////////////////////////
  {
    boost::filesystem::path input_file = boost::filesystem::path(
      json.get_string("input_file") + "_" + set + "_" + std::to_string(multiplier) + "_" + std::to_string(height) + "x" + std::to_string(width) + "x" + std::to_string(depth) + ".h5"
    );
    write_hdf5(input_file.string(), input);
    std::cout << "[Data] Wrote " << input_file.string() << std::endl;

    boost::filesystem::path input_sdf_file = boost::filesystem::path(
      json.get_string("input_sdf_file") + "_" + set + "_" + std::to_string(multiplier) + "_" + std::to_string(height) + "x" + std::to_string(width) + "x" + std::to_string(depth) + ".h5"
    );
    write_hdf5(input_sdf_file.string(), input_sdf);
    std::cout << "[Data] Wrote " << input_sdf_file.string() << std::endl;

    boost::filesystem::path full_space_file = boost::filesystem::path(
      json.get_string("full_space_file") + "_" + set + "_" + std::to_string(multiplier) + "_" + std::to_string(height) + "x" + std::to_string(width) + "x" + std::to_string(depth) + ".h5"
    );
    write_hdf5(full_space_file.string(), full_space);
    std::cout << "[Data] Wrote " << full_space_file.string() << std::endl;

    boost::filesystem::path part_space_file = boost::filesystem::path(
      json.get_string("part_space_file") + "_" + set + "_" + std::to_string(multiplier) + "_" + std::to_string(height) + "x" + std::to_string(width) + "x" + std::to_string(depth) + ".h5"
    );
    write_hdf5(part_space_file.string(), part_space);
    std::cout << "[Data] Wrote " << part_space_file.string() << std::endl;

    // !
    if (gt) {
      //boost::filesystem::path input_gt_file = boost::filesystem::path(
      //  json.get_string("input_gt_file") + "_" + set + "_" + std::to_string(multiplier) + "_" + std::to_string(height) + "x" + std::to_string(width) + "x" + std::to_string(depth) + ".h5"
      //);
      //write_hdf5(input_gt_file.string(), input_gt);
      //std::cout << "[Data] Wrote " << input_gt_file.string() << std::endl;

      //boost::filesystem::path input_sdf_gt_file = boost::filesystem::path(
      //  json.get_string("input_sdf_gt_file") + "_" + set + "_" + std::to_string(multiplier) + "_" + std::to_string(height) + "x" + std::to_string(width) + "x" + std::to_string(depth) + ".h5"
      //);
      //write_hdf5(input_sdf_gt_file.string(), input_sdf_gt);
      //std::cout << "[Data] Wrote " << input_sdf_gt_file.string() << std::endl;

      boost::filesystem::path input_combined_gt_file = boost::filesystem::path(
        json.get_string("input_combined_gt_file") + "_" + set + "_" + std::to_string(multiplier) + "_" + std::to_string(height) + "x" + std::to_string(width) + "x" + std::to_string(depth) + ".h5"
      );
      write_hdf5(input_combined_gt_file.string(), input_combined_gt);
      std::cout << "[Data] Wrote " << input_combined_gt_file.string() << std::endl;

      //boost::filesystem::path input_sdf_combined_gt_file = boost::filesystem::path(
      //  json.get_string("input_sdf_combined_gt_file") + "_" + set + "_" + std::to_string(multiplier) + "_" + std::to_string(height) + "x" + std::to_string(width) + "x" + std::to_string(depth) + ".h5"
      //);
      //write_hdf5(input_sdf_combined_gt_file.string(), input_sdf_combined_gt);
      //std::cout << "[Data] Wrote " << input_sdf_combined_gt_file.string() << std::endl;

      //boost::filesystem::path full_space_gt_file = boost::filesystem::path(
      //  json.get_string("full_space_gt_file") + "_" + set + "_" + std::to_string(multiplier) + "_" + std::to_string(height) + "x" + std::to_string(width) + "x" + std::to_string(depth) + ".h5"
      //);
      //write_hdf5(full_space_gt_file.string(), full_space_gt);
      //std::cout << "[Data] Wrote " << full_space_gt_file.string() << std::endl;

      //boost::filesystem::path part_space_gt_file = boost::filesystem::path(
      //  json.get_string("part_space_gt_file") + "_" + set + "_" + std::to_string(multiplier) + "_" + std::to_string(height) + "x" + std::to_string(width) + "x" + std::to_string(depth) + ".h5"
      //);
      //write_hdf5(part_space_gt_file.string(), part_space_gt);
      //std::cout << "[Data] Wrote " << part_space_gt_file.string() << std::endl;
    }

    // !
    boost::filesystem::path bounding_box_file = boost::filesystem::path(
      json.get_string("bounding_box_file") + "_" + set + "_" + std::to_string(multiplier) + "_" + std::to_string(height) + "x" + std::to_string(width) + "x" + std::to_string(depth) + ".txt"
    );
    write_bounding_boxes(bounding_box_file.string(), all_bounding_boxes);
    std::cout << "[Data] Wrote " << bounding_box_file.string() << std::endl;
  }

  std::cout << "[Data] could not match " << not_matched_gt << " ground truth boundng boxes" << std::endl;

  exit(0);
}
