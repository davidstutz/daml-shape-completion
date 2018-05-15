#ifndef KITTI_H_
#define KITTI_H_

// C++ includes
#include <string>
#include <map>

// Eigen
#include <Eigen/Dense>

// Boost
#include <boost/filesystem.hpp>

/** \brief Tools for evaluation on KITTI.
 */
namespace kitti {
  /** \brief Read all files in a directory matching the given extension.
   * \param[in] directory path to directory
   * \param[out] files read file paths
   * \param[in] extension extension to filter for
   */
  void read_directory(const std::string directory, std::map<int, std::string>& files, const std::string extension) {

    boost::filesystem::path dir(directory);
    boost::filesystem::directory_iterator end;

    files.clear();
    for (boost::filesystem::directory_iterator it(dir); it != end; ++it) {
      if (it->path().extension().string() == extension) {
        int number = std::stoi(it->path().stem().string());
        files.insert(std::pair<int, std::string>(number, it->path().string()));
      }
    }
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
      std::cout << "Could not read " << filepath << std::endl;
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

      // ignore "empty" ones, i.e. everything zero
      if (bounding_box.size.norm() > 1e-4) {
        bounding_boxes.push_back(bounding_box);
      }
    }

    return true;
  }
}

#endif