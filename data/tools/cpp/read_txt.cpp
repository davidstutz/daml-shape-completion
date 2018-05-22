#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cfloat>

// Eigen
#include <Eigen/Dense>

// Boost
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

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

    for (unsigned int i = 0; i < point_cloud.points.size(); i++) {
      this->points.push_back(point_cloud.points[i]);
    }
  }

  /** \brief Destructor. */
  ~PointCloud() {

  }

  /** \brief Read point cloud from txt file.
   * \param[in] filepath path to file to read
   * \param[out] point_cloud
   * \return success
   */
  static bool from_txt(const std::string &filepath, PointCloud &point_cloud) {
    std::ifstream file(filepath.c_str());
    std::string line;
    std::stringstream ss;

    std::getline(file, line);
    ss << line;

    int n_points = 0;
    ss >> n_points;

    if (n_points < 0) {
      return false;
    }

    for (int i = 0; i < n_points; i++) {
      std::getline(file, line);

      ss.clear();
      ss.str("");
      ss << line;

      Eigen::Vector3f point(0, 0, 0);
      ss >> point(0);
      ss >> point(1);
      ss >> point(2);

      point_cloud.add_point(point);
    }

    return true;
  }

  /** \brief Add a point to the point cloud.
   * \param[in] point point to add
   */
  void add_point(const Eigen::Vector3f &point) {
    this->points.push_back(point);
  }

  /** \brief Get number of points.
   * \return number of points
   */
  unsigned int num_points() const {
    return this->points.size();
  }

private:
  /** \brief The points of the point cloud. */
  std::vector<Eigen::Vector3f> points;

};

int main(int argc, char** argv) {
  boost::program_options::options_description desc("Allowed options");
  desc.add_options()
      ("help", "produce help message")
      ("input",  boost::program_options::value<std::string>(), "path to input OFF file");

  boost::program_options::positional_options_description positionals;
  positionals.add("input", 1);

  boost::program_options::variables_map parameters;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(desc).positional(positionals).run(), parameters);
  boost::program_options::notify(parameters);

  if (parameters.find("help") != parameters.end()) {
    std::cout << desc << std::endl;
    return 1;
  }

  boost::filesystem::path input(parameters["input"].as<std::string>());
  if (!boost::filesystem::is_regular_file(input)) {
    std::cout << "Input file does not exist." << std::endl;
    return 1;
  }

  PointCloud point_cloud;
  PointCloud::from_txt(input.string(), point_cloud);

  std::cout << "Read " << input << "." << std::endl;
  std::cout << point_cloud.num_points() << " points." << std::endl;

  return 0;
}