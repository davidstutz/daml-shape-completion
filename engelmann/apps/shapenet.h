#ifndef SHAPENET_H_
#define SHAPENET_H_

// C++ includes
#include <string>
#include <map>

// Eigen
#include <Eigen/Dense>

// Boost
#include <boost/filesystem.hpp>

/** \brief Tools for evaluation on KITTI.
 */
namespace shapenet {
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
}

#endif