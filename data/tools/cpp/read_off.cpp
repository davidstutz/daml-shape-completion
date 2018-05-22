#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cfloat>

// Boost
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

// Eigen
#include <Eigen/Dense>

/** \brief Just encapsulating vertices and faces. */
class Mesh {
public:
  /** \brief Empty constructor. */
  Mesh() {

  }

  /** \brief Reading an off file and returning the vertices x, y, z coordinates and the
   * face indices.
   * \param[in] filepath path to the OFF file
   * \param[out] mesh read mesh with vertices and faces
   * \return success
   */
  static bool from_off(const std::string filepath, Mesh& mesh) {

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
    (*out) << this->num_vertices() << " " << this->num_faces() << " 0" << std::endl;

    for (unsigned int v = 0; v < this->num_vertices(); v++) {
      (*out) << this->vertices[v](0) << " " << this->vertices[v](1) << " " << this->vertices[v](2) << std::endl;
    }

    for (unsigned int f = 0; f < this->num_faces(); f++) {
      (*out) << "3 " << this->faces[f](0) << " " << this->faces[f](1) << " " << this->faces[f](2) << std::endl;
    }

    out->close();
    delete out;

    return true;
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

private:

  /** \brief Vertices as (x,y,z)-vectors. */
  std::vector<Eigen::Vector3f> vertices;

  /** \brief Faces as list of vertex indices. */
  std::vector<Eigen::Vector3i> faces;
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

  Mesh mesh;
  Mesh::from_off(input.string(), mesh);

  std::cout << "Read " << input << "." << std::endl;
  std::cout << mesh.num_vertices() << " vertices and " << mesh.num_faces() << " faces." << std::endl;

  return 0;
}