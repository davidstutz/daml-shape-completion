#include <iostream>
#include <fstream>
#include <cerrno>
#include <cstring>

#include "helpers.h"
#include "exception.h"
#include "modelloader.h"

REMESHER_NAMESPACE_BEGIN

TriangleMeshPtr
ModelLoader::load_model (std::string const& filename)
{
  /* Detect model file. Use lazy detection with filename. */
  if (filename.size() <= 2)
    throw Exception("Cannot detect model format from filename");

  std::string extension = filename.substr(filename.size() - 2);

  if (extension == ".m")
    return ModelLoader::load_m_model(filename);

  if (filename.size() <= 4)
    throw Exception("Cannot detect model format from filename");

  extension = filename.substr(filename.size() - 4);
  if (extension == ".off")
    return ModelLoader::load_off_model(filename);
  else if (extension == ".ply")
    return ModelLoader::load_ply_model(filename);
  else if (extension == ".obj")
    return ModelLoader::load_obj_model(filename);

  throw Exception("Cannot detect model format from filename");
}

/* ---------------------------------------------------------------- */

TriangleMeshPtr
ModelLoader::load_off_model (std::string const& filename)
{
  /* Open file. */
  std::fstream input(filename.c_str());
  if (input.fail())
    throw Exception("Cannot open file: " + std::string(::strerror(errno)));

  /* Start parsing. */
  std::string buffer;
  bool parse_normals;
  input >> buffer; /* Read "OFF" file signature. */
  if (buffer == "NOFF")
    parse_normals = true;
  else if (buffer == "OFF")
    parse_normals = false;
  else
  {
    input.close();
    throw Exception("File format not recognized as OFF-model");
  }

  /* Create a new triangle mesh. */
  TriangleMeshPtr mesh = TriangleMesh::create();
  MeshVertexList& vertices = mesh->get_vertices();
  MeshFaceList& faces = mesh->get_faces();
  MeshNormalList& vertex_normals = mesh->get_vertex_normals();

  /* Clear the model data and init some values. */
  std::size_t num_vertices = 0;
  std::size_t num_faces = 0;
  std::size_t num_edges = 0;

  /* Read vertex, face and edge information. */
  /* TODO: Edges are ignored. */
  input >> num_vertices >> num_faces >> num_edges;

  vertices.reserve(num_vertices);
  faces.reserve(num_faces * 3);
  vertex_normals.reserve(num_vertices);

  /* Read vertices. */
  for (std::size_t i = 0; i < num_vertices; ++i)
  {
    float x, y, z;
    input >> x >> y >> z;
    vertices.push_back(Vec3f(x, y, z));

    /* Also read vertex normals if present. */
    if (parse_normals)
    {
      input >> x >> y >> z;
      vertex_normals.push_back(Vec3f(x, y, z));
    }
  }

  /* Read faces. */
  for (std::size_t i = 0; i < num_faces; ++i)
  {
    std::size_t n_vertices;
    input >> n_vertices;

    if (n_vertices == 3)
    {
      /* Polygon is a triangle, which is simply appended to the mesh. */
      MeshVIndex vidx[3];
      bool indices_good = true;
      for (int j = 0; j < 3; ++j)
      {
        input >> vidx[j];
        if (vidx[j] >= num_vertices)
        {
          std::cout << "OFF Loader: Warning: Face " << i << " has invalid vertex "
              << vidx[j] << ", skipping face." << std::endl;
          indices_good = false;
        }
      }

      if (indices_good)
        for (int j = 0; j < 3; ++j)
          faces.push_back(vidx[j]);
    }
    else if (n_vertices == 4)
    {
      /* Polygon is a quad and has to be converted to 2 triangles. */
      MeshVIndex vidx[4];
      bool indices_good = true;
      for (int j = 0; j < 4; ++j)
      {
        input >> vidx[j];
        if (vidx[j] >= num_vertices)
        {
          std::cout << "OFF Loader: Warning: Face " << i << " has invalid vertex "
              << vidx[j] << ", skipping face." << std::endl;
          indices_good = false;
        }
      }

      if (indices_good)
      {
        for (int j = 0; j < 3; ++j)
          faces.push_back(vidx[j]);
        for (int j = 0; j < 3; ++j)
          faces.push_back(vidx[(j + 2) % 4]);
      }
    }
    else
    {
      std::cout << "Verts: " << num_vertices << ", faces: "
          << num_faces << ", line: " << 2 + num_vertices + i << std::endl;
      std::cout << "Warning: Polygon with more than 4 edges!" << std::endl;
      ::exit(1);
    }
  }

  /* Close file stream. */
  input.close();

  //mesh->ensure_normals();
  //mesh->scale_and_center();
  return mesh;
}

/* ---------------------------------------------------------------- */

TriangleMeshPtr
ModelLoader::load_ply_model (std::string const& filename)
{
  enum PlyFormat
  {
    PLY_ASCII,
    PLY_BINARY_LE,
    PLY_BINARY_BE,
    PLY_UNKNOWN
  };

  enum PlyVertexElement
  {
    PLY_V_FLOAT_X = 0,
    PLY_V_FLOAT_Y,
    PLY_V_FLOAT_Z,
    PLY_V_UCHAR_R = 3,
    PLY_V_UCHAR_G,
    PLY_V_UCHAR_B,
    PLY_V_FLOAT_IGNORE,
    PLY_V_INT_IGNORE,
    PLY_V_BYTE_IGNORE
  };

  enum PlyFaceElement
  {
    PLY_F_VERTEX_INDICES,
    PLY_F_INT_IGNORE,
    PLY_F_BYTE_IGNORE
  };

  /* Open file. */
  std::fstream input(filename.c_str());
  if (input.fail())
    throw Exception("Cannot open file: " + std::string(::strerror(errno)));

  /* Start parsing. */
  std::string buffer;
  input >> buffer; /* Read "ply" file signature. */
  if (buffer != "ply")
  {
    input.close();
    throw Exception("File format not recognized as PLY-model");
  }

  /* Discard the rest of the line. */
  std::getline(input, buffer);

  PlyFormat ply_format = PLY_UNKNOWN;
  std::size_t num_faces = 0;
  std::size_t num_vertices = 0;
  std::vector<int> v_format;
  std::vector<int> f_format;

  bool critical = false;
  bool reading_verts = false;
  bool reading_faces = false;

  while (input.good())
  {
    std::getline(input, buffer);
    Helpers::clip_string(buffer);
    //std::cout << "Buffer: " << buffer << std::endl;

    if (buffer == "end_header")
      break;

    StringVector header = Helpers::split_string(buffer, ' ');

    /* Determine the format. */
    if (header[0] == "format")
    {
      if (header[1] == "ascii")
        ply_format = PLY_ASCII;
      else if (header[1] == "binary_little_endian")
        ply_format = PLY_BINARY_LE;
      else if (header[1] == "binary_big_endian")
        ply_format = PLY_BINARY_BE;
      else
        ply_format = PLY_UNKNOWN;
    }
    else if (header[0] == "comment")
    {
      std::cout << "PLY Loader: " << buffer << std::endl;
    }
    else if (header[0] == "element")
    {
      reading_faces = false;
      reading_verts = false;
      if (header[1] == "vertex")
      {
        reading_verts = true;
        num_vertices = Helpers::get_int_from_string(header[2]);
      }
      else if (header[1] == "face")
      {
        reading_faces = true;
        num_faces = Helpers::get_int_from_string(header[2]);
      }
      else
        std::cout << "PLY Loader: Element \"" << header[1]
            << "\" not recognized" << std::endl;
    }
    else if (header[0] == "property")
    {
      if (reading_verts)
      {
        /* List of accepted and handled attributes. */
        if (header[1] == "float" || header[1] == "float32")
        {
          if (header[2] == "x")
            v_format.push_back(PLY_V_FLOAT_X);
          else if (header[2] == "y")
            v_format.push_back(PLY_V_FLOAT_Y);
          else if (header[2] == "z")
            v_format.push_back(PLY_V_FLOAT_Z);
          else
            v_format.push_back(PLY_V_FLOAT_IGNORE);
        }
        else if (header[1] == "uchar" || header[1] == "uint8")
        {
          if (header[2] == "r" || header[2] == "red"
              || header[2] == "diffuse_red")
            v_format.push_back(PLY_V_UCHAR_R);
          else if (header[2] == "g" || header[2] == "green"
              || header[2] == "diffuse_green")
            v_format.push_back(PLY_V_UCHAR_G);
          else if (header[2] == "b" || header[2] == "blue"
              || header[2] == "diffuse_blue")
            v_format.push_back(PLY_V_UCHAR_B);
          else
            v_format.push_back(PLY_V_BYTE_IGNORE);
        }
        /* All other data types are ignored. */
        else if (header[1] == "int")
          v_format.push_back(PLY_V_INT_IGNORE);
        else
        {
          std::cout << "PLY Loader: Unrecognized vertex property \""
              << header[1] << " " << header[2] << "\"" << std::endl;
          critical = true;
          continue;
        }
      }
      else if (reading_faces)
      {
        if (header[1] == "list")
          f_format.push_back(PLY_F_VERTEX_INDICES);
        else if (header[1] == "int")
          f_format.push_back(PLY_F_INT_IGNORE);
        else if (header[1] == "uchar")
          f_format.push_back(PLY_F_BYTE_IGNORE);
        else
        {
          std::cout << "PLY Loader: Unrecognized face property \""
              << header[1] << " " << header[2] << "\"" << std::endl;
          critical = true;
        }
      }
    }
  }

  if (critical || num_vertices == 0)
  {
    #if 0
    std::cout << "Critical: " << critical << ", Vertices: " << num_vertices
        << ", Faces: " << num_faces << std::endl
        << "  has x,y,z: " << vert_has_x << vert_has_y << vert_has_z
        << std::endl;
    #endif
    input.close();
    throw Exception("File headers not recognized as PLY format");
  }

  if (ply_format == PLY_UNKNOWN)
  {
    input.close();
    throw Exception("PLY file encoding not recognized by parser");
  }

  /* Create a new triangle mesh. */
  TriangleMeshPtr mesh = TriangleMesh::create();
  MeshVertexList& vertices = mesh->get_vertices();
  MeshFaceList& faces = mesh->get_faces();
  MeshVertexColorList& vcolors = mesh->get_vertex_colors();

  bool want_colors = false;
  for (std::size_t i = 0; i < v_format.size(); ++i)
    if (v_format[i] == PLY_V_UCHAR_R || v_format[i] == PLY_V_UCHAR_G
        || v_format[i] == PLY_V_UCHAR_B)
    {
      want_colors = true;
      break;
    }

  /* Start reading the vertex data. */
  vertices.reserve(num_vertices);
  for (std::size_t i = 0; i < num_vertices; ++i)
  {
    Vec3f vertex(0.0f, 0.0f, 0.0f);
    Vec3uc color(0, 0, 0);
    for (std::size_t n = 0; n < v_format.size(); ++n)
    {
      PlyVertexElement elem = (PlyVertexElement)v_format[n];
      if (elem == PLY_V_FLOAT_X || elem == PLY_V_FLOAT_Y
          || elem == PLY_V_FLOAT_Z)
      {
        float f;
        if (ply_format == PLY_ASCII)
          input >> f;
        else if (ply_format == PLY_BINARY_LE || ply_format == PLY_BINARY_BE)
          input.read((char*)&f, sizeof(float));

        if (ply_format == PLY_BINARY_BE)
          Helpers::swap_endianess(f);

        vertex[(int)elem - PLY_V_FLOAT_X] = f;
      }
      else if (elem == PLY_V_UCHAR_R || elem == PLY_V_UCHAR_G
          || elem == PLY_V_UCHAR_B)
      {
        unsigned char f;
        if (ply_format == PLY_ASCII)
        {
          int tmp_f;
          input >> tmp_f;
          f = (unsigned char)tmp_f;
        }
        else if (ply_format == PLY_BINARY_LE || ply_format == PLY_BINARY_BE)
          input.read((char*)&f, sizeof(unsigned char));

        color[(int)elem - PLY_V_UCHAR_R] = f;
      }
      else if (elem == PLY_V_FLOAT_IGNORE)
      {
        float f;
        if (ply_format == PLY_ASCII)
          input >> f;
        else if (ply_format == PLY_BINARY_LE || ply_format == PLY_BINARY_BE)
          input.read((char*)&f, sizeof(float));
      }
      else if (elem == PLY_V_INT_IGNORE)
      {
        int f;
        if (ply_format == PLY_ASCII)
          input >> f;
        else if (ply_format == PLY_BINARY_LE || ply_format == PLY_BINARY_BE)
          input.read((char*)&f, sizeof(int));
      }
      else /* if (elem == PLY_V_BYTE_IGNORE) */
      {
        char c;
        if (ply_format == PLY_ASCII)
        {
          int tmp_c;
          input >> tmp_c;
          c = (char)tmp_c;
        }
        else if (ply_format == PLY_BINARY_LE || ply_format == PLY_BINARY_BE)
          input.read((char*)&c, sizeof(char));
      }
    }
    vertices.push_back(vertex);
    if (want_colors)
      vcolors.push_back(color);

    /*
    if (i < 10 || i > num_vertices - 10)
      std::cout << "Read: " << vertex << " / ("
          << (int)color[0] << "," << (int)color[1]
          << "," << (int)color[2] << ")" << std::endl;
    */
  }

  /* Start reading the face data. */
  faces.reserve(num_faces * 3);
  for (std::size_t i = 0; i < num_faces; ++i)
  {
    for (std::size_t n = 0; n < f_format.size(); ++n)
    {
      PlyFaceElement elem = (PlyFaceElement)f_format[n];
      if (elem == PLY_F_VERTEX_INDICES)
      {
        unsigned char n_verts;
        if (ply_format == PLY_ASCII)
        {
          input >> n_verts;
          n_verts = (unsigned char)(n_verts - '0');
        }
        else if (ply_format == PLY_BINARY_LE || ply_format == PLY_BINARY_BE)
          input.read((char*)&n_verts, sizeof(unsigned char));

        if (n_verts == 3)
        {
          unsigned int a, b, c;
          if (ply_format == PLY_ASCII)
            input >> a >> b >> c;
          else if (ply_format == PLY_BINARY_LE || ply_format == PLY_BINARY_BE)
          {
            input.read((char*)&a, sizeof(unsigned int));
            input.read((char*)&b, sizeof(unsigned int));
            input.read((char*)&c, sizeof(unsigned int));
          }

          if (ply_format == PLY_BINARY_BE)
          {
            Helpers::swap_endianess(a);
            Helpers::swap_endianess(b);
            Helpers::swap_endianess(c);
          }

          faces.push_back((MeshVIndex)a);
          faces.push_back((MeshVIndex)b);
          faces.push_back((MeshVIndex)c);
        }
        else if (n_verts == 4)
        {
          unsigned int a, b, c, d;
          if (ply_format == PLY_ASCII)
            input >> a >> b >> c >> d;
          else if (ply_format == PLY_BINARY_LE || ply_format == PLY_BINARY_BE)
          {
            input.read((char*)&a, sizeof(unsigned int));
            input.read((char*)&b, sizeof(unsigned int));
            input.read((char*)&c, sizeof(unsigned int));
            input.read((char*)&d, sizeof(unsigned int));
          }

          if (ply_format == PLY_BINARY_BE)
          {
            Helpers::swap_endianess(a);
            Helpers::swap_endianess(b);
            Helpers::swap_endianess(c);
            Helpers::swap_endianess(d);
          }

          faces.push_back((MeshVIndex)a);
          faces.push_back((MeshVIndex)b);
          faces.push_back((MeshVIndex)c);

          faces.push_back((MeshVIndex)c);
          faces.push_back((MeshVIndex)d);
          faces.push_back((MeshVIndex)a);
        }
        else
        {
          std::cout << "PLY Loader: Ignoring face with "
              << n_verts << " vertices!" << std::endl;
        }
      }
      else if (elem == PLY_F_INT_IGNORE)
      {
        int f;
        if (ply_format == PLY_ASCII)
          input >> f;
        else if (ply_format == PLY_BINARY_LE || ply_format == PLY_BINARY_BE)
          input.read((char*)&f, sizeof(int));
      }
      else /* if (elem == PLY_F_BYTE_IGNORE) */
      {
        char c;
        if (ply_format == PLY_ASCII)
          input >> c;
        else if (ply_format == PLY_BINARY_LE || ply_format == PLY_BINARY_BE)
          input.read((char*)&c, sizeof(char));
      }
    }
  }

  /* Close the file stream. */
  input.close();

  //mesh->ensure_normals();
  //mesh->scale_and_center();

  return mesh;
}

/* ---------------------------------------------------------------- */

TriangleMeshPtr
ModelLoader::load_obj_model (std::string const& filename)
{
  /* Open file. */
  std::fstream input(filename.c_str());
  if (input.fail())
    throw Exception("Cannot open file: " + std::string(::strerror(errno)));

  TriangleMeshPtr mesh = TriangleMesh::create();
  MeshVertexList& vertices = mesh->get_vertices();
  MeshFaceList& faces = mesh->get_faces();

  /* Start parsing. The OBJ file format does not have a unique header. */
  while (input.good())
  {
    std::string buffer;
    std::getline(input, buffer);

    Helpers::chop_string(buffer);

    if (!buffer.empty() && buffer[0] == '#')
    {
      std::cout << "OBJ comment: " << buffer << std::endl;
      continue;
    }

    if (buffer.empty())
      continue;

    StringVector sv = Helpers::split_string(buffer, ' ');
    if (sv[0] == "v")
    {
      /* Read a vertex. */
      if (sv.size() != 4)
      {
        std::cout << "OBJ loader: Vertex format not recognized" << std::endl;
        continue;
      }

      Vec3f v(0.0f, 0.0f, 0.0f);
      v[0] = Helpers::get_float_from_string(sv[1]);
      v[1] = Helpers::get_float_from_string(sv[2]);
      v[2] = Helpers::get_float_from_string(sv[3]);
      vertices.push_back(v);

      //std::cout << "Adding vertex: " << vertices.back() << std::endl;
    }
    else if (sv[0] == "f")
    {
      /* Read a face. */
      if (sv.size() != 4)
      {
        std::cout << "OBJ loader: Face format not recognized" << std::endl;
        continue;
      }

      //std::cout << "Adding face " << buffer << std::endl;
      for (int i = 0; i < 3; ++i)
      {
        int vidx = Helpers::get_int_from_string(sv[1 + i]);
        if (vidx < 0)
          vidx = (int)vertices.size() + vidx;
        else
          vidx = vidx - 1;

        faces.push_back(vidx);
      }
    }
    else if (sv[0].empty())
    {
      continue;
    }
    else
    {
      std::cout << "OBJ loader: Element \"" << sv[0]
          << "\" not recognized" << std::endl;
      continue;
    }
  }

  /* Close the file stream. */
  input.close();

  //mesh->ensure_normals();
  //mesh->scale_and_center();

  #if 0
  std::cout << "Mesh vertices: " << vertices.size()
    << ", faces: " << faces.size() / 3 << std::endl;
  #endif

  return mesh;
}

/* ---------------------------------------------------------------- */

TriangleMeshPtr
ModelLoader::load_m_model (std::string const& filename)
{
  /* Open file. */
  std::fstream input(filename.c_str());
  if (input.fail())
    throw Exception("Cannot open file: " + std::string(::strerror(errno)));

  TriangleMeshPtr mesh = TriangleMesh::create();
  MeshVertexList& vertices = mesh->get_vertices();
  MeshFaceList& faces = mesh->get_faces();

  std::vector<int> vidx;

  /* Start parsing. The OBJ file format does not have a unique header. */
  std::size_t line_cnt = 0;
  while (input.good())
  {
    std::string buffer;
    std::getline(input, buffer);
    Helpers::chop_string(buffer);
    Helpers::normalize_string(buffer);
    line_cnt += 1;

    if (buffer.empty())
      continue;

    if (buffer[0] == '#')
      continue;

    StringVector toks = Helpers::split_string(buffer, ' ');

    if (!toks.empty() && toks[0] == "Edge")
        continue;

    if (toks.size() != 5)
    {
      input.close();
      throw Exception("Invalid line format at line "
          + Helpers::get_string(line_cnt));
    }

    if (toks[0] == "Vertex")
    {
      int id = Helpers::get_int_from_string(toks[1]);
      float x = Helpers::get_float_from_string(toks[2]);
      float y = Helpers::get_float_from_string(toks[3]);
      float z = Helpers::get_float_from_string(toks[4]);
      if (vidx.size() <= (std::size_t)id)
        vidx.resize(id + 1, 0);
      vidx[id] = (int)vertices.size();
      vertices.push_back(Vec3f(x, y, z));
    }
    else if (toks[0] == "Face")
    {
      int v1 = Helpers::get_int_from_string(toks[2]);
      int v2 = Helpers::get_int_from_string(toks[3]);
      int v3 = Helpers::get_int_from_string(toks[4]);
      faces.push_back(vidx[v1]);
      faces.push_back(vidx[v2]);
      faces.push_back(vidx[v3]);
    }
    else
    {
      input.close();
      throw Exception("Invalid line identifier");
    }
  }

  input.close();
  return mesh;
}

REMESHER_NAMESPACE_END
