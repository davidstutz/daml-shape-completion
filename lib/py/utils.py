#!/usr/bin/env python
"""
Some I/O utilities.
"""

import os
import re
import json
import h5py
import numpy as np
import zipfile

def write_hdf5(file, tensor, key = 'tensor'):
    """
    Write a simple tensor, i.e. numpy array ,to HDF5.

    :param file: path to file to write
    :type file: str
    :param tensor: tensor to write
    :type tensor: numpy.ndarray
    :param key: key to use for tensor
    :type key: str
    """

    assert type(tensor) == np.ndarray, 'file %s not found' % file

    h5f = h5py.File(file, 'w')

    chunks = list(tensor.shape)
    if len(chunks) > 2:
        chunks[2] = 1
        if len(chunks) > 3:
            chunks[3] = 1
            if len(chunks) > 4:
                chunks[4] = 1

    h5f.create_dataset(key, data = tensor, chunks = tuple(chunks), compression = 'gzip')
    h5f.close()

def read_hdf5(file, key = 'tensor'):
    """
    Read a tensor, i.e. numpy array, from HDF5.

    :param file: path to file to read
    :type file: str
    :param key: key to read
    :type key: str
    :return: tensor
    :rtype: numpy.ndarray
    """

    assert os.path.exists(file), 'file %s not found' % file

    h5f = h5py.File(file, 'r')
    tensor = h5f[key][()]
    h5f.close()

    return tensor

def read_json(file):
    """
    Read a JSON file.

    :param file: path to file to read
    :type file: str
    :return: parsed JSON as dict
    :rtype: dict
    """

    assert os.path.exists(file), 'file %s not found' % file

    with open(file, 'r') as fp:
        return json.load(fp)

def write_json(file, data):
    """
    Read a JSON file.

    :param file: path to file to read
    :type file: str
    :param data: data to write
    :type data: mixed
    :return: parsed JSON as dict
    :rtype: dict
    """

    with open(file, 'w') as fp:
        json.dump(data, fp)

def write_off(file, vertices, faces):
    """
    Writes the given vertices and faces to OFF.

    :param vertices: vertices as tuples of (x, y, z) coordinates
    :type vertices: [(float)]
    :param faces: faces as tuples of (num_vertices, vertex_id_1, vertex_id_2, ...)
    :type faces: [(int)]
    """

    num_vertices = len(vertices)
    num_faces = len(faces)

    assert num_vertices > 0
    assert num_faces > 0

    with open(file, 'w') as fp:
        fp.write('OFF\n')
        fp.write(str(num_vertices) + ' ' + str(num_faces) + ' 0\n')

        for vertex in vertices:
            assert len(vertex) == 3, 'invalid vertex with %d dimensions found (%s)' % (len(vertex), file)
            fp.write(str(vertex[0]) + ' ' + str(vertex[1]) + ' ' + str(vertex[2]) + '\n')

        for face in faces:
            assert face[0] == 3, 'only triangular faces supported (%s)' % file
            assert len(face) == 4, 'faces need to have 3 vertices, but found %d (%s)' % (len(face), file)

            for i in range(len(face)):
                assert face[i] >= 0 and face[i] < num_vertices, 'invalid vertex index %d (of %d vertices) (%s)' % (face[i], num_vertices, file)

                fp.write(str(face[i]))
                if i < len(face) - 1:
                    fp.write(' ')

            fp.write('\n')

        # add empty line to be sure
        fp.write('\n')

def read_off(file):
    """
    Reads vertices and faces from an off file.

    :param file: path to file to read
    :type file: str
    :return: vertices and faces as lists of tuples
    :rtype: [(float)], [(int)]
    """

    assert os.path.exists(file), 'file %s not found' % file

    with open(file, 'r') as fp:
        lines = fp.readlines()
        lines = [line.strip() for line in lines]

        assert lines[0] == 'OFF' or lines[0] == 'off', 'invalid OFF file %s' % file

        parts = lines[1].split(' ')
        assert len(parts) == 3

        num_vertices = int(parts[0])
        assert num_vertices > 0

        num_faces = int(parts[1])
        assert num_faces > 0

        vertices = []
        for i in range(num_vertices):
            vertex = lines[2 + i].split(' ')
            vertex = [float(point) for point in vertex]
            assert len(vertex) == 3

            vertices.append(vertex)

        faces = []
        for i in range(num_faces):
            face = lines[2 + num_vertices + i].split(' ')
            face = [index.strip() for index in face]

            # check to be sure
            for index in face:
                assert index != '', 'found empty vertex index: %s (%s)' % (lines[2 + num_vertices + i], file)

            face = [int(index) for index in face]

            assert face[0] == len(face) - 1, 'face should have %d vertices but as %d (%s)' % (face[0], len(face) - 1, file)
            assert face[0] == 3, 'only triangular meshes supported (%s)' % file
            for index in face:
                assert index >= 0 and index < num_vertices, 'vertex %d (of %d vertices) does not exist (%s)' % (index, num_vertices, file)

            assert len(face) > 1

            faces.append(face)

        return vertices, faces

    assert False, 'could not open %s' % file

def write_obj(file, vertices, faces):
    """
    Writes the given vertices and faces to OBJ.

    :param vertices: vertices as tuples of (x, y, z) coordinates
    :type vertices: [(float)]
    :param faces: faces as tuples of (num_vertices, vertex_id_1, vertex_id_2, ...)
    :type faces: [(int)]
    """

    num_vertices = len(vertices)
    num_faces = len(faces)

    assert num_vertices > 0
    assert num_faces > 0

    with open(file, 'w') as fp:
        for vertex in vertices:
            assert len(vertex) == 3, 'invalid vertex with %d dimensions found (%s)' % (len(vertex), file)
            fp.write('v' + ' ' + str(vertex[0]) + ' ' + str(vertex[1]) + ' ' + str(vertex[2]) + '\n')

        for face in faces:
            assert len(face) == 3, 'only triangular faces supported (%s)' % file
            fp.write('f ')

            for i in range(len(face)):
                assert face[i] >= 0 and face[i] < num_vertices, 'invalid vertex index %d (of %d vertices) (%s)' % (face[i], num_vertices, file)

                # face indices are 1-based
                fp.write(str(face[i] + 1))
                if i < len(face) - 1:
                    fp.write(' ')

            fp.write('\n')

        # add empty line to be sure
        fp.write('\n')

def read_obj(file):
    """
    Reads vertices and faces from an obj file.

    :param file: path to file to read
    :type file: str
    :return: vertices and faces as lists of tuples
    :rtype: [(float)], [(int)]
    """

    assert os.path.exists(file), 'file %s not found' % file

    with open(file, 'r') as fp:
        lines = fp.readlines()
        lines = [line.strip() for line in lines if line.strip()]

        vertices = []
        faces = []
        for line in lines:
            parts = line.split(' ')
            parts = [part.strip() for part in parts if part]

            if parts[0] == 'v':
                assert len(parts) == 4, \
                    'vertex should be of the form v x y z, but found %d parts instead (%s)' % (len(parts), file)
                assert parts[1] != '', 'vertex x coordinate is empty (%s)' % file
                assert parts[2] != '', 'vertex y coordinate is empty (%s)' % file
                assert parts[3] != '', 'vertex z coordinate is empty (%s)' % file

                vertices.append([float(parts[1]), float(parts[2]), float(parts[3])])
            elif parts[0] == 'f':
                assert len(parts) == 4, \
                    'face should be of the form f v1/vt1/vn1 v2/vt2/vn2 v2/vt2/vn2, but found %d parts (%s) instead (%s)' % (len(parts), line, file)

                components = parts[1].split('/')
                assert len(components) >= 1 and len(components) <= 3, \
                   'face component should have the forms v, v/vt or v/vt/vn, but found %d components instead (%s)' % (len(components), file)
                assert components[0].strip() != '', \
                    'face component is empty (%s)' % file
                v1 = int(components[0])

                components = parts[2].split('/')
                assert len(components) >= 1 and len(components) <= 3, \
                    'face component should have the forms v, v/vt or v/vt/vn, but found %d components instead (%s)' % (len(components), file)
                assert components[0].strip() != '', \
                    'face component is empty (%s)' % file
                v2 = int(components[0])

                components = parts[3].split('/')
                assert len(components) >= 1 and len(components) <= 3, \
                    'face component should have the forms v, v/vt or v/vt/vn, but found %d components instead (%s)' % (len(components), file)
                assert components[0].strip() != '', \
                    'face component is empty (%s)' % file
                v3 = int(components[0])

                #assert v1 != v2 and v2 != v3 and v3 != v2, 'degenerate face detected: %d %d %d (%s)' % (v1, v2, v3, file)
                if v1 == v2 or v2 == v3 or v1 == v3:
                    print('[Info] skipping degenerate face in %s' % file)
                else:
                    faces.append([v1 - 1, v2 - 1, v3 - 1]) # indices are 1-based!
            else:
                assert False, 'expected either vertex or face but got line: %s (%s)' % (line, file)

        return vertices, faces

    assert False, 'could not open %s' % file

def write_ply(file, vertices, faces):
    """
    Writes the given vertices and faces to PLY.

    :param vertices: vertices as tuples of (x, y, z) coordinates
    :type vertices: [(float)]
    :param faces: faces as tuples of (num_vertices, vertex_id_1, vertex_id_2, ...)
    :type faces: [(int)]
    """

    num_vertices = len(vertices)
    num_faces = len(faces)

    assert num_vertices > 0
    assert num_faces > 0

    with open(file, 'w') as fp:
        fp.write('ply\n')
        fp.write('format ascii 1.0\n')
        fp.write('element vertex ' + str(len(vertices)))
        fp.write('property float x')
        fp.write('property float y')
        fp.write('property float z')
        fp.write('element face ' + str(len(faces)))
        fp.write('property list uchar int vertex_indices')
        fp.write('end_header')

        for vertex in vertices:
            assert len(vertex) == 3, 'invalid vertex with %d dimensions found (%s)' % (len(vertex), file)
            fp.write(str(vertex[0]) + ' ' + str(vertex[1]) + ' ' + str(vertex[2]) + '\n')

        for face in faces:
            assert len(face) == 3, 'only triangular faces supported (%s)' % file
            fp.write('3 ')

            for i in range(len(face)):
                assert face[i] >= 0 and face[i] < num_vertices, 'invalid vertex index %d (of %d vertices) (%s)' % (face[i], num_vertices, file)

                # face indices are 0-based
                fp.write(str(face[i]))
                if i < len(face) - 1:
                    fp.write(' ')

            fp.write('\n')

        # add empty line to be sure
        fp.write('\n')

def read_ordered_directory(dir):
    """
    Gets a list of file names ordered by integers (if integers are found
    in the file names).

    :param dir: path to directory
    :type dir: str
    :return: list of file names
    :rtype: [str]
    """

    # http://stackoverflow.com/questions/4623446/how-do-you-sort-files-numerically
    def get_int(value):
        """
        Convert the input value to integer if possible.

        :param value: mixed input value
        :type value: mixed
        :return: value as integer, or value
        :rtype: mixed
        """

        try:
            return int(value)
        except:
            return value

    def alphanum_key(string):
        """
        Turn a string into a list of string and number chunks,
        e.g. "z23a" -> ["z", 23, "a"].

        :param string: input string
        :type string: str
        :return: list of elements
        :rtype: [int|str]
        """

        return [get_int(part) for part in re.split('([0-9]+)', string)]

    def sort_filenames(filenames):
        """
        Sort the given list by integers if integers are found in the element strings.

        :param filenames: file names to sort
        :type filenames: [str]
        """

        filenames.sort(key = alphanum_key)

    assert os.path.exists(dir), 'directory %s not found' % dir

    filenames = [dir + '/' + filename for filename in os.listdir(dir)]
    sort_filenames(filenames)

    return filenames

def extract_zip(zip_file, out_dir):
    """
    Extract a ZIP file.

    :param zip_file: path to ZIP file
    :type zip_file: str
    :param out_dir: path to extract ZIP file to
    :type out_dir: str
    """

    zip_ref = zipfile.ZipFile(zip_file, 'r')
    zip_ref.extractall(out_dir)
    zip_ref.close()

def read_point_cloud(file):
    """
    Reads point cloud.

    :param file: path to file to read
    :type file: str
    :return: points
    :rtype: [numpy.ndarray]
    """

    assert os.path.exists(file), 'file %s not found' % file

    with open(file, 'r') as fp:
        lines = fp.readlines()
        lines = [line.strip() for line in lines if line.strip()]

        num_points = int(lines[0])
        lines = lines[1:]

        points = np.zeros((num_points, 3))

        i = 0
        for line in lines:
            parts = line.split(' ')
            assert len(parts) == 3, 'expected 3 points but got %d (%s)' % (len(parts), file)

            points[i, 0] = float(parts[0])
            points[i, 1] = float(parts[1])
            points[i, 2] = float(parts[2])
            i = i + 1

        return points

def makedir(dir):
    """
    Creates directory if it does not exist.

    :param dir: directory path
    :type dir: str
    """

    if not os.path.exists(dir):
        os.makedirs(dir)