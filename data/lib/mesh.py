#!/usr/bin/env python
"""
Mesh utilities.
"""

import math
import numpy as np
import utils

class Mesh:
    """
    Represents a mesh.
    """

    def __init__(self, vertices = [[]], faces = [[]]):
        """
        Construct a mesh from vertices and faces.

        :param vertices: list of vertices, or numpy array
        :type vertices: [[float]] or numpy.ndarray
        :param faces: list of faces or numpy array, i.e. the indices of the corresponding vertices per triangular face
        :type faces: [[int]] fo rnumpy.ndarray
        """

        self.vertices = np.array(vertices, dtype = float)
        """ (numpy.ndarray) Vertices. """

        self.faces = np.array(faces, dtype = int)
        """ (numpy.ndarray) Faces. """

        assert self.vertices.shape[1] == 3
        assert self.faces.shape[1] == 3

    def extents(self):
        """
        Get the extents.

        :return: (min_x, min_y, min_z), (max_x, max_y, max_z)
        :rtype: (float, float, float), (float, float, float)
        """

        min = [0]*3
        max = [0]*3

        for i in range(3):
            min[i] = np.min(self.vertices[:, i])
            max[i] = np.max(self.vertices[:, i])

        return tuple(min), tuple(max)

    def switch_axes(self, axis_1, axis_2):
        """
        Switch the two axes, this is usually useful for switching y and z axes.

        :param axis_1: index of first axis
        :type axis_1: int
        :param axis_2: index of second axis
        :type axis_2: int
        """

        temp = np.copy(self.vertices[:, axis_1])
        self.vertices[:, axis_1] = self.vertices[:, axis_2]
        self.vertices[:, axis_2] = temp

    def mirror(self, axis):
        """
        Mirror given axis.

        :param axis: axis to mirror
        :type axis: int
        """

        self.vertices[:, axis] *= -1

    def scale(self, scales):
        """
        Scale the mesh in all dimensions.

        :param scales: tuple of length 3 with scale for (x, y, z)
        :type scales: (float, float, float)
        """

        assert len(scales) == 3

        for i in range(3):
            self.vertices[:, i] *= scales[i]

    def translate(self, translation):
        """
        Translate the mesh.

        :param translation: translation as (x, y, z)
        :type translation: (float, float, float)
        """

        assert len(translation) == 3

        for i in range(3):
            self.vertices[:, i] += translation[i]

    def _rotate(self, R):

        self.vertices = np.dot(R, self.vertices.T)
        self.vertices = self.vertices.T

    def rotate(self, rotation):
        """
        Rotate the mesh.

        :param rotation: rotation in (angle_x, angle_y, angle_z); angles in radians
        :type rotation: (float, float, float
        :return:
        """
        assert len(rotation) == 3

        x = rotation[0]
        y = rotation[1]
        z = rotation[2]

        # rotation around the x axis
        R = np.array([[1, 0, 0], [0, math.cos(x), -math.sin(x)], [0, math.sin(x), math.cos(x)]])
        self._rotate(R)

        # rotation around the y axis
        R = np.array([[math.cos(y), 0, math.sin(y)], [0, 1, 0], [-math.sin(y), 0, math.cos(y)]])
        self._rotate(R)

        # rotation around the z axis
        R = np.array([[math.cos(z), -math.sin(z), 0], [math.sin(z), math.cos(z), 0], [0, 0, 1]])
        self._rotate(R)

    def inv_rotate(self, rotation):
        """
        Rotate the mesh.

        :param rotation: rotation in (angle_x, angle_y, angle_z); angles in radians
        :type rotation: (float, float, float
        :return:
        """
        assert len(rotation) == 3

        x = rotation[0]
        y = rotation[1]
        z = rotation[2]

        # rotation around the x axis
        R = np.array([[1, 0, 0], [0, math.cos(x), -math.sin(x)], [0, math.sin(x), math.cos(x)]])
        R = R.T
        self._rotate(R)

        # rotation around the y axis
        R = np.array([[math.cos(y), 0, math.sin(y)], [0, 1, 0], [-math.sin(y), 0, math.cos(y)]])
        R = R.T
        self._rotate(R)

        # rotation around the z axis
        R = np.array([[math.cos(z), -math.sin(z), 0], [math.sin(z), math.cos(z), 0], [0, 0, 1]])
        R = R.T
        self._rotate(R)

    def copy(self):
        """
        Copy the mesh.

        :return: copy of the mesh
        :rtype: Mesh
        """

        mesh = Mesh(self.vertices.copy(), self.faces.copy())

        return mesh

    @staticmethod
    def from_off(filepath):
        """
        Read a mesh from OFF.

        :param filepath: path to OFF file
        :type filepath: str
        :return: mesh
        :rtype: Mesh
        """

        vertices, faces = utils.read_off(filepath)

        real_faces = []
        for face in faces:
            assert len(face) == 4
            real_faces.append([face[1], face[2], face[3]])

        return Mesh(vertices, real_faces)

    def to_off(self, filepath):
        """
        Write mesh to OFF.

        :param filepath: path to write file to
        :type filepath: str
        """

        faces = np.ones((self.faces.shape[0], 4), dtype = int)*3
        faces[:, 1:4] = self.faces[:, :]

        utils.write_off(filepath, self.vertices.tolist(), faces.tolist())

    @staticmethod
    def from_volume(volume):
        """
        Create a mesh from a voxel grid/volume.

        :param volume: volume
        :type volume: numpy.ndarray
        :return: mesh
        :rtype: Mesh
        """

        height = volume.shape[0]
        width = volume.shape[1]
        depth = volume.shape[2]

        vertices = []
        faces = []

        xx, yy, zz = np.where(volume > 0.5)

        for i in range(len(xx)):
            v000 = (yy[i], xx[i], zz[i])                # 0
            v001 = (yy[i], xx[i], zz[i] + 1)            # 1
            v010 = (yy[i], xx[i] + 1, zz[i])            # 2
            v011 = (yy[i], xx[i] + 1, zz[i] + 1)        # 3
            v100 = (yy[i] + 1, xx[i], zz[i])            # 4
            v101 = (yy[i] + 1, xx[i], zz[i] + 1)        # 5
            v110 = (yy[i] + 1, xx[i] + 1, zz[i])        # 6
            v111 = (yy[i] + 1, xx[i] + 1, zz[i] + 1)    # 7

            n = len(vertices)
            f1 = [n + 0, n + 2, n + 4]
            f2 = [n + 4, n + 2, n + 6]
            f3 = [n + 1, n + 3, n + 5]
            f4 = [n + 5, n + 3, n + 7]
            f5 = [n + 0, n + 1, n + 2]
            f6 = [n + 1, n + 3, n + 2]
            f7 = [n + 4, n + 5, n + 7]
            f8 = [n + 4, n + 7, n + 6]
            f9 = [n + 4, n + 0, n + 1]
            f10 = [n + 4, n + 5, n + 1]
            f11 = [n + 2, n + 3, n + 6]
            f12 = [n + 3, n + 7, n + 6]

            vertices.append(v000)
            vertices.append(v001)
            vertices.append(v010)
            vertices.append(v011)
            vertices.append(v100)
            vertices.append(v101)
            vertices.append(v110)
            vertices.append(v111)

            faces.append(f1)
            faces.append(f2)
            faces.append(f3)
            faces.append(f4)
            faces.append(f5)
            faces.append(f6)
            faces.append(f7)
            faces.append(f8)
            faces.append(f9)
            faces.append(f10)
            faces.append(f11)
            faces.append(f12)

        return Mesh(vertices, faces)

    @staticmethod
    def from_obj(filepath):
        """
        Read a mesh from OBJ.

        :param filepath: path to OFF file
        :type filepath: str
        :return: mesh
        :rtype: Mesh
        """

        vertices, faces = utils.read_obj(filepath)
        return Mesh(vertices, faces)

    def to_obj(self, filepath):
        """
        Write mesh to OBJ file.

        :param filepath: path to OBJ file
        :type filepath: str
        """

        utils.write_obj(filepath, self.vertices.tolist(), self.faces.tolist())

    def to_ply(self, filepath):
        """
        Write mesh to OBJ file.

        :param filepath: path to OBJ file
        :type filepath: str
        """

        utils.write_ply(filepath, self.vertices.tolist(), self.faces.tolist())