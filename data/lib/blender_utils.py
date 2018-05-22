# befor eimporting anything!
import sys
import os
sys.path.insert(1, os.path.realpath('/work/dev-box/blender-2.79-linux-glibc219-x86_64/2.79/python/lib/python3.5/site-packages/'))

from enum import Enum

class LogColor:
    INFO = '\033[94m'
    WARNING = '\033[93m'
    ERROR = '\033[91m\033[1m'
    ENDC = '\033[0m'

class LogLevel(Enum):
    INFO = 1
    WARNING = 2
    ERROR = 3

def log(output,level=LogLevel.INFO):
    if level == LogLevel.INFO:
        sys.stderr.write(LogColor.INFO)
    elif level == LogLevel.WARNING:
        sys.stderr.write(LogColor.WARNING)
    elif level == LogLevel.ERROR:
        sys.stderr.write(LogColor.ERROR)
    sys.stderr.write(str(output))
    sys.stderr.write(LogColor.ENDC)
    sys.stderr.write("\n")
    sys.stderr.flush()

import bpy
import bmesh
import math
import numpy as np

import binvox_rw
import import_off
import_off.register()

from bpy_extras.io_utils import axis_conversion
from bpy.props import EnumProperty

sphere_base_mesh = None
cube_base_mesh = None
circle_base_mesh = None


def initialize(width=512, height=448):

    bpy.ops.mesh.primitive_ico_sphere_add()
    global sphere_base_mesh
    sphere_base_mesh = bpy.context.scene.objects.active.data.copy()
    for face in sphere_base_mesh.polygons:
        face.use_smooth = True

    bpy.ops.mesh.primitive_cube_add()
    global cube_base_mesh
    cube_base_mesh = bpy.context.scene.objects.active.data.copy()

    bpy.ops.mesh.primitive_circle_add(vertices=1024, radius=1, fill_type='NGON')
    global circle_base_mesh
    circle_base_mesh = bpy.context.scene.objects.active.data.copy()

    # Delete the scene, except for the camera and the lamp
    for obj in bpy.data.objects:
        if str(obj.name) in ['Camera']:
            continue
        obj.select = True
        bpy.ops.object.delete()

    scene = bpy.context.scene

    # set the camera and its constraint
    cam = scene.objects['Camera']
    cam.location = (0, 3.0, 1.0)
    cam.data.lens = 35
    cam.data.sensor_width = 32
    cam.data.sensor_height = 32
    cam_constraint = cam.constraints.new(type='TRACK_TO')
    cam_constraint.track_axis = 'TRACK_NEGATIVE_Z'
    cam_constraint.up_axis = 'UP_Y'

    def parent_obj_to_camera(b_camera):
        origin = (0, 0, 0)
        b_empty = bpy.data.objects.new('Empty', None)
        b_empty.location = origin
        b_camera.parent = b_empty  # setup parenting

        scn = bpy.context.scene
        scn.objects.link(b_empty)
        scn.objects.active = b_empty
        return b_empty

    camera_target = parent_obj_to_camera(cam)
    cam_constraint.target = camera_target

    locations = [
        (-0.98382, 0.445997, 0.526505),
        (-0.421806, -0.870784, 0.524944),
        (0.075576, -0.960128, 0.816464),
        (0.493553, -0.57716, 0.928208),
        (0.787275, -0.256822, 0.635172),
        (1.01032, 0.148764, 0.335078)
    ]

    for i in range(len(locations)):
        lamp_data = bpy.data.lamps.new(name='Point Lamp ' + str(i), type='POINT')
        lamp_data.shadow_method = 'RAY_SHADOW'
        lamp_data.shadow_ray_sample_method = 'CONSTANT_QMC'
        lamp_data.use_shadow = True
        lamp_data.shadow_soft_size = 1e6
        lamp_data.distance = 2
        lamp_data.energy = 0.1
        lamp_data.use_diffuse = True
        lamp_data.use_specular = True
        lamp_data.falloff_type = 'CONSTANT'

        lamp_object = bpy.data.objects.new(name='Spot Lamp ' + str(i), object_data=lamp_data)
        scene.objects.link(lamp_object)
        lamp_object.location[0] = locations[i][0]
        lamp_object.location[1] = locations[i][1]
        lamp_object.location[2] = locations[i][2]
        lamp_object.rotation_euler[0] = 0
        lamp_object.rotation_euler[1] = 0
        lamp_object.rotation_euler[2] = 0
        lamp_object.parent = camera_target

    try:
        if (2, 78, 0) <= bpy.app.version:
            # https://blender.stackexchange.com/questions/5281/blender-sets-compute-device-cuda-but-doesnt-use-it-for-actual-render-on-ec2
            bpy.context.user_preferences.addons['cycles'].preferences.compute_device_type = 'CUDA'
            bpy.context.user_preferences.addons['cycles'].preferences.devices[0].use = True
        else:
            bpy.context.user_preferences.system.compute_device_type = 'CUDA'
    except TypeError:
        pass

    scene.render.use_file_extension = False
    scene.render.resolution_x = width
    scene.render.resolution_y = height
    scene.render.resolution_percentage = 100
    scene.render.use_antialiasing = True
    scene.render.use_shadows = True
    world = bpy.context.scene.world
    world.zenith_color = [1.0, 1.0, 1.0]
    world.horizon_color = [1.0, 1.0, 1.0]
    scene.render.alpha_mode = 'SKY'
    world.light_settings.use_environment_light = True
    world.light_settings.environment_color = 'PLAIN'
    world.light_settings.environment_energy = 0.5

    return camera_target


def make_material(name, diffuse, alpha, shadow=False):
    material = bpy.data.materials.new(name)
    material.diffuse_color = diffuse
    material.diffuse_shader = 'LAMBERT'
    material.diffuse_intensity = 1
    material.specular_color = (1, 1, 1)
    material.specular_shader = 'COOKTORR'
    material.specular_intensity = 2
    material.alpha = alpha
    material.use_transparency = True
    material.ambient = 1.0

    material.use_cast_shadows = shadow
    material.use_shadows = shadow

    return material


def shadow_plane(material, offset = (0, 0, 0), scale = 1):

    global circle_base_mesh
    ob = bpy.data.objects.new("BRC_Shadow_Plane", circle_base_mesh)

    ob.location = offset
    ob.scale = (scale, scale, scale)
    bpy.context.scene.objects.link(ob)

    mat = material
    mat.use_shadows = True
    mat.use_transparent_shadows = True
    mat.use_only_shadow = True
    mat.use_raytrace = True
    mat.ambient = 0

    ob.data.materials.append(mat)
    ob.active_material_index = 0
    ob.active_material = mat


def _load_mesh(name, vertices, faces):

    # vertices should be list of lists
    # faces should be list of lists
    edges = []
    mesh = bpy.data.meshes.new(name=name)
    mesh.from_pydata(vertices, edges, faces)
    # mesh.vertices.add(len(verts))
    # mesh.vertices.foreach_set("co", unpack_list(verts))

    # mesh.faces.add(len(facets))
    # mesh.faces.foreach_set("vertices", unpack_face_list(facets))

    mesh.validate()
    mesh.update()

    scene = bpy.context.scene
    obj = bpy.data.objects.new(mesh.name, mesh)
    scene.objects.link(obj)
    scene.objects.active = obj
    obj.select = True

    axis_forward = EnumProperty(
        name="Forward",
        items=(('X', "X Forward", ""),
               ('Y', "Y Forward", ""),
               ('Z', "Z Forward", ""),
               ('-X', "-X Forward", ""),
               ('-Y', "-Y Forward", ""),
               ('-Z', "-Z Forward", ""),
               ),
        default='Y',
    )
    axis_up = EnumProperty(
        name="Up",
        items=(('X', "X Up", ""),
               ('Y', "Y Up", ""),
               ('Z', "Z Up", ""),
               ('-X', "-X Up", ""),
               ('-Y', "-Y Up", ""),
               ('-Z', "-Z Up", ""),
               ),
        default='Z',
    )

    global_matrix = axis_conversion(from_forward=axis_forward, from_up=axis_up).to_4x4()

    obj.matrix_world = global_matrix
    scene.update()

    return mesh


def load_mesh(name, vertices, faces, material, offset=(0, 0, 0), scale=1, axes='xyz'):
    _load_mesh(name, vertices, faces)

    assert len(offset) == 3
    assert scale > 0
    assert len(axes) == 3

    x_index = axes.find('x')
    y_index = axes.find('y')
    z_index = axes.find('z')

    assert x_index >= 0 and x_index < 3
    assert y_index >= 0 and y_index < 3
    assert z_index >= 0 and z_index < 3
    assert x_index != y_index and x_index != z_index and y_index != z_index

    for obj in bpy.context.scene.objects:

        # obj.name contains the group name of a group of faces, see http://paulbourke.net/dataformats/obj/
        # every mesh is of type 'MESH', this works not only for ShapeNet but also for 'simple'
        # obj files
        if obj.type == 'MESH' and not 'BRC' in obj.name:

            # change color
            # this is based on https://stackoverflow.com/questions/4644650/blender-how-do-i-add-a-color-to-an-object
            # but needed changing a lot of attributes according to documentation
            obj.data.materials.append(material)

            for vertex in obj.data.vertices:
                # make a copy, otherwise axes switching does not work
                vertex_copy = (vertex.co[0], vertex.co[1], vertex.co[2])

                vertex.co[0] = vertex_copy[x_index]
                vertex.co[1] = vertex_copy[y_index]
                vertex.co[2] = vertex_copy[z_index]

                vertex.co[0] = vertex.co[0] * scale + offset[0]
                vertex.co[1] = vertex.co[1] * scale + offset[1]
                vertex.co[2] = vertex.co[2] * scale + offset[2]

            obj.name = 'BRC_' + obj.name


def load_off(off_file, material, offset=(0, 0, 0), scale=1, axes='xyz'):
    bpy.ops.import_mesh.off(filepath=off_file)

    assert len(offset) == 3
    assert scale > 0
    assert len(axes) == 3

    x_index = axes.find('x')
    y_index = axes.find('y')
    z_index = axes.find('z')

    assert x_index >= 0 and x_index < 3
    assert y_index >= 0 and y_index < 3
    assert z_index >= 0 and z_index < 3
    assert x_index != y_index and x_index != z_index and y_index != z_index

    for obj in bpy.context.scene.objects:

        # obj.name contains the group name of a group of faces, see http://paulbourke.net/dataformats/obj/
        # every mesh is of type 'MESH', this works not only for ShapeNet but also for 'simple'
        # obj files
        if obj.type == 'MESH' and not 'BRC' in obj.name:

            # change color
            # this is based on https://stackoverflow.com/questions/4644650/blender-how-do-i-add-a-color-to-an-object
            # but needed changing a lot of attributes according to documentation
            obj.data.materials.append(material)

            for vertex in obj.data.vertices:
                # make a copy, otherwise axes switching does not work
                vertex_copy = (vertex.co[0], vertex.co[1], vertex.co[2])

                vertex.co[0] = vertex_copy[x_index]
                vertex.co[1] = vertex_copy[y_index]
                vertex.co[2] = vertex_copy[z_index]

                vertex.co[0] = vertex.co[0] * scale + offset[0]
                vertex.co[1] = vertex.co[1] * scale + offset[1]
                vertex.co[2] = vertex.co[2] * scale + offset[2]

            obj.name = 'BRC_' + obj.name


def load_txt(txt_file, radius, material, offset=(0, 0, 0), scale=1, axes='xyz'):
    global sphere_base_mesh

    assert len(offset) == 3
    assert scale > 0
    assert len(axes) == 3

    x_index = axes.find('x')
    y_index = axes.find('y')
    z_index = axes.find('z')

    assert x_index >= 0 and x_index < 3
    assert y_index >= 0 and y_index < 3
    assert z_index >= 0 and z_index < 3
    assert x_index != y_index and x_index != z_index and y_index != z_index

    voxel_file = open(txt_file, 'r')
    voxel_lines = voxel_file.readlines()
    voxel_file.close()

    mesh = bmesh.new()
    for line in voxel_lines:
        vals = line.split(' ')
        if not line.startswith('#') and line.strip() != '' and len(vals) >= 3:
            location = (
                float(vals[x_index]) * scale + offset[0],
                float(vals[y_index]) * scale + offset[1],
                float(vals[z_index]) * scale + offset[2]
            )

            m = sphere_base_mesh.copy()
            for vertex in m.vertices:
                vertex.co[0] = vertex.co[0] * radius + location[0]
                vertex.co[1] = vertex.co[1] * radius + location[1]
                vertex.co[2] = vertex.co[2] * radius + location[2]

            mesh.from_mesh(m)

    mesh2 = bpy.data.meshes.new('Mesh')
    mesh.to_mesh(mesh2)

    obj = bpy.data.objects.new('BRC_Point_Cloud', mesh2)
    obj.data.materials.append(material)
    obj.active_material_index = 0
    obj.active_material = material

    bpy.context.scene.objects.link(obj)


def load_binvox(binvox_file, radius, material, offset, scale, axes):
    global cube_base_mesh

    assert len(offset) == 3
    assert len(scale) == 3
    assert len(axes) == 3

    x_index = axes.find("x")
    y_index = axes.find("y")
    z_index = axes.find("z")

    assert x_index >= 0 and x_index < 3
    assert y_index >= 0 and y_index < 3
    assert z_index >= 0 and z_index < 3
    assert x_index != y_index and x_index != z_index and y_index != z_index

    with open(binvox_file, 'rb') as f:
        model = binvox_rw.read_as_3d_array(f)

    points = np.where(model.data)
    locations = np.zeros((points[0].shape[0], 3), dtype=float)
    locations[:, 0] = (points[x_index][:] + 0.5) / model.data.shape[x_index]
    locations[:, 1] = (points[y_index][:] + 0.5) / model.data.shape[y_index]
    locations[:, 2] = (points[z_index][:] + 0.5) / model.data.shape[z_index]
    locations[:, 0] -= 0.5
    locations[:, 1] -= 0.5
    locations[:, 2] -= 0.5

    locations[:, 0] = locations[:, 0] * scale[0] + offset[0]
    locations[:, 1] = locations[:, 1] * scale[1] + offset[1]
    locations[:, 2] = locations[:, 2] * scale[2] + offset[2]

    mesh = bmesh.new()
    for i in range(locations.shape[0]):
            m = cube_base_mesh.copy()
            for vertex in m.vertices:
                vertex.co[0] = vertex.co[0] * radius + locations[i, 0]
                vertex.co[1] = vertex.co[1] * radius + locations[i, 1]
                vertex.co[2] = vertex.co[2] * radius + locations[i, 2]

            mesh.from_mesh(m)

    mesh2 = bpy.data.meshes.new('Mesh')
    mesh.to_mesh(mesh2)

    obj = bpy.data.objects.new('BRC_Occupancy', mesh2)
    obj.data.materials.append(material)
    obj.active_material_index = 0
    obj.active_material = material

    bpy.context.scene.objects.link(obj)


def load_volume(volume, radius, material, offset, scale, axes):
    global cube_base_mesh

    assert len(offset) == 3
    assert len(scale) == 3
    assert len(axes) == 3

    x_index = axes.find("x")
    y_index = axes.find("y")
    z_index = axes.find("z")

    assert x_index >= 0 and x_index < 3
    assert y_index >= 0 and y_index < 3
    assert z_index >= 0 and z_index < 3
    assert x_index != y_index and x_index != z_index and y_index != z_index

    points = np.where(volume > 0)
    locations = np.zeros((points[0].shape[0], 3), dtype=float)
    locations[:, 0] = (points[x_index][:] + 0.5) / volume.shape[x_index]
    locations[:, 1] = (points[y_index][:] + 0.5) / volume.shape[y_index]
    locations[:, 2] = (points[z_index][:] + 0.5) / volume.shape[z_index]
    locations[:, 0] -= 0.5
    locations[:, 1] -= 0.5
    locations[:, 2] -= 0.5

    locations[:, 0] = locations[:, 0] * scale[0] + offset[0]
    locations[:, 1] = locations[:, 1] * scale[1] + offset[1]
    locations[:, 2] = locations[:, 2] * scale[2] + offset[2]

    mesh = bmesh.new()
    for i in range(locations.shape[0]):
            m = cube_base_mesh.copy()
            for vertex in m.vertices:
                vertex.co[0] = vertex.co[0] * radius + locations[i, 0]
                vertex.co[1] = vertex.co[1] * radius + locations[i, 1]
                vertex.co[2] = vertex.co[2] * radius + locations[i, 2]

            mesh.from_mesh(m)

    mesh2 = bpy.data.meshes.new('Mesh')
    mesh.to_mesh(mesh2)

    obj = bpy.data.objects.new('BRC_Occupancy', mesh2)
    obj.data.materials.append(material)
    obj.active_material_index = 0
    obj.active_material = material

    bpy.context.scene.objects.link(obj)


def render(camera_target, output_file, rotation, distance):
    bpy.context.scene.render.filepath = output_file

    camera_target.rotation_euler[0] = math.radians(rotation[0])
    camera_target.rotation_euler[1] = math.radians(rotation[1])
    camera_target.rotation_euler[2] = math.radians(rotation[2])

    cam = bpy.context.scene.objects['Camera']
    cam.location = (0, 3.0 * distance, 1.0 * distance)

    bpy.ops.render.render(animation=False, write_still=True)
