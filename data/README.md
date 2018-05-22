# Weakly-Supervised Shape Completion

The data corresponds to

    @inproceedings{Stutz2018CVPR,
        title = {Learning 3D Shape Completion from Laser Scan Data with Weak Supervision },
        author = {Stutz, David and Geiger, Andreas},
        booktitle = {IEEE Conference on Computer Vision and Pattern Recognition (CVPR)},
        publisher = {IEEE Computer Society},
        year = {2018}
    }

Please cite this work, as well as [3] and [4] when using the data.

    [3] Andreas Geiger, Philip Lenz, Raquel Urtasun:
        Are we ready for autonomous driving? The KITTI vision benchmark suite. CVPR 2012: 3354-3361
    [4] Angel X. Chang, Thomas A. Funkhouser, Leonidas J. Guibas, Pat Hanrahan, Qi-Xing Huang, Zimo Li, Silvio Savarese, Manolis Savva, Shuran Song, Hao Su, Jianxiong Xiao, Li Yi, Fisher Yu:
        ShapeNet: An Information-Rich 3D Model Repository. CoRR abs/1512.03012 (2015)

![Examples of the data.](screenshot.jpg?raw=true "Examples of the data.")

## Data

The data is derived from [ShapeNet](https://www.shapenet.org/terms)
and [KITTI](http://www.cvlibs.net/datasets/kitti/). For ShapeNet,
two datasets, in the paper referred to as SN-clean and SN-noisy,
were created. We also provide the simplified cars from ShapeNet, simplified
using [this semi-convex hull algorithm](http://www.cvlibs.net/software/semi_convex_hull/):

* [KITTI (2.7GB)](https://s3.eu-central-1.amazonaws.com/avg-shape-completion/cvpr2018_shape_completion_kitti.zip)
* [SN-clean (5.7GB)](https://s3.eu-central-1.amazonaws.com/avg-shape-completion/cvpr2018_shape_completion_clean.zip)
* [SN-noisy (4.0GB)](https://s3.eu-central-1.amazonaws.com/avg-shape-completion/cvpr2018_shape_completion_noisy.zip)
* [Simplified ShapeNet Cars (38.5MB)](https://s3.eu-central-1.amazonaws.com/avg-shape-completion/cvpr2018_shape_completion_shapenet_models.zip)

**If these links do not work anymore, please let us know!**

Details on the included files:

| File | Description |
| ---  | --- |
| training_inference_off_gt_10_48x64_24x54x24_clean_large | Training meshes as OFF files, scaled to `[0,24] x [0,54] x [0,24]`. |
| training_inference_txt_gt_10_48x64_24x54x24_clean_large | Training point clouds as TXT files, same scale. |
| training_prior_off_gt_10_48x64_24x54x24_clean_large | Meshes as OFF files for training a shape prior if applicable, same scale. |
| validation_off_gt_6_48x64_24x54x24_clean_large | Validation meshes as OFF files, same scale. |
| validation_txt_gt_6_48x64_24x54x24_clean_large | Validation point clouds as OFF files, same scale. |
| training_inference_inputs_10_48x64_24x54x24_clean_large.h5 | Training point cloud occupancy grids, voxelized into `24 x 54 x 24`, resulting in a `N x 1 x 24 x 54 x 24` tensor. |
| training_inference_inputs_lsdf_10_48x64_24x54x24_clean_large.h5 | Training point clouds log distance functions, same size. |
| training_inference_outputs_10_48x64_24x54x24_clean_large.h5 | Training shapes as occupancy grids, voxelized into same size. |
| training_inference_outputs_lsdf_10_48x64_24x54x24_clean_large.h5 | Training shapes as log signed distance functions, same size. |
| training_inference_space_10_48x64_24x54x24_clean_large.h5 | Training free space (i.e. free space derived from point clouds), voxelized into same size. |
| training_prior_outputs_10_48x64_24x54x24_clean_large.h5 | Occupancy grids of shapes for training a shape prior (if applicable), same size. |
| training_prior_outputs_lsdf_10_48x64_24x54x24_clean_large.h5 | Log signed distance functions of shapes for training a shape prior (if applicable), same size. |
| validation_inputs_6_48x64_24x54x24_clean_large.h5 | Validation point cloud occupancy grids, voxelized into same size. |
| validation_inputs_lsdf_6_48x64_24x54x24_clean_large.h5 |Validation point clouds voxelized into log distance functions, same size. |
| validation_outputs_6_48x64_24x54x24_clean_large.h5 | Validation shapes as occupancy grids, same size. |
| validation_outputs_lsdf_6_48x64_24x54x24_clean_large.h5 | Validation shapes as log signed distance functions, same size. |
| validation_space_6_48x64_24x54x24_clean_large.h5 | Validation free space, same size. |
| real_space_statistics_training_prior.h5 | Occupancy probabilities computed on the shape prior training set; indicates the probability of a voxel being occupied. |

| File | Description |
| --- | --- |
| bounding_boxes_txt_training_padding_corrected_1_24x54x24 | Training point clouds as TXT files, scaled to `[0,24] x [0,54] x [0,24]`. |
| bounding_boxes_txt_validation_gt_padding_corrected_1_24x54x24 | Validation point clouds as TXT files, same scale. |
| velodyne_gt_txt_training_padding_corrected_1_24x54x24 | Training ground truth point clouds as TXT files, same scale. |
| velodyne_gt_txt_validation_gt_padding_corrected_1_24x54x24 | Validation ground truth point clouds as TXT files, same scale. |
| bounding_boxes_training_padding_corrected_1_24x54x24.txt | List of training bounding boxes as TXT files, each bounding box being described by `(width_x, height_y, depth_z, translation_x, translation_y, translation_z, rotation_x, rotation_y, rotation_z)`, not scaled. | 
| bounding_boxes_validation_gt_padding_corrected_1_24x54x24.txt | List of validation bounding boxes as above. |
| input_training_padding_corrected_1_24x54x24_f.h5 | Training point clouds as occupancy grids, as `24 x 54 x 24` resulting in a `N x 1 x 24 x 54 x 24` tensor. |
| input_validation_gt_padding_corrected_1_24x54x24_f.h5 | Validation point clouds as occupancy grids, same size. |
| input_lsdf_training_padding_corrected_1_24x54x24_f.h5 | Training point clouds as log distance functions, same size. |
| input_lsdf_validation_gt_padding_corrected_1_24x54x24_f.h5 | Validation point clouds as log distance functions, same size. |
| part_space_training_padding_corrected_1_24x54x24_f.h5 | Training free space as occupancy grids, same size. |
| part_space_validation_gt_padding_corrected_1_24x54x24_f.h5 | Validation free space as occupancy grids, same size. |
| real_space_statistics_training_prior.h5 | Occupancy probabilities computed on the shape prior training set; indicates the probability of a voxel being occupied. |

Also see the [supplementary material](http://davidstutz.de/wordpress/wp-content/uploads/2018/04/shape-completion-cvpr2018-suppmat.pdf);
some important notes:

* For `.h5` files, the shapes/point clouds are stored `height x width x depth`.
* For `.off` and `.txt` files the axes are: x = right, y = up, z = forward. Note that this is
  different from many other tools (including KITTI's tools) use different axes.

## Tools

Some tools for I/O can be found in `tools/`. Tools for visualization can be
found in [davidstutz/bpy-visualization-utils](https://github.com/bpy-visualization-utils).

The tools are mostly self-explanatory, but include:

* Python:
    * Conversion between [OFF](http://shape.cs.princeton.edu/benchmark/documentation/off_format.html)
      and [OBJ](http://paulbourke.net/dataformats/obj/) formats.
    * Conversion between our TXT point cloud format and [PLY](http://paulbourke.net/dataformats/ply/)
      format.
    * I/O code for KITTI's bounding boxes (note that the format is not the
      same as used for KITTI).
* LUA:
    * Example for HDF% I/O.
* C++:
    * Examples for reading OFF, HDF5 and TXT point clouds files.

## Data Generation

The code provided in this directory comes untested and without warranty.
It is provided merely as reference of how the data was derived.

The code should be self-explanatory, but may be tricky to understand in the beginning.
The C++ tools can be compiled with CMake, given that Boost, Eigen3, HDF5 and JSONCpp are installed.
For rendering [griegler/pyrender](https://github.com/griegler/pyrender) is required.
The easiest way to get started is to look into the configuration files outlining
the individual steps of the data pipeline and then look into the corresponding python scripts.

**Note that the data provided above only represents the final results of the
full pipeline.** You need to download the original datasets yourself.

Check the [top-level README](README.md) for license information. Additionally, see:

* [Semi-Convex-Hull Algorithm](http://www.cvlibs.net/software/semi_convex_hull/)
* [ICP](http://www.cvlibs.net/software/libicp/)
* [dimatura/binvox-rw-py](https://github.com/dimatura/binvox-rw-py)
* [alextsui05/blender-off-addon](https://github.com/alextsui05/blender-off-addon)
* [christopherbatty/SDFGen](https://github.com/christopherbatty/SDFGen)
* [Tomas_Akenine-Möller](http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/)
* [pyrender](https://github.com/griegler/pyrender)

## License

Note that the data is based on [ShapeNet](https://www.shapenet.org/terms) [1],
and [KITTI](http://www.cvlibs.net/datasets/kitti/) [2].
Check the corresponding websites for licenses.
The derived benchmarks are licensed as
[CC BY-NC-SA 3.0](Attribution-https://creativecommons.org/licenses/by-nc-sa/3.0/).

The code includes snippets from the following repositories:

* [pyrender](https://github.com/griegler/pyrender)
* [pyfusion](https://github.com/griegler/pyfusion)
* [Tomas_Akenine-Möller](http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/)
* [KITTI](http://www.cvlibs.net/datasets/kitti/)
* [Box-Ray Intersection](http://www.cs.utah.edu/~awilliam/box/)
* [Tomas Akenine-Möller Code](http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/)
* [griegler/pyrender](https://github.com/griegler/pyrender)
* [christopherbatty/SDFGen](https://github.com/christopherbatty/SDFGen)
* [High-Resolution Timer](http://www.songho.ca/misc/timer/timer.html)
* [Tronic/cmake-modules](https://github.com/Tronic/cmake-modules)
* [dimatura/binvox-rw-py](https://github.com/dimatura/binvox-rw-py)
* [alextsui05/blender-off-addon](https://github.com/alextsui05/blender-off-addon)

The remaining code is licensed as follows:

Copyright (c) 2018 David Stutz

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.