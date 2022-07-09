# Weakly Supervised Shape Completion

This repository contains the code for the weakly-supervised shape completion
method, called amortized maximum likelihood (AML), described in:

    @inproceedings{Stutz2018CVPR,
        title = {Learning 3D Shape Completion from Laser Scan Data with Weak Supervision },
        author = {Stutz, David and Geiger, Andreas},
        booktitle = {IEEE Conference on Computer Vision and Pattern Recognition (CVPR)},
        publisher = {IEEE Computer Society},
        year = {2018}
    }
    @misc{Stutz2017,
        author = {David Stutz},
        title = {Learning Shape Completion from Bounding Boxes with CAD Shape Priors},
        month = {September},
        year = {2017},
        institution = {RWTH Aachen University},
        address = {Aachen, Germany},
        howpublished = {http://davidstutz.de/},
    }

If you use this code for your research, please cite the above papers.

See [davidstutz/cvpr2018-shape-completion](https://github.com/davidstutz/cvpr2018-shape-completion)
for the LaTeX source of the paper and additional repositories.

![Illustration of the proposed approach.](screenshot.png?raw=true "Illustration of the proposed approach.")

## Overview

The paper proposes a weakly-supervised approach to shape completion on real
data, specifically [KITTI](http://www.cvlibs.net/datasets/kitti/).
In particular, a variational auto-encoder is
trained to learn a shape prior - on a set of cars extracted from
[ShapeNet](https://www.shapenet.org/). The generative model, this means
the decoder, is then fixed and a new encoder is trained to embed
observations - from point clouds in KITTI or synthetsized using depth images
on ShapeNet - in the same latent space. The encoder can be trained in an
unsupervised fashion - as we know the object category and bounding boxes on KITTI
the approach can be described as weakly-supervised. In particular, the encoder
predicts codes which match the prior on the latent space, a unit Gaussian,
and the loss between the generated shape and the observations makes sure
that the shape fits the observations. The overall approach is called
amortized maximum likelihood loss, as the encoder is trained
to minimize a maximum likelihood loss. As shape representation, occupancy
grids and signed distance functions are used.

In this repository we provide our implementation of the
amortized maximum likelihood approach, the maximum likelihood baseline,
the supervised baseline and related work [1].
The repository also contains an adapted version of

* [VisualComputingInstitute/ShapePriors_GCPR16](https://github.com/VisualComputingInstitute/ShapePriors_GCPR16)

## Installation

**[sumanthrao1997/daml_shape_completion_docker](https://github.com/sumanthrao1997/daml_shape_completion_docker) includes a docker to run the code in this repository.**

LUA/Torch requirements:

* Torch ([torch/distro](https://github.com/torch/distro) recommended);
* [deepmind/torch-hdf5](https://github.com/deepmind/torch-hdf5);
* [harningt/luajson](https://github.com/harningt/luajson);
* [luafilesystem](http://keplerproject.github.io/luafilesystem);
* [clementfarabet/lua---nnx](https://github.com/clementfarabet/lua---nnx) (already included in [torch/distro](https://github.com/torch/distro));
* [nicholas-leonard/cunnx](https://github.com/nicholas-leonard/cunnx);
* [davidstutz/torch-volumetric-nnup](https://github.com/davidstutz/torch-volumetric-nnup) (see the installation instructions in the repository).

Installing [deepmind/torch-hdf5](https://github.com/deepmind/torch-hdf5)
might be tricky. After building orch-hdf5,

    git clone https://github.com/deepmind/torch-hdf5
    cd torch-hdf5
    luarocks make hdf5-0-0.rockspec
    cd ..

it might be necessary to adapt the configuration in case you installed
HDF5 locally. For example, when installing hdf5 locally in `DB_PATH`,
`torch/install/share/lua/5.1/hdf5/config.lua` might need
to be adapted as follows:

    require('os')
  
    db_path = os.getenv("DB_PATH")
    hdf5._config = {
        HDF5_INCLUDE_PATH = db_path .. "/hdf5/hdf5/include/",
        HDF5_LIBRARIES = db_path .. "/hdf5/hdf5/lib/libhdf5_cpp.so;" .. db_path .. "/hdf5/hdf5/lib/libhdf5.so;/usr/lib/x86_64-linux-gnu/libpthread.so;/usr/lib/x86_64-linux-gnu/libz.so;/usr/lib/x86_64-linux-gnu/libdl.so;/usr/lib/x86_64-linux-gnu/libm.so"
    }

Make sure that `nnx`, `cunnx` and the volumetric nearest neighbor
upsampling layer works by following the instructions in
[davidstutz/torch-volumetric-nnup](https://github.com/davidstutz/torch-volumetric-nnup).

The remaining packages can easily be installed using luarocks. You can run

    th check_requirements.lua

to check the packages listed above.

Pyton requirements:

* NumPy;
* h5py;
* [PyMCubes](https://github.com/davidstutz/PyMCubes) (**make sure to use the `voxel_center` branch**).

For installing PyMCubes, follow the instructions [here](https://github.com/davidstutz/PyMCubes);
NumPy and h5py can be installed using `pip` and might themselves
have dependencies.

We also include an implementation of the method by Engelmann et al. [1].

    [1] Francis Engelmann, Jörg Stückler, Bastian Leibe:
        Joint Object Pose Estimation and Shape Reconstruction in Urban Street Scenes Using 3D Shape Priors. GCPR 2016: 219-230

First, make sure to install:

* [OpenCV 2.4.13.x](https://opencv.org/);
* [VTK 7.1.x](https://www.vtk.org/download/);
* [Ceres](http://ceres-solver.org/installation.html):
  * Eigen3;
  * [GLog](https://github.com/google/glog);
  * [GFlags](https://github.com/gflags/gflags);
  * ([Suitesparse](http://faculty.cse.tamu.edu/davis/suitesparse.html));

The installation of Ceres might be a bit annoying; so we provide our
installation scripts for some of the dependencies in `rw/dependencies` for further details.
These still need to be adapted (e.g. they assume that dependencies as well
as Ceres are installed in `$WORK/dev-box` where `$WORK` is a defined base directory).
Note that SuiteSparse is optional but significantly reduces runtime
(by roughly factor 2-4).

When all dependencies are installed, make sure to adapt the corresponding
CMake files in `rw/cmake_modules`. This means removing `NO_CMAKE_SYSTEM_PATH`
if necessary and inserting the correct paths to the installations.

Then:

    cd rw/external/viz
    mkdir build
    cd build
    cmake ..
    make
    # make sure VIZ is built correctly
    cd ../../
    mkdir build
    cd build
    cmake ..
    make
    # make sure KittiShapePrior and ShapeNetShapePrior are built correctly

Also make sure to download the pre-trained PCA shape prior
from [VisualComputingInstitute/ShapePriors_GCPR16](https://github.com/VisualComputingInstitute/ShapePriors_GCPR16).

For the C++ implementation of the evaluation tool (mesh-to-mesh and point-to-mesh distances),
follow the instructions here: [davidstutz/mesh-evaluation](https://github.com/davidstutz/mesh-evaluation);
essentially, the tool requires:

* CMake;
* Boost;
* Eigen;
* OpenMP;
* C++11.

Make sure to adapt the corresponding CMake modules, then run:

    cd mesh-evaluation
    mkdir build
    cd build
    cmake ..
    make

## Data

The data is derived from [ShapeNet](https://www.shapenet.org/terms)
and [KITTI](http://www.cvlibs.net/datasets/kitti/). For ShapeNet,
two datasets, in the paper referred to as SN-clean and SN-noisy,
were created. We also provide the simplified cars from ShapeNet, simplified
using [this semi-convex hull algorithm](http://www.cvlibs.net/software/semi_convex_hull/).

Download links:

* KITTI (2.7GB): [Amazon AWS](https://s3.eu-central-1.amazonaws.com/avg-shape-completion/cvpr2018_shape_completion_kitti.zip) [MPI-INF](https://datasets.d2.mpi-inf.mpg.de/cvpr2018-shape-completion/cvpr2018_shape_completion_kitti.zip)
* SN-clean (5.7GB): [Amazon AWS](https://s3.eu-central-1.amazonaws.com/avg-shape-completion/cvpr2018_shape_completion_clean.zip) [MPI-INF](https://datasets.d2.mpi-inf.mpg.de/cvpr2018-shape-completion/cvpr2018_shape_completion_clean.zip)
* SN-noisy (4.0GB): [Amazon AWS](https://s3.eu-central-1.amazonaws.com/avg-shape-completion/cvpr2018_shape_completion_noisy.zip) [MPI-INF](https://datasets.d2.mpi-inf.mpg.de/cvpr2018-shape-completion/cvpr2018_shape_completion_noisy.zip)
* Simplified ShapeNet Cars (38.5MB): [Amazon AWS](https://s3.eu-central-1.amazonaws.com/avg-shape-completion/cvpr2018_shape_completion_shapenet_models.zip) [MPI-INF](https://datasets.d2.mpi-inf.mpg.de/cvpr2018-shape-completion/cvpr2018_shape_completion_shapenet_models.zip)

**See [Data](data/README.md) for more details.**

Make sure to cite [3] and [4] in addition to this paper when using the data.

    [3] Andreas Geiger, Philip Lenz, Raquel Urtasun:
        Are we ready for autonomous driving? The KITTI vision benchmark suite. CVPR 2012: 3354-3361
    [4] Angel X. Chang, Thomas A. Funkhouser, Leonidas J. Guibas, Pat Hanrahan, Qi-Xing Huang, Zimo Li, Silvio Savarese, Manolis Savva, Shuran Song, Hao Su, Jianxiong Xiao, Li Yi, Fisher Yu:
        ShapeNet: An Information-Rich 3D Model Repository. CoRR abs/1512.03012 (2015)

## Models

The models can be downloaded from: [Models](https://s3.eu-central-1.amazonaws.com/avg-shape-completion/cvpr2018_shape_completion_models.zip).
The model for [1] is available at [VisualComputingInstitute/ShapePriors_GCPR16](https://github.com/VisualComputingInstitute/ShapePriors_GCPR16).

The downloaded ZIP-archive contains the pre-trained models as
`.dat` files for the shape prior (`vae`), the proposed weakly-supervised
approach (`daml`) and the supervised baseline (`sup`):

    clean/
    |- vae/
       |- prior_model.dat
    |- daml/
       |- inference_model.dat
    |- sup/
       |- inference_model.dat
    noisy/
    |- daml/
       |- inference_model.dat
    |- sup/
       |- inference_model.dat
    kitti/
    |- daml/
       |- inference_model.dat

These models have been saved using `data/tools/lua/compress_model_dat.lua`
in order to reduce their size.

For running a model, it is sufficient to extract the corresponding
model in the correct `base_directory` as specified in the configuration files.
For example, for running the shape prior, create the folder `vae/clean` (as
determined by `vae/clean.json`) and put `prior_model.dat` inside this folder.
Then:

    th vae_run.lua clean.json

For more details, see the training instructions below.

## Experiments

Make sure that data is downloaded and all requirements are met, for example run
`check_requirements.lua` and `check_requirements.py`; also
build the evaluation tool from
[davidstutz/mesh-evaluation](https://github.com/davidstutz/mesh-evaluation)
as described above.

### Shape Prior

For training and evaluating the proposed amortized maximum likelihood approach:

* Adapt the configuration file `vae/clean.json`; specifically, set the
  `data_directory` to the correct location of the downloaded "clean" dataset.
* Change to the `vae` directory and run

    ```
    th vae_train.lua clean.json
    ```

This will train a shape prior, to test whether training works, configuration options
such as `epochs` can be adapted to only run a few iterations. More details can
be found in `vae/clean.json` where all configuration options are commented.

### Shape Inference

The next step is to train a new encoder:

* Change to the `daml` directory.
* Create a directory `clean` and copy `prior_model.dat` from `vae/clean`
  (after training the shape prior!).
* Adapt the configuration file `daml/clean.json`; specifically,
  set the `data_directory` to the correct location of the downloaded
  "clean" dataset.
* Run

    ```
    th encoder_train.lua clean.json
    ```

Afterwards, a prior model and an inference model (in the form of the
corresponding `.dat` files) are available.

### Evaluation

For evaluating the shape prior:

* Change to the `vae` directory.
* Check that the `clean` directory contains `prior_model.dat` (indicating that
  training succeeded) and a set of `*_predictions.h5` files.
* Run

    ```
    python vae_test.py clean.json
    ```

* The tool will create a directory `clean/off` corresponding to the predicted
  meshes and write the evaluation, specifically the Hamming distance
  of the predicted occupancy, in `clean/results.txt`.
* Use [davidstutz/mesh-evaluation](https://github.com/davidstutz/mesh-evaluation) to evaluate the meshes in `clean/off`
  against the ground truth meshes as downloaded with the data.

For evaluating the inference model:

* Change to the `daml` directory.
* Check that the `clean` directory contains `prior_model.dat` and
  `inference_model.dat` (indicating that training succeeded) and a set of
  `*_predictions.h5` files.
* Run

    ```
    python encoder_test.py clean.json
    ```

* The tool will create a `clean/off` directory and evaluate occupancy
  in `clean/results.txt`.
* Use [davidstutz/mesh-evaluation](https://github.com/davidstutz/mesh-evaluation) to evaluate the meshes in `clean/off`
  against the ground truth meshes as downloaded with the data.

_Note that the occupancy predictions (i.e. the volumes) are in
`H x W x D = 24 x 54 x 24` (corresponding to height, width and depth) while
the predicted meshes are `W x H x D` coordinates. In contrast to many
other publications, we use height as the second dimension and depth as the
third, in detail, the axes are x = right, y = up, z = forward._

#### Maximum Likelihood Baseline

After training a shape prior, shape completion can be performed
using standard maximum likelihood - this was used as baseline in the paper.
Instructions:

* Change the directory to `ml`.
* Adapt `clean.json` and set `data_directory` to the correct
  location of the downloaded data.*
* Run

    ```
    th ml_train.lua clean.json
    ```
    
* For evaluation, run

    ```
    python ml_test.py clean.json
    ```

* Run [davidstutz/mesh-evaluation](https://github.com/davidstutz/mesh-evaluation) on the created `clean/off`
  file.

#### Supervised Baseline

To train the supervised baseline:

* Change the directory to `sup` and adapt `clean.json` as
  described above.
* Run

    ```
    th vae_train.lua clean.json
    ```

* For evaluation, run

    ```
    python vae_test.py clean.json
    ```

* Finally, use [davidstutz/mesh-evaluation](https://github.com/davidstutz/mesh-evaluation) to evaluate the meshes in `clean/off`
  against the ground truth meshes as downloaded with the data.

#### Engelmann et al. Baseline

First, make sure that the work by Engelmann et al. [1] can be compiled
as outlined in [Installation](#installation).

Then, two command line tools are provided:

* `KittiShapePrior` for running the approach on KITTI; arguments are
  the input directory with point clouds as `.txt` files, a `.txt` file containing
  the correpsonding bounding boxes, and the output directory.
* `ShapeNetShapePrior` for running the approach on ShapeNet ("clean" and "noisy");
  arguments are the input directory containing the point clouds as `.txt` files,
  and the output directory.

Running the approach on ShapeNet might look as follows:

    ./ShapeNetShapePrior /path/to/training_inference_txt_gt_10_48x64_24x54x24_clean_large output_directory

Subsequently, `rw/tools/shapenet_marching_cubes.py` can be used to obtain
meshes from the predicted signed distance functions:

    python shapenet_marching_cubes.py output_directory off_directory

For KITTI, the approach is similar; however, in addition to the input points,
the bounding boxes are required:

    ./KittiShapePrior /path/to/bounding_boxes_txt_validation_gt_padding_corrected_1_24x54x24/ /path/to/bounding_boxes_validation_gt_padding_corrected_1_24x54x24.txt output_directory

Similarly, `rw/tools/kitti_marching_cubes.py` expects the bounding boxes
as second argument as well.

### Visualization

The original meshes included in the data downloads
as well as the predicted meshes (of both the shape prior and shape inference)
can be visualized using [MeshLab](http://www.meshlab.net/).

The [davidstutz/bpy-visualization-utils](https://github.com/davidstutz/davidstutz/bpy-visualization-utils)
repository also provides utilities for visualization using Blender and Python.

## License

Licenses for source code and data corresponding to:

D. Stutz, A. Geiger. **Learning 3D Shape Completion from Laser Scan Data with Weak Supervision.** IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 2018.

Note that the source code and/or data is based on the following projects for which separate licenses apply:

* [ShapeNet](https://www.shapenet.org/terms)
* [KITTI](http://www.cvlibs.net/datasets/kitti/)
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

### Source Code

Copyright (c) 2018 David Stutz, Max-Planck-Gesellschaft

**Please read carefully the following terms and conditions and any accompanying documentation before you download and/or use this software and associated documentation files (the "Software").**

The authors hereby grant you a non-exclusive, non-transferable, free of charge right to copy, modify, merge, publish, distribute, and sublicense the Software for the sole purpose of performing non-commercial scientific research, non-commercial education, or non-commercial artistic projects.

Any other use, in particular any use for commercial purposes, is prohibited. This includes, without limitation, incorporation in a commercial product, use in a commercial service, or production of other artefacts for commercial purposes.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

You understand and agree that the authors are under no obligation to provide either maintenance services, update services, notices of latent defects, or corrections of defects with regard to the Software. The authors nevertheless reserve the right to update, modify, or discontinue the Software at any time.

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. You agree to cite the corresponding papers (see above) in documents and papers that report on research using the Software.

### Data

Copyright (c) 2018 David Stutz, Max-Planck-Gesellschaft

**Please read carefully the following terms and conditions and any accompanying documentation before you download and/or use the data (the "Dataset").**

The authors grant you a non-exclusive, non-transferable, free of charge right: To download the Dataset and use it on computers owned, leased or otherwise controlled by you and/or your organisation; To use the Dataset for the sole purpose of performing non-commercial scientific research, non-commercial education, or non-commercial artistic projects.

Any other use, in particular any use for commercial purposes, is prohibited. This includes, without limitation, incorporation in a commercial product, use in a commercial service, or production of other artefacts for commercial purposes.

Without prior written approval from the authors, the Dataset, in whole or in part, shall not be further distributed, published, copied, or disseminated in any way or form whatsoever, whether for profit or not. This includes further distributing, copying or disseminating to a different facility or organizational unit in the requesting university, organization, or company.

THE DATASET IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE DATASET OR THE USE OR OTHER DEALINGS IN THE DATASET.

You understand and agree that the authors are under no obligation to provide either maintenance services, update services, notices of latent defects, or corrections of defects with regard to the Dataset. The authors nevertheless reserve the right to update, modify, or discontinue the Dataset at any time.

You agree to cite the corresponding papers (see above) in documents and papers that report on research using the Dataset.
