// Implementation of the GCPR 2016 Paper "Joint Object Pose Estimation and
// Shape Reconstruction in Urban Street Scenes Using 3D Shape Priors" by Engelmann et al.
// Copyright (C) 2016  Francis Engelmann - Visual Computing Institute RWTH Aachen University
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// C/C++ includes
#include <iostream>
#include <string>
#include <iostream>
#include <fstream>
#include <time.h>
#include <cmath>
#include <memory>
#include <unordered_map>
#include <tuple>
#include <algorithm>
#include <ctime>

// VIZ includes
#include <viz/viz.h>

// Own includes
#include "io.h"
#include "utils.h"
#include "tracking.h"
#include "optimization/shape.h"
#include "optimization/pose.h"

// KITTI includes
#include "kitti.h"

struct parameters {
  /** \brief File of point cloud in TXT format holding the point cloud extracted from a bounding box
   * and scaled to [0,widht] x [0,height] x [0,depth].
   */
  std::string input_directory;
  /** \brief Output file to save the computed TSDF as vector.
   */
  std::string output_directory;
  /** \brief Whether to visualize results.
   */
  bool interactive = false;
  /** \brief Only process this index.
   */
  int index = -1;

  int height = 24; // y axis upwards
  int width = 54; // x axis right
  int depth = 24; // z axis forward

  parameters(int argc, char** argv) {
    if (argc < 2+1) { std::cerr << "Not enough parameters!" << std::endl; }
    input_directory = std::string{argv[1]} + "/";
    output_directory = std::string{argv[2]} + "/";
    if (!boost::filesystem::is_directory(output_directory)) { boost::filesystem::create_directories(output_directory); }
    if (argc > 2+1) { interactive = std::string{argv[3]} == "1"; }
    if (argc > 3+1) { index = atoi(argv[4]); }
  }

  void print(void) {
    std::cout << "input_directory: " << input_directory << std::endl;
    std::cout << "output_directory: " << output_directory << std::endl;
  }
};

int main(int argc, char** argv) {

  // Read parameters
  parameters params(argc, argv);
  bool print_params = true;
  if (print_params) params.print();

  // Read SVD decompositions
  Eigen::MatrixXd& S = gvl::ShapeSpacePCA::instance().S;
  Eigen::MatrixXd& V = gvl::ShapeSpacePCA::instance().V;
  Eigen::MatrixXd& mean = gvl::ShapeSpacePCA::instance().mean;
  gvl::readMatrixFromBinaryFile(path_S, S);
  gvl::readMatrixFromBinaryFile(path_V, V);
  gvl::readMatrixFromBinaryFile(path_mean_shape, mean);

  // Read point cloud files
  std::map<int, std::string> txt_files;
  kitti::read_directory(params.input_directory, txt_files, ".txt");
  std::cout << "Read "  << txt_files.size() << " point clouds" << std::endl;

  std::vector<int> indices;
  if (params.index >= 0) {
    indices.push_back(params.index);
  }
  else {
    for (unsigned int i = 0; i < txt_files.size(); i++) {
      indices.push_back(i);
    }
  }

  double total_time = 0;
  for (unsigned int i = 0; i < indices.size(); i++) {
    int index = indices[i];

    std::cout << "Processing " << txt_files[index] << std::endl;
    auto pointcloud = gvl::read_pointcloudTxt(txt_files[index]);
    // center at (0,0,0)
    Eigen::Matrix4d trafo = Eigen::Matrix4d::Identity();
    trafo(0,3) -= params.width/2;
    trafo(1,3) -= params.height/2;
    trafo(2,3) -= params.depth/2;
    pointcloud->transform(trafo);
    // scale to meters and translate
    trafo = Eigen::Matrix4d::Identity()*6/params.width;
    trafo(1,3) = 1;
    pointcloud->transform(trafo);
    // axes swap
    trafo = Eigen::Matrix4d::Zero();
    trafo(0,2) = 1;
    trafo(1,1) = -1; // y points downward in KITTI camera coordinates
    trafo(2,0) = -1;
    pointcloud->transform(trafo);

    //for (unsigned int i = 0; i < pointcloud->points.size(); i++) {
    //  std::cout << pointcloud->points.at(i).x << " " << pointcloud->points.at(i).y << " " << pointcloud->points.at(i).z << std::endl;
    //}

    auto det = std::make_shared<gvl::Detection>();
    det->boundingbox = std::make_shared<gvl::BoundingBox>();
    det->frame_id = 0;
    det->translation = Eigen::Vector3d::Zero();
    det->rotation_y = 0;
    det->pose = Eigen::Matrix4d::Identity();
    det->z = Eigen::VectorXd::Zero(r,1); // r = 5, subspace dim
    det->shape = mean;
    det->pointcloud = pointcloud;

    // Optimization - iterate alternativly until convergence
    std::clock_t start = std::clock();
    double total_cost_curr=DBL_MAX/2;
    double total_cost_prev=DBL_MAX;
    double delta = 0.1;
    int t = 0;
    while (total_cost_prev - total_cost_curr > delta && t<1000) {
      total_cost_prev = total_cost_curr;
      double pose_cost = gvl::optimize_pose(*det, false, &t);
      double shape_cost = gvl::optimize_shape(*det, false, &t);
      total_cost_curr = pose_cost + shape_cost;
      std::cout << "Iteration: " << t++ << " Cost: " << total_cost_curr << std::endl;
    }

    std::clock_t end = std::clock();
    double elapsed = double(end - start) / CLOCKS_PER_SEC;
    std::cout << "Optimization took " << elapsed << " seconds" << std::endl;
    total_time += elapsed;

    std::string result_file = params.output_directory + std::to_string(index) + ".txt";
    gvl::writeMatrixToFile(det->shape, result_file);

    std::string pose_file = params.output_directory + std::to_string(index) + "_pose.txt";
    gvl::writeMatrixToFile(det->pose, pose_file);

    if (params.interactive) { // show 3d of result
      viz::Visualization v(1600,800, 2);
      v.connectCameras(0,1);
      v.lookFromAt(Eigen::Vector3d(0,-7,-20),Eigen::Vector3d(0,0,0));

      Eigen::Vector3f text_color(255,255,255);

      // Left viewer
      v.setActiveRenderer(0);
      v.setBackgroundColor(0,0,0);
      v.addText2d("Input", 30, 30, text_color);
      v.addPointcloud(pointcloud->getVertices(),Eigen::Vector3f(255, 255, 255), 3);

      // Right viewer
      v.setActiveRenderer(1);
      v.setBackgroundColor(0,0,0);
      v.addText2d("Our result", 30, 30, text_color);
      v.addPointcloud(pointcloud->getVertices(),Eigen::Vector3f(255, 255, 255), 3);

      v.addZeroSurfaceFromGrid(det->shape, det->pose, Eigen::Vector3f(255,255,255), false, 0,0.1,1);
      v.addZeroSurfaceFromGrid(det->shape, det->pose, Eigen::Vector3f(0,0,0), true, 0,0.1,2);
      v.show();
    }
  }

  std::cout << "Total optimization time: " << total_time << std::endl;

  return 0;
}
