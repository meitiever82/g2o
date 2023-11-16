// g2o - General Graph Optimization
// Copyright (C) 2012 R. Kümmerle
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/core/eigen.hpp>

#include "g2o/core/auto_differentiation.h"
#include "g2o/core/base_unary_edge.h"
#include "g2o/core/base_vertex.h"
#include "g2o/core/optimization_algorithm_factory.h"
#include "g2o/core/sparse_optimizer.h"
#include "g2o/stuff/command_args.h"
#include "g2o/stuff/sampler.h"

#include "data.h"
#include "cyRepresentation.h"

using namespace std;

G2O_USE_OPTIMIZATION_LIBRARY(dense);

int main(int argc, char** argv) {
  int numPoints;
  int maxIterations;
  bool verbose;
  g2o::CommandArgs arg;
  arg.param("numPoints", numPoints, 100,
            "number of points sampled from the circle");
  arg.param("i", maxIterations, 10, "perform n iterations");
  arg.param("v", verbose, false, "verbose output of the optimization process");
  arg.parseArgs(argc, argv);

  std::vector<Eigen::Vector3d> landmarks;
  std::vector<Eigen::Matrix4d> v_Twc;
  std::vector<Eigen::Matrix4d> v_noisyTwc;

  createLandmarks(landmarks);
  std::vector<Eigen::Vector3d> noisyLandmarks = addLandmarksNoise(landmarks);
  createCameraPose(v_Twc, v_noisyTwc);
  const size_t pose_num = v_Twc.size();

  //cv::Mat cv_K = (cv::Mat_<double>(3, 3) << 480, 0, 320, 0, 480, 240, 0, 0, 1);  
  cv::Mat cv_K(3, 3, CV_64FC1);
  // Fill the matrix with values
  cv_K.at<double>(0, 0) = 480;
  cv_K.at<double>(0, 1) = 0;
  cv_K.at<double>(0, 2) = 320;
  cv_K.at<double>(1, 0) = 0;
  cv_K.at<double>(1, 1) = 480;
  cv_K.at<double>(1, 2) = 240;
  cv_K.at<double>(2, 0) = 0;
  cv_K.at<double>(2, 1) = 0;
  cv_K.at<double>(2, 2) = 1;
  //cv::Mat cv_K(480, 0, 320, 0, 480, 240, 0, 0, 1);
  
  Eigen::Tensor<double, 3, Eigen::RowMajor> K;
  cv::cv2eigen(cv_K, K);

  std::vector<Eigen::Vector2i> features_curr;
  // Setup optimizer
  g2o::SparseOptimizer optimizer;

  return 0;
}
