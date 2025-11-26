// g2o - General Graph Optimization
// Copyright (C) 2011 R. Kuemmerle, G. Grisetti, W. Burgard
// All rights reserved.
//
// Edge for projecting 3D line to 2D normalized plane observation

#ifndef G2O_EDGE_SE3_LINE3D_PROJECTION_H
#define G2O_EDGE_SE3_LINE3D_PROJECTION_H

#include "g2o/config.h"
#include "g2o/core/base_binary_edge.h"
#include "g2o/types/slam3d/vertex_se3.h"
#include "g2o_types_slam3d_addons_api.h"
#include "line3d.h"
#include "vertex_line3d_tangent.h"

namespace g2o {

/**
 * @brief 2D line observation on normalized plane (theta, rho)
 */
class Line2D : public Vector2 {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Line2D() { setZero(); }

  Line2D(double theta, double rho) {
    (*this)(0) = theta;
    (*this)(1) = rho;
  }

  double theta() const { return (*this)(0); }
  double rho() const { return (*this)(1); }
};

/**
 * @brief Edge connecting SE3 pose and 3D line, with 2D projection observation
 *
 * Projects 3D line to normalized image plane using Plücker coordinates.
 * Observation: Line2D (theta, rho) - 2D line in polar form
 * Vertices: [0] VertexSE3 (camera pose), [1] VertexLine3D (3D line)
 */
class G2O_TYPES_SLAM3D_ADDONS_API EdgeSE3Line3DProjection
    : public BaseBinaryEdge<2, Line2D, VertexSE3, VertexLine3DTangent> {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  EdgeSE3Line3DProjection();

  /**
   * @brief Compute reprojection error
   *
   * 1. Transform world line to camera frame
   * 2. Project to normalized plane using Plücker cross product
   * 3. Compute angle and distance error
   */
  virtual void computeError();

  /**
   * @brief Compute Jacobian matrices
   *
   * Jacobian chain: error -> l_pred -> L_C -> L_W -> parameters
   * - J[0]: 2x6 for SE3 pose
   * - J[1]: 2x4 for Line3D parameters
   */
  virtual void linearizeOplus();

  virtual bool read(std::istream& is);
  virtual bool write(std::ostream& os) const;

  virtual bool setMeasurementData(const double* d) {
    _measurement(0) = d[0];
    _measurement(1) = d[1];
    return true;
  }

  virtual bool getMeasurementData(double* d) const {
    d[0] = _measurement(0);
    d[1] = _measurement(1);
    return true;
  }

  virtual int measurementDimension() const { return 2; }

 private:
  /**
   * @brief Project 3D line to 2D normalized plane
   * @param line_C Line in camera frame
   * @return [n1, n2, d] - normalized plane line parameters
   */
  Vector3 projectToNormalizedPlane(const Line3D& line_C) const;

  /**
   * @brief Compute skew-symmetric matrix from vector
   */
  static Matrix3 skew(const Vector3& v);

  /**
   * @brief Compute orthonormal basis perpendicular to direction
   */
  static void computeOrthonormalBasis(const Vector3& d, Vector3& e1, Vector3& e2);

  /**
   * @brief Check if line parameters are degenerate
   */
  bool isDegenerateLine(const Vector3& d, const Vector3& w) const;

  /**
   * @brief Check if projection is degenerate
   */
  bool isDegenerateProjection(const Vector3& l_pred) const;
};

}  // namespace g2o

#endif
