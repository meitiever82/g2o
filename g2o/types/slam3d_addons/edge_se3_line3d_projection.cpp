// g2o - General Graph Optimization
// Copyright (C) 2011 R. Kuemmerle, G. Grisetti, W. Burgard
// All rights reserved.
//
// Fully analytical Jacobians implementation for SE3-Line3D projection edge

#include "edge_se3_line3d_projection.h"

#include <cmath>
#include <limits>

#include "g2o/stuff/misc.h"

namespace g2o {

// ============================================================================
// Constructor
// ============================================================================
EdgeSE3Line3DProjection::EdgeSE3Line3DProjection()
    : BaseBinaryEdge<2, Line2D, VertexSE3, VertexLine3DTangent>() {}

// ============================================================================
// Static helper: skew-symmetric matrix
// ============================================================================
Matrix3 EdgeSE3Line3DProjection::skew(const Vector3& v) {
  Matrix3 S;
  S <<     0, -v.z(),  v.y(),
       v.z(),      0, -v.x(),
      -v.y(),  v.x(),      0;
  return S;
}

// ============================================================================
// Static helper: compute orthonormal basis perpendicular to d
// ============================================================================
void EdgeSE3Line3DProjection::computeOrthonormalBasis(
    const Vector3& d, Vector3& e1, Vector3& e2) {
  // Choose axis least aligned with d for numerical stability
  if (std::abs(d.x()) < std::abs(d.y()) && std::abs(d.x()) < std::abs(d.z())) {
    e1 = d.cross(Vector3::UnitX()).normalized();
  } else if (std::abs(d.y()) < std::abs(d.z())) {
    e1 = d.cross(Vector3::UnitY()).normalized();
  } else {
    e1 = d.cross(Vector3::UnitZ()).normalized();
  }
  e2 = d.cross(e1);  // Already unit since d and e1 are orthonormal
}

// ============================================================================
// Degenerate line detection (relaxed thresholds for robustness)
// ============================================================================
bool EdgeSE3Line3DProjection::isDegenerateLine(
    const Vector3& d, const Vector3& w) const {
  // Direction near zero
  if (d.squaredNorm() < 1e-10) return true;

  // Moment explosion
  if (w.squaredNorm() > 1e10) return true;

  // Note: We don't check orthogonality or unit norm here because
  // g2o's Line3D transformation may introduce small numerical errors.
  // The projection function handles these gracefully.

  return false;
}

// ============================================================================
// Degenerate projection detection
// ============================================================================
bool EdgeSE3Line3DProjection::isDegenerateProjection(
    const Vector3& l_pred) const {
  if (!std::isfinite(l_pred(0))) return true;
  double norm_sq = l_pred(0) * l_pred(0) + l_pred(1) * l_pred(1);
  return norm_sq < 1e-12;
}

// ============================================================================
// Project 3D line to normalized plane
// ============================================================================
Vector3 EdgeSE3Line3DProjection::projectToNormalizedPlane(
    const Line3D& line_C) const {
  Vector3 d_C = line_C.d();  // Direction (unit vector)
  Vector3 w_C = line_C.w();  // Moment

  // Project to normalized plane z=1: l = (n1, n2, rho)
  // Using Plücker line projection formula
  double n1_raw = d_C(2) * w_C(1) - d_C(1) * w_C(2);
  double n2_raw = d_C(0) * w_C(2) - d_C(2) * w_C(0);
  double rho_raw = d_C(1) * w_C(0) - d_C(0) * w_C(1);

  double norm = std::sqrt(n1_raw * n1_raw + n2_raw * n2_raw);

  // Degenerate: line parallel to optical axis
  if (norm < 1e-6) {
    return Vector3(std::numeric_limits<double>::quiet_NaN(), 0, 0);
  }

  return Vector3(n1_raw / norm, n2_raw / norm, rho_raw / norm);
}

// ============================================================================
// computeError
// ============================================================================
void EdgeSE3Line3DProjection::computeError() {
  const VertexSE3* v_pose = static_cast<const VertexSE3*>(_vertices[0]);
  const VertexLine3DTangent* v_line = static_cast<const VertexLine3DTangent*>(_vertices[1]);

  Isometry3 T_WtoC = v_pose->estimate().inverse();
  Line3D line_W = v_line->estimate();

  if (isDegenerateLine(line_W.d(), line_W.w())) {
    _error.setZero();
    return;
  }

  Line3D line_C = T_WtoC * line_W;

  if (isDegenerateLine(line_C.d(), line_C.w())) {
    _error.setZero();
    return;
  }

  Vector3 l_pred = projectToNormalizedPlane(line_C);

  if (isDegenerateProjection(l_pred)) {
    _error.setZero();
    return;
  }

  double n1_pred = l_pred(0);
  double n2_pred = l_pred(1);
  double rho_pred = l_pred(2);

  double theta_obs = _measurement(0);
  double rho_obs = _measurement(1);
  double n1_obs = std::cos(theta_obs);
  double n2_obs = std::sin(theta_obs);

  double dot = n1_pred * n1_obs + n2_pred * n2_obs;
  double cross = n1_pred * n2_obs - n2_pred * n1_obs;
  _error(0) = std::atan2(cross, dot);
  _error(1) = rho_pred - rho_obs;

  if (!std::isfinite(_error(0))) _error(0) = 0.0;
  if (!std::isfinite(_error(1))) _error(1) = 0.0;
}

#if 0
// ============================================================================
// Numerical Jacobian computation (for verification)
// ============================================================================
void EdgeSE3Line3DProjection::linearizeOplus() {
  // 使用数值微分计算 Jacobian
  const double delta = 1e-6;

  // 保存当前误差
  computeError();
  Vector2 error0 = _error;

  // ========================================================================
  // 数值计算位姿 Jacobian (2×6)
  // ========================================================================
  VertexSE3* v_pose = static_cast<VertexSE3*>(_vertices[0]);
  Isometry3 pose_backup = v_pose->estimate();

  for (int i = 0; i < 6; ++i) {
    Vector6 delta_vec = Vector6::Zero();
    delta_vec(i) = delta;

    // 正向扰动
    v_pose->oplus(delta_vec.data());
    computeError();
    Vector2 error_plus = _error;

    // 恢复
    v_pose->setEstimate(pose_backup);

    // 负向扰动
    delta_vec(i) = -delta;
    v_pose->oplus(delta_vec.data());
    computeError();
    Vector2 error_minus = _error;

    // 恢复
    v_pose->setEstimate(pose_backup);

    // 中心差分
    _jacobianOplusXi.col(i) = (error_plus - error_minus) / (2.0 * delta);
  }

  // ========================================================================
  // 数值计算线 Jacobian (2×4)
  // ========================================================================
  VertexLine3DTangent* v_line = static_cast<VertexLine3DTangent*>(_vertices[1]);
  Line3D line_backup = v_line->estimate();

  for (int i = 0; i < 4; ++i) {
    Eigen::Vector4d delta_vec = Eigen::Vector4d::Zero();
    delta_vec(i) = delta;

    // 正向扰动
    v_line->oplus(delta_vec.data());
    computeError();
    Vector2 error_plus = _error;

    // 恢复
    v_line->setEstimate(line_backup);

    // 负向扰动
    delta_vec(i) = -delta;
    v_line->oplus(delta_vec.data());
    computeError();
    Vector2 error_minus = _error;

    // 恢复
    v_line->setEstimate(line_backup);

    // 中心差分
    _jacobianOplusXj.col(i) = (error_plus - error_minus) / (2.0 * delta);
  }

  // 恢复原始误差
  computeError();
}
#endif

// ============================================================================
// Analytical Jacobian computation
// ============================================================================
void EdgeSE3Line3DProjection::linearizeOplus() {
  const VertexSE3* v_pose = static_cast<const VertexSE3*>(_vertices[0]);
  const VertexLine3DTangent* v_line = static_cast<const VertexLine3DTangent*>(_vertices[1]);

  Isometry3 T_CtoW = v_pose->estimate();
  Isometry3 T_WtoC = T_CtoW.inverse();
  Matrix3 R_WtoC = T_WtoC.rotation();
  Vector3 t_WtoC = T_WtoC.translation();

  Line3D line_W = v_line->estimate();
  Vector3 d_W = line_W.d();
  Vector3 w_W = line_W.w();

  // Degenerate case: set Jacobians to zero
  if (isDegenerateLine(d_W, w_W)) {
    _jacobianOplusXi.setZero();
    _jacobianOplusXj.setZero();
    return;
  }

  Line3D line_C = T_WtoC * line_W;
  Vector3 d_C = line_C.d();
  Vector3 w_C = line_C.w();

  if (isDegenerateLine(d_C, w_C)) {
    _jacobianOplusXi.setZero();
    _jacobianOplusXj.setZero();
    return;
  }

  // Project to normalized plane
  double n1_raw = d_C(2) * w_C(1) - d_C(1) * w_C(2);
  double n2_raw = d_C(0) * w_C(2) - d_C(2) * w_C(0);
  double rho_raw = d_C(1) * w_C(0) - d_C(0) * w_C(1);

  double norm = std::sqrt(n1_raw * n1_raw + n2_raw * n2_raw);
  if (norm < 1e-6) {
    _jacobianOplusXi.setZero();
    _jacobianOplusXj.setZero();
    return;
  }

  double inv_norm = 1.0 / norm;
  double n1 = n1_raw * inv_norm;
  double n2 = n2_raw * inv_norm;
  double rho = rho_raw * inv_norm;

  double theta_obs = _measurement(0);
  double n1_obs = std::cos(theta_obs);
  double n2_obs = std::sin(theta_obs);

  double dot = n1 * n1_obs + n2 * n2_obs;
  double cross = n1 * n2_obs - n2 * n1_obs;
  double denom = dot * dot + cross * cross;

  // ========================================================================
  // Layer 1: ∂error/∂l_pred (2×3)
  // ========================================================================
  Eigen::Matrix<double, 2, 3> de_dl;
  de_dl(0, 0) = (n2_obs * dot - n1_obs * cross) / denom;
  de_dl(0, 1) = (-n1_obs * dot - n2_obs * cross) / denom;
  de_dl(0, 2) = 0.0;
  de_dl(1, 0) = 0.0;
  de_dl(1, 1) = 0.0;
  de_dl(1, 2) = 1.0;

  // ========================================================================
  // Layer 2: ∂l_pred/∂l_raw (3×3) - normalization
  // ========================================================================
  Matrix3 dl_draw;
  dl_draw(0, 0) = (1.0 - n1 * n1) * inv_norm;
  dl_draw(0, 1) = (-n1 * n2) * inv_norm;
  dl_draw(0, 2) = 0.0;
  dl_draw(1, 0) = (-n1 * n2) * inv_norm;
  dl_draw(1, 1) = (1.0 - n2 * n2) * inv_norm;
  dl_draw(1, 2) = 0.0;
  dl_draw(2, 0) = (-rho * n1) * inv_norm;
  dl_draw(2, 1) = (-rho * n2) * inv_norm;
  dl_draw(2, 2) = inv_norm;

  // ========================================================================
  // Layer 3: ∂l_raw/∂L_C (3×6) - cross product
  // ========================================================================
  Eigen::Matrix<double, 3, 6> draw_dLC;
  draw_dLC(0, 0) = 0;        draw_dLC(0, 1) = -w_C(2);  draw_dLC(0, 2) = w_C(1);
  draw_dLC(1, 0) = w_C(2);   draw_dLC(1, 1) = 0;        draw_dLC(1, 2) = -w_C(0);
  draw_dLC(2, 0) = -w_C(1);  draw_dLC(2, 1) = w_C(0);   draw_dLC(2, 2) = 0;
  draw_dLC(0, 3) = 0;        draw_dLC(0, 4) = d_C(2);   draw_dLC(0, 5) = -d_C(1);
  draw_dLC(1, 3) = -d_C(2);  draw_dLC(1, 4) = 0;        draw_dLC(1, 5) = d_C(0);
  draw_dLC(2, 3) = d_C(1);   draw_dLC(2, 4) = -d_C(0);  draw_dLC(2, 5) = 0;

  Eigen::Matrix<double, 2, 6> de_dLC = de_dl * dl_draw * draw_dLC;

  // ========================================================================
  // Pose Jacobian: ∂error/∂ξ_pose (2×6)
  // g2o: T_CtoW' = T_CtoW * exp(δt, δq), δω = 2*δq
  // ========================================================================
  Matrix3 d_C_skew = skew(d_C);
  Matrix3 w_C_skew = skew(w_C);

  Eigen::Matrix<double, 6, 6> dLC_dpose;
  dLC_dpose.setZero();
  dLC_dpose.block<3, 3>(3, 0) = d_C_skew;       // ∂w_C/∂t = [d_C]×
  dLC_dpose.block<3, 3>(0, 3) = 2.0 * d_C_skew; // ∂d_C/∂q = 2*[d_C]×
  dLC_dpose.block<3, 3>(3, 3) = 2.0 * w_C_skew; // ∂w_C/∂q = 2*[w_C]×

  _jacobianOplusXi = de_dLC * dLC_dpose;

  // ========================================================================
  // Line Jacobian: ∂error/∂ξ_line (2×4)
  // ========================================================================
  Eigen::Matrix<double, 6, 6> dLC_dLW;
  dLC_dLW.setZero();
  dLC_dLW.block<3, 3>(0, 0) = R_WtoC;
  Matrix3 t_skew = skew(t_WtoC);
  dLC_dLW.block<3, 3>(3, 0) = t_skew * R_WtoC;
  dLC_dLW.block<3, 3>(3, 3) = R_WtoC;

  Eigen::Matrix<double, 2, 6> de_dLW = de_dLC * dLC_dLW;

  Vector3 e1_W, e2_W;
  computeOrthonormalBasis(d_W, e1_W, e2_W);

  Eigen::Matrix<double, 6, 4> dLW_dxi;
  dLW_dxi.setZero();
  dLW_dxi.block<3, 1>(0, 0) = -e2_W;           // ∂d/∂u1
  dLW_dxi.block<3, 1>(0, 1) =  e1_W;           // ∂d/∂u2
  dLW_dxi.block<3, 1>(3, 0) = e1_W.cross(w_W); // ∂w/∂u1
  dLW_dxi.block<3, 1>(3, 1) = e2_W.cross(w_W); // ∂w/∂u2
  dLW_dxi.block<3, 1>(3, 2) = e1_W;            // ∂w/∂v1
  dLW_dxi.block<3, 1>(3, 3) = e2_W;            // ∂w/∂v2

  _jacobianOplusXj = de_dLW * dLW_dxi;

  // NaN/Inf protection
  for (int i = 0; i < 6; ++i) {
    if (!std::isfinite(_jacobianOplusXi(0, i)) || !std::isfinite(_jacobianOplusXi(1, i))) {
      _jacobianOplusXi.col(i).setZero();
    }
  }
  for (int i = 0; i < 4; ++i) {
    if (!std::isfinite(_jacobianOplusXj(0, i)) || !std::isfinite(_jacobianOplusXj(1, i))) {
      _jacobianOplusXj.col(i).setZero();
    }
  }
}

// ============================================================================
// I/O
// ============================================================================
bool EdgeSE3Line3DProjection::read(std::istream& is) {
  double theta, rho;
  is >> theta >> rho;
  _measurement = Line2D(theta, rho);
  for (int i = 0; i < 2; ++i)
    for (int j = i; j < 2; ++j) {
      is >> information()(i, j);
      if (i != j) information()(j, i) = information()(i, j);
    }
  return true;
}

bool EdgeSE3Line3DProjection::write(std::ostream& os) const {
  os << _measurement(0) << " " << _measurement(1);
  for (int i = 0; i < 2; ++i)
    for (int j = i; j < 2; ++j)
      os << " " << information()(i, j);
  return os.good();
}

}  // namespace g2o