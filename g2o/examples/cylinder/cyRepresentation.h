#ifndef VO_SIMULATE_CYLINDER_H
#define VO_SIMULATE_CYLINDER_H

#include <cmath>
#include <chrono>
#include <Eigen/Core>
#include <opencv2/core/core.hpp>

#include <g2o/core/base_vertex.h>
#include <g2o/core/base_unary_edge.h>
#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/core/optimization_algorithm_gauss_newton.h>
#include <g2o/core/optimization_algorithm_dogleg.h>
#include <g2o/solvers/dense/linear_solver_dense.h>
#include <g2o/core/base_binary_edge.h>
#include <g2o/types/sba/types_sba.h>

#define pi 3.1415926535

class CylinderFittingVertex : public g2o::BaseVertex<5, Eigen::Matrix<double,5,1>> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    // 重置
    virtual void setToOriginImpl() override {
      _estimate << 0.0, 0.0, 0.0, 0.0, 10.0;
    }

    // 更新
    virtual void oplusImpl(const double *update) override {
      Eigen::Matrix<double,5,1> last=_estimate;
      _estimate += Eigen::Matrix<double,5,1>(update);
    }

    // 存盘和读盘：留空
    bool read(std::istream& /*is*/) override { return false; }
    bool write(std::ostream& /*os*/) const override { return false; }
};

class CylinderFittingEdge : public g2o::BaseBinaryEdge<1, double, g2o::VertexPointXYZ, CylinderFittingVertex> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    bool read(std::istream& /*is*/) override {
      std::cerr << __PRETTY_FUNCTION__ << " not implemented yet" << std::endl;
      return false;
    }

    bool write(std::ostream& /*os*/) const override {
      std::cerr << __PRETTY_FUNCTION__ << " not implemented yet" << std::endl;
      return false;
    }

    // 计算曲线模型误差
    void computeError() override {
      const g2o::VertexPointXYZ *vPoint = static_cast<const g2o::VertexPointXYZ *> (_vertices[0]);
      const CylinderFittingVertex *v = static_cast<const CylinderFittingVertex *> (_vertices[1]);
      const Eigen::Matrix<double,5,1> abc = v->estimate();
    }

    void linearizeOplus() override{ 
      const g2o::VertexPointXYZ* vPoint = static_cast<const g2o::VertexPointXYZ*>(_vertices[0]);
      const CylinderFittingVertex* v = static_cast<const CylinderFittingVertex*>(_vertices[1]);
      const Eigen::Matrix<double, 5, 1> abc = v->estimate();
    }

};

#endif //VO_SIMULATE_CYLINDER_H

