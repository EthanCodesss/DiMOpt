#ifndef STATE_H
#define STATE_H
#pragma once

#include <casadi/casadi.hpp>
#include <list>
namespace mropt::Problem{ class Robot;}

namespace mropt::StateSpace {
class State {
protected:
  // 友元类
  friend class mropt::Problem::Robot;
  // 状态变量的数量
  int nx_{0};
  // 时间步数
  int Nx_{0};
  // 存储状态变量的CasADi矩阵
  casadi::MX X_;
  // 状态变量的下边界和上边界
  std::vector<double> lb_, ub_;

  std::list<casadi::MX> get_constraints();
public:
  // 枚举类型, 表示状态变量的位置
  enum class POS { x = 0, y = 1, o = 2 };
  casadi::Slice all;

  // casadi::Slice 是CasADi库中用于表示矩阵或者向量中连续子区间的类
  casadi::Slice xy_slice;
  // 内部结构体, 表示状态的三个分量
  struct state {
    double x, y, o;
  };
  // 在state初始化时, 初始化了一个切片, 这个切片用于后续访问和操作状态变量中的特定分量
  State() { xy_slice = casadi::Slice((int) POS::x, (int) POS::y + 1); }

  // 接收一个右值引用参数, 使用std::move 实现资源的转移
  void setX(casadi::MX &&X) {
    int _nx = X.size1();
    Nx_ = X.size2();
    if (_nx != nx_) {
      std::cerr << "[State] State space requires dim " << nx_ << std::endl;
    }
    X_ = std::move(X);
  }

  casadi::MX X() { return X_; };
  virtual int nx() = 0;
  int Nx() { return Nx_; };
  virtual casadi::SX X_ode() = 0;
  virtual casadi::SX XY_ode() = 0;
  //virtual casadi::MX gX() = 0;
  virtual ~State() = default;

  virtual casadi::SX get_weights() const = 0;
  virtual casadi::SX get_std_values() const = 0;
  virtual State &set_weights_std_values(std::vector<double> weight, std::vector<double> vars_std) = 0;

  virtual casadi::MX xy() = 0;
  virtual casadi::MX x() = 0;
  virtual casadi::MX y() = 0;
  virtual casadi::MX o() = 0;

  const casadi::MX get_X_0() { return X_(all, 0); }
  const casadi::MX get_X_f() { return X_(all, Nx_ - 1); }

  void set_bounds(casadi::Opti &ocp);

  void set_lower_bounds(std::vector<double> lb);
  void set_upper_bounds(std::vector<double> ub);

  virtual std::shared_ptr<State> clone() const = 0;
};
}
#endif