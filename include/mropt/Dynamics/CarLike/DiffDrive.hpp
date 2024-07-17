#ifndef DIFFDRIVE_H
#define DIFFDRIVE_H
#pragma once

#include "../ode.hpp"
#include "mropt/StateSpace/SE2/SE2.hpp"
#include "mropt/ControlSpace/R2/VW.hpp"
#include "mropt/RobotShape/CircleRobot.h"
namespace mropt::Dynamics::CarLike {
class DiffDrive : public ode {
private:
  // 表示差分驱动机器人的轮距
  double L_;
protected:
  DiffDrive() : ode() {}
public:

  //调用了ode的默认构造函数
  DiffDrive(const std::shared_ptr<mropt::ControlSpace::VW> &vw,
            const std::shared_ptr<mropt::StateSpace::SE2> &se2,
            const std::shared_ptr<mropt::RobotShape::CircleRobot> &cr)
      : ode(vw, se2, cr), L_(cr->getL()) {
    //Model
    // SX::vertcat是一个静态成员函数, 用于垂直连接多个符号表达式, 形成一个符号向量
    SX X_dot = SX::vertcat(
        {vw->v_ode * cos(se2->o_ode),
         vw->v_ode * sin(se2->o_ode),
         (1 / L_) * vw->w_ode});
    f_ = Function("f", {se2->X_ode(), vw->U_ode()}, {X_dot});
    jac_f_ = f_.jacobian();

  }
  std::shared_ptr<mropt::Dynamics::ode> clone(
      const std::shared_ptr<mropt::StateSpace::State> &state,
      const std::shared_ptr<mropt::ControlSpace::Control> &control,
      const std::shared_ptr<mropt::RobotShape::Footprint> &footprint) const override;
  ~DiffDrive() = default;
};
}
#endif