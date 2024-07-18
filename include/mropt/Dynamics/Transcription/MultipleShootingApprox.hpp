#ifndef MULTIPLESHOOTINGAPPROX_H
#define MULTIPLESHOOTINGAPPROX_H
#pragma once

#include <casadi/casadi.hpp>
#include "../Transcription.hpp"
#include "../OdeApprox.hpp"

namespace mropt::Dynamics {
class MultipleShootingApprox : public Transcription {

private:
public:
  explicit MultipleShootingApprox(const std::shared_ptr<OdeApprox> &ode_approx)
      : Transcription(ode_approx) {}
  ~MultipleShootingApprox() = default;

  std::vector<MX> get_constraints() override {
    // 建立一个vector
    std::vector<MX> constraints{};
    MX g_sum{0.0};
    for (int k = 0; k < N; ++k) {
      // 计算下一个状态
      const auto &x_next = ode_->state_space_->X()(all, k) +
          integrator(
              *(ode_approx_->fapprox(k)),
              (1 / (double) N) * tf,
              ode_->state_space_->X()(all, k),
              ode_->control_space_->U()(all, k));
      // 计算累积误差
      g_sum = g_sum + sum1(MX::abs(ode_->state_space_->X()(all, k + 1) - x_next));
      constraints.push_back(ode_->state_space_->X()(all, k + 1) - x_next);
    }
    J_model_ = Function("J_model", {ode_->state_space_->X(), ode_->control_space_->U(), tf}, {g_sum});
    return constraints;
  }

  void set_J_real() override {
    MX g_sum{0.0};
    MX g_max{0.0};
    for (int k = 0; k < N; ++k) {
      const auto &x_next = ode_->state_space_->X()(all, k) + integrator(
          ode_->f(),
          (1 / (double) N) * tf,
          ode_->state_space_->X()(all, k),
          ode_->control_space_->U()(all, k));
      g_sum = g_sum + sum1(MX::abs(ode_->state_space_->X()(all, k + 1) - x_next));
      g_max = mmax(MX::vertcat({g_max, MX::abs(ode_->state_space_->X()(all, k + 1) - x_next)}));
    }
    J_max_ = Function("J_max", {ode_->state_space_->X(), ode_->control_space_->U(), tf}, {g_max*N/tf}); //TODO: put it back to g_max
    J_real_ = Function("J_real", {ode_->state_space_->X(), ode_->control_space_->U(), tf}, {g_sum});//TODO: it wa g_sum  before
  }

private:
  std::shared_ptr<Transcription> clone(const std::shared_ptr<OdeApprox> ode_approx) const override {
    return std::make_shared<MultipleShootingApprox>(ode_approx);
  }
};
}
#endif