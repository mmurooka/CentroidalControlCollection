/* Author: Masaki Murooka */

#include <gtest/gtest.h>

#include <chrono>
#include <fstream>
#include <iostream>

#include <qp_solver_collection/QpSolverCollection.h>

#include <CCC/Constants.h>
#include <CCC/DdpZmp.h>

#include "FootstepManager.h"
#include "SimModels.h"

TEST(TestDdpZmp, Test1)
{
  double horizon_duration = 2.0; // [sec]
  double horizon_dt = 0.02; // [sec]
  int horizon_steps = static_cast<int>(horizon_duration / horizon_dt);
  double sim_dt = 0.005; // [sec]
  double mass = 100.0; // [Kg]
  double ref_com_height = 1.0; // [m]
  std::vector<double> disturb_time_list = {4.5, 8.5}; // [sec]
  Eigen::Vector2d disturb_impulse_per_mass = Eigen::Vector2d(0.05, 0.05); // [m/s]

  // Setup DDP
  std::vector<double> computation_duration_list;
  CCC::DdpZmp ddp(mass, horizon_dt, horizon_steps);
  ddp.ddp_solver_->config().max_iter = 3;

  // Setup footstep
  FootstepManager footstep_manager;
  double transit_duration = 0.2; // [sec]
  double swing_duration = 0.8; // [sec]
  footstep_manager.appendFootstep(
      Footstep(Foot::Left, Eigen::Vector2d(0.2, 0.1), 2.0, transit_duration, swing_duration));
  footstep_manager.appendFootstep(
      Footstep(Foot::Right, Eigen::Vector2d(0.4, -0.1), 3.0, transit_duration, swing_duration));
  footstep_manager.appendFootstep(
      Footstep(Foot::Left, Eigen::Vector2d(0.6, 0.1), 4.0, transit_duration, swing_duration));
  footstep_manager.appendFootstep(
      Footstep(Foot::Right, Eigen::Vector2d(0.8, -0.1), 5.0, transit_duration, swing_duration));
  footstep_manager.appendFootstep(
      Footstep(Foot::Left, Eigen::Vector2d(0.6, 0.1), 6.0, transit_duration, swing_duration));
  footstep_manager.appendFootstep(
      Footstep(Foot::Right, Eigen::Vector2d(0.6, -0.1), 7.0, transit_duration, swing_duration));
  std::function<CCC::DdpZmp::RefData(double)> ref_data_func = [&](double t) {
    CCC::DdpZmp::RefData ref_data;
    ref_data.zmp = footstep_manager.refZmp(t);
    ref_data.com_z = ref_com_height;
    return ref_data;
  };

  // Setup simulation
  ComZmpSim3d sim(mass, sim_dt);
  sim.state_.z[0] = ref_com_height;

  // Setup dump file
  std::string file_path = "/tmp/TestDdpZmp.txt";
  std::ofstream ofs(file_path);
  ofs << "time com_pos_x com_pos_y com_pos_z planned_zmp_x planned_zmp_y planned_force_z ref_zmp_x ref_zmp_y ref_com_z "
         "capture_point_x capture_point_y capture_point_ds_end_x capture_point_ds_end_y "
         "capture_point_ds_end_constrained_x capture_point_ds_end_constrained_y ddp_iter computation_time"
      << std::endl;

  // Setup control loop
  CCC::DdpZmp::PlannedData planned_data;

  // Run control loop
  constexpr double end_time = 10.0; // [sec]
  double t = 0;
  bool first_iter = true;
  while(t < end_time)
  {
    // Plan
    auto start_time = std::chrono::system_clock::now();
    footstep_manager.update(t);
    CCC::DdpZmp::InitialParam initial_param;
    initial_param.pos = sim.state_.pos();
    initial_param.vel = sim.state_.vel();
    if(first_iter)
    {
      first_iter = false;
      initial_param.u_list.assign(horizon_steps,
                                  CCC::DdpZmp::DdpProblem::InputDimVector(sim.state_.pos().x(), sim.state_.pos().y(),
                                                                          mass * CCC::constants::g));
    }
    else
    {
      initial_param.u_list = ddp.ddp_solver_->controlData().u_list;
    }
    planned_data = ddp.planOnce(ref_data_func, initial_param, t);
    computation_duration_list.push_back(
        1e3
        * std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::system_clock::now() - start_time)
              .count());

    ////////////////
    // Calculate capture point
    Eigen::Vector2d capture_point =
        (initial_param.pos + std::sqrt(sim.state_.pos().z() / CCC::constants::g) * initial_param.vel).head<2>();

    // Calculate ds_end_time
    // End time of double-support phase (i.e., next swing_start_time)
    double ds_end_time = t + 1.0;
    if(footstep_manager.footstep_list_.size() >= 1 && t < footstep_manager.footstep_list_[0].swing_start_time)
    {
      ds_end_time = footstep_manager.footstep_list_[0].swing_start_time;
    }
    else if(footstep_manager.footstep_list_.size() >= 2)
    {
      ds_end_time = footstep_manager.footstep_list_[1].swing_start_time;
    }

    // Calculate capture_point_ds_end
    int idx_ds_end = std::clamp(static_cast<int>((ds_end_time - t) / horizon_dt), 0, horizon_steps - 1);
    CCC::DdpZmp::DdpProblem::StateDimVector x_ds_end = ddp.ddp_solver_->controlData().x_list[idx_ds_end];
    Eigen::Vector2d capture_point_ds_end =
        Eigen::Vector2d(x_ds_end[0] + std::sqrt(sim.state_.pos().z() / CCC::constants::g) * x_ds_end[1],
                        x_ds_end[2] + std::sqrt(sim.state_.pos().z() / CCC::constants::g) * x_ds_end[3]);

    Eigen::Vector2d capture_point_ds_end_constrained = Eigen::Vector2d::Zero();
    if(footstep_manager.footstep_list_.size() >= 2 && footstep_manager.footstep_list_[0].swing_start_time <= t
       && t < footstep_manager.footstep_list_[0].swing_end_time)
    {
      // Calculate variables
      double omega = std::sqrt(CCC::constants::g / sim.state_.pos().z());
      // End time of single-support phase (i.e., current swing_end_time)
      double ss_end_time = footstep_manager.footstep_list_[0].swing_end_time;
      double exp_ss = std::exp(omega * (ss_end_time - t));
      double exp_ds = std::exp(omega * (ds_end_time - ss_end_time));

      // Calculate coefficients of least-square
      Eigen::Matrix<double, 2, 4> ls_coeff = Eigen::Matrix<double, 2, 4>::Zero();
      ls_coeff.leftCols<2>().diagonal().setConstant(exp_ds * (1 - exp_ss));
      ls_coeff.rightCols<2>().diagonal().setConstant(1 - exp_ds);
      Eigen::Vector2d ls_const = exp_ss * exp_ds * capture_point - capture_point_ds_end;

      // Set QP coefficients
      QpSolverCollection::QpCoeff qp_coeff;
      qp_coeff.setup(4, 0, 0);
      qp_coeff.obj_mat_ = ls_coeff.transpose() * ls_coeff;
      qp_coeff.obj_vec_ = ls_coeff.transpose() * ls_const;
      const Eigen::Vector2d & foot_size = footstep_manager.foot_size_;
      Footstance tmp_footstance = footstep_manager.footstance_;
      tmp_footstance.erase(footstep_manager.footstep_list_[0].foot);
      auto support_region_ss = tmp_footstance.supportRegion();
      qp_coeff.x_min_.head<2>() = support_region_ss[0] - 0.5 * foot_size;
      qp_coeff.x_max_.head<2>() = support_region_ss[1] + 0.5 * foot_size;
      tmp_footstance.emplace(footstep_manager.footstep_list_[0].foot, footstep_manager.footstep_list_[0].pos);
      auto support_region_ds = tmp_footstance.supportRegion();
      qp_coeff.x_min_.tail<2>() = support_region_ds[0] - 0.5 * foot_size;
      qp_coeff.x_max_.tail<2>() = support_region_ds[1] + 0.5 * foot_size;

      // Solve QP
      auto qp_solver = QpSolverCollection::allocateQpSolver(QpSolverCollection::QpSolverType::QLD);
      Eigen::VectorXd qp_sol = qp_solver->solve(qp_coeff);
      Eigen::Vector2d zmp_ss = qp_sol.head<2>();
      Eigen::Vector2d zmp_ds = qp_sol.tail<2>();

      // Calculate capture point with support region constraints
      Eigen::Vector2d capture_point_ss_end_constrained = exp_ss * (capture_point - zmp_ss) + zmp_ss;
      capture_point_ds_end_constrained = exp_ds * (capture_point_ss_end_constrained - zmp_ds) + zmp_ds;
    }
    ////////////////

    // Dump
    const auto & ref_data = ref_data_func(t);
    ofs << t << " " << sim.state_.pos().transpose() << " " << planned_data.zmp.transpose() << " "
        << planned_data.force_z << " " << ref_data.zmp.transpose() << " " << ref_com_height << " "
        << capture_point.transpose() << " " << capture_point_ds_end.transpose() << " "
        << capture_point_ds_end_constrained.transpose() << " " << ddp.ddp_solver_->traceDataList().back().iter << " "
        << computation_duration_list.back() << std::endl;

    // Check
    EXPECT_LT((planned_data.zmp - ref_data.zmp).norm(), 0.1); // [m]
    EXPECT_LT(std::abs(sim.state_.pos().z() - ref_com_height), 0.1); // [m]

    // Simulate
    t += sim_dt;
    ComZmpSim3d::Input sim_input;
    sim_input.zmp = planned_data.zmp;
    sim_input.force_z = planned_data.force_z;
    sim.update(sim_input);

    // Add disturbance
    for(double disturb_time : disturb_time_list)
    {
      if(disturb_time <= t && t < disturb_time + sim_dt)
      {
        sim.addDisturb(disturb_impulse_per_mass);
        break;
      }
    }
  }

  // Final check
  const auto & ref_data = ref_data_func(t);
  EXPECT_LT((planned_data.zmp - ref_data.zmp).norm(), 1e-2);
  EXPECT_LT(std::abs(sim.state_.pos().z() - ref_com_height), 1e-2);
  EXPECT_LT((sim.state_.pos().head<2>() - ref_data.zmp).norm(), 1e-2);
  EXPECT_LT(sim.state_.vel().norm(), 1e-2);

  Eigen::Map<Eigen::VectorXd> computation_duration_vec(computation_duration_list.data(),
                                                       computation_duration_list.size());
  std::cout << "Computation time per control cycle:\n"
            << "  mean: " << computation_duration_vec.mean() << " [ms], stdev: "
            << std::sqrt((computation_duration_vec.array() - computation_duration_vec.mean()).square().mean())
            << " [ms], max: " << computation_duration_vec.maxCoeff() << " [ms]" << std::endl;

  std::cout << "Run the following commands in gnuplot:\n"
            << "  set key autotitle columnhead\n"
            << "  set key noenhanced\n"
            << "  plot \"" << file_path
            << "\" u 1:2 w lp, \"\" u 1:5 w lp, \"\" u 1:8 w l lw 2, \"\" u 1:11 w lp, \"\" u 1:13 w lp # x\n"
            << "  plot \"" << file_path
            << "\" u 1:3 w lp, \"\" u 1:6 w lp, \"\" u 1:9 w l lw 2, \"\" u 1:12 w lp, \"\" u 1:14 w lp # y\n"
            << "  plot \"" << file_path << "\" u 1:4 w lp, \"\" u 1:10 w l lw 2 # z\n"
            << "  plot \"" << file_path << "\" u 1:13 w lp # ddp_iter\n"
            << "  plot \"" << file_path << "\" u 1:14 w lp # computation_time\n";
}

TEST(TestDdpZmp, CheckDerivatives)
{
  constexpr double deriv_eps = 1e-6;

  double horizon_dt = 0.005; // [sec]
  double mass = 100.0; // [Kg]
  CCC::DdpZmp::WeightParam weight_param;
  auto ddp_problem = std::make_shared<CCC::DdpZmp::DdpProblem>(horizon_dt, mass, weight_param);

  std::function<CCC::DdpZmp::RefData(double)> ref_data_func = [](double t) {
    CCC::DdpZmp::RefData ref_data;
    ref_data.zmp << 0.1, -0.2; // [m]
    ref_data.com_z = 1.0; // [m]
    return ref_data;
  };
  ddp_problem->setRefDataFunc(ref_data_func);

  double t = 0;
  CCC::DdpZmp::DdpProblem::StateDimVector x;
  x << 1.0, 0.1, -2.1, -0.5, 1.1, 0.5; // [m], [m/s]
  CCC::DdpZmp::DdpProblem::InputDimVector u;
  u << 1.0, -2.0, 1000.0; // [m], [N]

  // Check state_eq_deriv
  {
    CCC::DdpZmp::DdpProblem::StateStateDimMatrix state_eq_deriv_x_analytical;
    CCC::DdpZmp::DdpProblem::StateInputDimMatrix state_eq_deriv_u_analytical;
    ddp_problem->calcStateEqDeriv(t, x, u, state_eq_deriv_x_analytical, state_eq_deriv_u_analytical);

    CCC::DdpZmp::DdpProblem::StateStateDimMatrix state_eq_deriv_x_numerical;
    CCC::DdpZmp::DdpProblem::StateInputDimMatrix state_eq_deriv_u_numerical;
    for(int i = 0; i < ddp_problem->stateDim(); i++)
    {
      state_eq_deriv_x_numerical.col(i) =
          (ddp_problem->stateEq(t, x + deriv_eps * CCC::DdpZmp::DdpProblem::StateDimVector::Unit(i), u)
           - ddp_problem->stateEq(t, x - deriv_eps * CCC::DdpZmp::DdpProblem::StateDimVector::Unit(i), u))
          / (2 * deriv_eps);
    }
    for(int i = 0; i < ddp_problem->inputDim(); i++)
    {
      state_eq_deriv_u_numerical.col(i) =
          (ddp_problem->stateEq(t, x, u + deriv_eps * CCC::DdpZmp::DdpProblem::InputDimVector::Unit(i))
           - ddp_problem->stateEq(t, x, u - deriv_eps * CCC::DdpZmp::DdpProblem::InputDimVector::Unit(i)))
          / (2 * deriv_eps);
    }

    EXPECT_LT((state_eq_deriv_x_analytical - state_eq_deriv_x_numerical).norm(), 1e-6);
    EXPECT_LT((state_eq_deriv_u_analytical - state_eq_deriv_u_numerical).norm(), 1e-6);
  }

  // Check running_cost_deriv
  {
    CCC::DdpZmp::DdpProblem::StateDimVector running_cost_deriv_x_analytical;
    CCC::DdpZmp::DdpProblem::InputDimVector running_cost_deriv_u_analytical;
    ddp_problem->calcRunningCostDeriv(t, x, u, running_cost_deriv_x_analytical, running_cost_deriv_u_analytical);

    CCC::DdpZmp::DdpProblem::StateDimVector running_cost_deriv_x_numerical;
    CCC::DdpZmp::DdpProblem::InputDimVector running_cost_deriv_u_numerical;
    for(int i = 0; i < ddp_problem->stateDim(); i++)
    {
      running_cost_deriv_x_numerical[i] =
          (ddp_problem->runningCost(t, x + deriv_eps * CCC::DdpZmp::DdpProblem::StateDimVector::Unit(i), u)
           - ddp_problem->runningCost(t, x - deriv_eps * CCC::DdpZmp::DdpProblem::StateDimVector::Unit(i), u))
          / (2 * deriv_eps);
    }
    for(int i = 0; i < ddp_problem->inputDim(); i++)
    {
      running_cost_deriv_u_numerical[i] =
          (ddp_problem->runningCost(t, x, u + deriv_eps * CCC::DdpZmp::DdpProblem::InputDimVector::Unit(i))
           - ddp_problem->runningCost(t, x, u - deriv_eps * CCC::DdpZmp::DdpProblem::InputDimVector::Unit(i)))
          / (2 * deriv_eps);
    }

    EXPECT_LT((running_cost_deriv_x_analytical - running_cost_deriv_x_numerical).norm(), 1e-6);
    EXPECT_LT((running_cost_deriv_u_analytical - running_cost_deriv_u_numerical).norm(), 1e-6);
  }

  // Check terminal_cost_deriv
  {
    CCC::DdpZmp::DdpProblem::StateDimVector terminal_cost_deriv_x_analytical;
    ddp_problem->calcTerminalCostDeriv(t, x, terminal_cost_deriv_x_analytical);

    CCC::DdpZmp::DdpProblem::StateDimVector terminal_cost_deriv_x_numerical;
    for(int i = 0; i < ddp_problem->stateDim(); i++)
    {
      terminal_cost_deriv_x_numerical[i] =
          (ddp_problem->terminalCost(t, x + deriv_eps * CCC::DdpZmp::DdpProblem::StateDimVector::Unit(i))
           - ddp_problem->terminalCost(t, x - deriv_eps * CCC::DdpZmp::DdpProblem::StateDimVector::Unit(i)))
          / (2 * deriv_eps);
    }

    EXPECT_LT((terminal_cost_deriv_x_analytical - terminal_cost_deriv_x_numerical).norm(), 1e-6);
  }
}

int main(int argc, char ** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
