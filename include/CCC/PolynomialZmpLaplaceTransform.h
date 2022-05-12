/* Author: Masaki Murooka */

#pragma once

#include <functional>

#include <qp_solver_collection/QpSolverCollection.h>

namespace CCC
{
/** \brief Walking control based on analytical solution of DCM using Laplace transform of ZMP trajectory for
   one-dimensional motion.

    See the following for a detailed formulation.
      - T Kamioka, et al. Simultaneous optimization of ZMP and footsteps based on the analytical solution of divergent
   component of motion. ICRA, 2018.
 */
class PolynomialZmpLaplaceTransform1d
{
  friend class PolynomialZmpLaplaceTransform;

public:
  /** \brief Reference data. */
  struct RefData
  {
    //! Time list of contact switching
    std::vector<double> contact_switch_time_list;

    //! Capture point of last contact switching
    double last_capture_point = 0;

    //! Function to return min/max limits of ZMP [m]
    std::function<std::array<double, 2>(double)> zmp_limits_func;
  };

  /** \brief Initial parameter.

      Initial parameter is capture point only.
  */
  using InitialParam = double;

public:
  /** \brief Constructor.
      \param com_height height of robot CoM [m]
      \param qp_solver_type QP solver type
   */
  PolynomialZmpLaplaceTransform1d(
      double com_height,
      QpSolverCollection::QpSolverType qp_solver_type = QpSolverCollection::QpSolverType::Any);

  /** \brief Plan one step.
      \param ref_data reference data
      \param initial_param initial parameter (current capture point)
      \param current_time current time (i.e., start time of horizon) [sec]
      \returns planned ZMP
   */
  double planOnce(const RefData & ref_data, const InitialParam & initial_param, double current_time);

protected:
  /** \brief Equation (9) \f$\xi^{mon}_{p,j}(t)\f$ in the paper. */
  double xiMonP(int j, double t) const;

protected:
  //! QP solver
  std::shared_ptr<QpSolverCollection::QpSolver> qp_solver_;

  //! QP coefficients
  QpSolverCollection::QpCoeff qp_coeff_;

  //! Time constant for inverted pendulum dynamics
  double omega_ = 0;

  //! Polynomial degree of ZMP trajectory
  int polynomial_degree_ = 9;

  //! Matrix dependent on polynomial degree
  Eigen::MatrixXd H_bar_ = Eigen::MatrixXd::Zero(0, 0);
};

/** \brief Walking control based on analytical solution of DCM using Laplace transform of ZMP trajectory.

    See the following for a detailed formulation.
      - T Kamioka, et al. Simultaneous optimization of ZMP and footsteps based on the analytical solution of divergent
   component of motion. ICRA, 2018.
 */
class PolynomialZmpLaplaceTransform
{
public:
  /** \brief Reference data.

      \todo It is assumed that ZMP limits are independent for the x and y components. This assumption is not valid if
     the foot is placed diagonally during the single-support phase or if the feet are not aligned during the
     double-support phase.

      \todo Footstep online planning is not supported.
   */
  struct RefData
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Time list of contact switching
    std::vector<double> contact_switch_time_list;

    //! Capture point of last contact switching
    Eigen::Vector2d last_capture_point = Eigen::Vector2d::Zero();

    //! Function to return min/max limits of ZMP [m]
    std::function<std::array<Eigen::Vector2d, 2>(double)> zmp_limits_func;
  };

  /** \brief Initial parameter.

      Initial parameter is capture point only.
  */
  using InitialParam = Eigen::Vector2d;

public:
  /** \brief Constructor.
      \param com_height height of robot CoM [m]
      \param qp_solver_type QP solver type
   */
  PolynomialZmpLaplaceTransform(double com_height,
                                QpSolverCollection::QpSolverType qp_solver_type = QpSolverCollection::QpSolverType::Any)
  : pzlt_1d_(std::make_shared<PolynomialZmpLaplaceTransform1d>(com_height, qp_solver_type))
  {
  }

  /** \brief Plan one step.
      \param ref_data reference data
      \param initial_param initial parameter (current capture point)
      \param current_time current time (i.e., start time of horizon) [sec]
      \returns planned ZMP
   */
  Eigen::Vector2d planOnce(const RefData & ref_data, const InitialParam & initial_param, double current_time);

protected:
  //! One-dimensional linear MPC
  std::shared_ptr<PolynomialZmpLaplaceTransform1d> pzlt_1d_;
};
} // namespace CCC
