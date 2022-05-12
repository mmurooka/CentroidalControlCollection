/* Author: Masaki Murooka */

#include <CCC/Constants.h>
#include <CCC/PolynomialZmpLaplaceTransform.h>

using namespace CCC;

PolynomialZmpLaplaceTransform1d::PolynomialZmpLaplaceTransform1d(double com_height,
                                                                 QpSolverCollection::QpSolverType qp_solver_type)
: omega_(std::sqrt(constants::g / com_height))
{
  qp_solver_ = QpSolverCollection::allocateQpSolver(qp_solver_type);
}

double PolynomialZmpLaplaceTransform1d::planOnce(const RefData & ref_data,
                                                 const InitialParam & initial_param,
                                                 double current_time)
{
  // Set segment_time_list relative to current_time
  std::vector<double> segment_time_list = {0.0};
  for(const auto & time : ref_data.contact_switch_time_list)
  {
    segment_time_list.push_back(time - current_time);
  }

  // Set variables
  int e = segment_time_list.size() - 1;
  int n = polynomial_degree_ + 1; // Increment by 1 to include a constant term
  int dim_var = e * n;

  // Equation (15) in the paper
  Eigen::MatrixXd A_cyc(1, dim_var);
  Eigen::VectorXd b_cyc(1);
  for(int i = 0; i < e; i++)
  {
    double t_i = segment_time_list[i];
    double T_i = segment_time_list[i + 1] - segment_time_list[i];
    for(int j = 0; j < n; j++)
    {
      A_cyc(0, i * n + j) = std::exp(-1 * omega_ * t_i) * std::pow(T_i, -1 * j) * xiMonP(j, T_i);
    }
  }
  b_cyc << -1 * initial_param + std::exp(-1 * omega_ * segment_time_list.back()) * ref_data.last_capture_point;

  // todo: Equation (20)

  // Equation (21) in the paper
  Eigen::MatrixXd A_continuity_pos = Eigen::MatrixXd::Zero(e - 1, dim_var);
  for(int i = 0; i < e - 1; i++)
  {
    A_continuity_pos.block(i, i * n, 1, n).setConstant(1);
    A_continuity_pos(i, (i + 1) * n) = -1;
  }

  // Equation (22) in the paper
  Eigen::MatrixXd A_continuity_vel = Eigen::MatrixXd::Zero(e - 1, dim_var);
  for(int i = 0; i < e - 1; i++)
  {
    double T_i_prev = segment_time_list[i + 1] - segment_time_list[i];
    double T_i = segment_time_list[i + 2] - segment_time_list[i + 1];
    for(int j = 1; j < n; j++)
    {
      A_continuity_vel(i, i * n + j) = static_cast<double>(j) / T_i_prev;
    }
    A_continuity_vel(i, (i + 1) * n + 1) = -1 / T_i;
  }

  // Equation (23) in the paper
  Eigen::MatrixXd A_continuity_acc = Eigen::MatrixXd::Zero(e - 1, dim_var);
  for(int i = 0; i < e - 1; i++)
  {
    double T_i_prev = segment_time_list[i + 1] - segment_time_list[i];
    double T_i = segment_time_list[i + 2] - segment_time_list[i + 1];
    for(int j = 2; j < n; j++)
    {
      A_continuity_acc(i, i * n + j) = static_cast<double>(j) * static_cast<double>(j - 1) / std::pow(T_i_prev, 2);
    }
    A_continuity_acc(i, (i + 1) * n + 2) = -2 / std::pow(T_i, 2);
  }

  // todo: Equation (24)

  // Equation (25) in the paper
  if(H_bar_.rows() != n || H_bar_.cols() != n)
  {
    H_bar_.setZero(n, n);
    for(int r = 0; r < n; r++)
    {
      if(r < 2)
      {
        continue;
      }
      for(int c = 0; c < n; c++)
      {
        if(c < 2)
        {
          continue;
        }
        H_bar_(r, c) = r * (r - 1) * c * (c - 1) / (r + c - 3);
      }
    }
  }

  // Set QP coefficients
  int dim_eq = A_cyc.rows() + A_continuity_pos.rows() + A_continuity_vel.rows() + A_continuity_acc.rows();
  int dim_ineq = 0;
  qp_coeff_.setup(dim_var, dim_eq, dim_ineq);
  for(int i = 0; i < e; i++)
  {
    double T_i = segment_time_list[i + 1] - segment_time_list[i];
    qp_coeff_.obj_mat_.block(i * n, i * n, n, n) = H_bar_ / std::pow(T_i, 3);
  }
  qp_coeff_.eq_mat_ << A_cyc, A_continuity_pos, A_continuity_vel, A_continuity_acc;
  qp_coeff_.eq_vec_.segment(0, b_cyc.size()) = b_cyc;

  // Solve QP
  Eigen::VectorXd qp_sol = qp_solver_->solve(qp_coeff_);

  // Calculate ZMP
  double zmp = qp_sol[0]; // Current ZMP is the constant term of the first segment

  return zmp;
}

namespace
{
int permutation(int n, int k)
{
  int ret = 1;
  for(int i = 0; i < k; i++)
  {
    ret *= n - i;
  }
  return ret;
}

int factorial(int n)
{
  int ret = 1;
  for(int i = 1; i <= n; i++)
  {
    ret *= i;
  }
  return ret;
}
} // namespace

double PolynomialZmpLaplaceTransform1d::xiMonP(int j, double t) const
{
  // Equation (9) in the paper
  double sum = 0;
  for(int l = 0; l < j; l++)
  {
    sum += permutation(j, j - l) * std::pow(t, l) * std::pow(omega_, -1 * (j - l));
  }
  return -1 * (factorial(j) * std::pow(omega_, -1 * j) - std::exp(-1 * omega_ * t) * sum);
}

Eigen::Vector2d PolynomialZmpLaplaceTransform::planOnce(const RefData & ref_data,
                                                        const InitialParam & initial_param,
                                                        double current_time)
{
  // Calculate ZMP
  Eigen::Vector2d planned_zmp;

  PolynomialZmpLaplaceTransform1d::RefData ref_data_x;
  ref_data_x.contact_switch_time_list = ref_data.contact_switch_time_list;
  ref_data_x.last_capture_point = ref_data.last_capture_point.x();
  ref_data_x.zmp_limits_func = [&](double t) -> std::array<double, 2> {
    const auto & zmp_limits = ref_data.zmp_limits_func(t);
    return std::array<double, 2>{zmp_limits[0].x(), zmp_limits[1].x()};
  };
  PolynomialZmpLaplaceTransform1d::InitialParam initial_param_x = initial_param.x();
  planned_zmp.x() = pzlt_1d_->planOnce(ref_data_x, initial_param_x, current_time);

  PolynomialZmpLaplaceTransform1d::RefData ref_data_y;
  ref_data_y.contact_switch_time_list = ref_data.contact_switch_time_list;
  ref_data_y.last_capture_point = ref_data.last_capture_point.y();
  ref_data_y.zmp_limits_func = [&](double t) -> std::array<double, 2> {
    const auto & zmp_limits = ref_data.zmp_limits_func(t);
    return std::array<double, 2>{zmp_limits[0].y(), zmp_limits[1].y()};
  };
  PolynomialZmpLaplaceTransform1d::InitialParam initial_param_y = initial_param.y();
  planned_zmp.y() = pzlt_1d_->planOnce(ref_data_y, initial_param_y, current_time);

  return planned_zmp;
}
