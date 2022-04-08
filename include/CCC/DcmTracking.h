/* Author: Masaki Murooka */

#pragma once

#include <map>

#include <Eigen/Core>

#include <CCC/Constants.h>

namespace CCC
{
/** \brief Walking control based on tracking of divergent component of motion (DCM)

    See the following for a detailed formulation.
      - J Englsberger and C Ott. Integration of vertical CoM motion and angular momentum in an extended capture point
   tracking controller for bipedal walking. Humanoids, 2012.
      - J Englsberger, et al. Three-dimensional bipedal walking control using divergent component of motion. IROS, 2013.
 */
class DcmTracking
{
public:
  /** \brief Initial parameter. */
  struct InitialParam
  {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Current DCM [m]
    Eigen::Vector2d current_dcm = Eigen::Vector2d::Zero();

    //! Reference ZMP [m]
    Eigen::Vector2d ref_zmp = Eigen::Vector2d::Zero();
  };

public:
  /** \brief Constructor.
      \param com_height height of robot CoM [m]
      \param feedback_gain feedback gain to calculate ZMP
   */
  DcmTracking(double com_height, double feedback_gain = 1.0)
  : omega_(std::sqrt(CCC::constants::g / com_height)), feedback_gain_(feedback_gain)
  {
  }

  /** \brief Plan one step.
    \param time_zmp_list list of pairs of future ZMP and switching time
    \param initial_param initial parameter
    \param current_time current time (i.e., start time of horizon) [sec]
    \returns planned ZMP

    In the referenced paper, the length of time_zmp_list is three.
 */
  Eigen::Vector2d planOnce(const std::map<double, Eigen::Vector2d> & time_zmp_list,
                           const InitialParam & initial_param,
                           double current_time) const;

public:
  //! Feedback gain to calculate control ZMP
  double feedback_gain_ = 0;

protected:
  //! Time constant for inverted pendulum dynamics
  double omega_ = 0;
};
} // namespace CCC
