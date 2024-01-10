#define EIGEN_RUNTIME_NO_MALLOC
#define eigen_assert(x) \
  if(!(x)) { throw std::runtime_error("Fail: " #x); }

#include "../include/pendulum_feasibility_solver/feasibility_solver.h"
#include "../include/pendulum_feasibility_solver/polygons.h"

void feasibility_solver::timings_constraints(Eigen::Ref<Eigen::MatrixXd> A_out,
                                             Eigen::Ref<Eigen::VectorXd> b_out,
                                             const int NStepsTimings)
{
  Eigen::internal::set_is_malloc_allowed(false);
  const Eigen::Index NVariables = A_out.cols();
  const Eigen::Index N_mu = NVariables - N_.rows(); // N var minus the slack variables

  const double mu_ss_max = exp(-eta_ * t_ss_range_.y());
  const double mu_ss_min = exp(-eta_ * t_ss_range_.x());

  const double mu_ds_max = exp(-eta_ * t_ds_range_.y());
  const double mu_ds_min = exp(-eta_ * t_ds_range_.x());

  const double mu_s_max = exp(-eta_ * t_s_range_.y());
  const double mu_s_min = exp(-eta_ * t_s_range_.x());

  // Number of constraints:
  // - mu must decrease over the horizon (N_mu - 1 constraints)
  // - Step must be bounded and Sg supp must be bounded (5 * NStepsTimings constraints)
  // - Mu must be positive (N_mu constraints)
  const Eigen::Index N_cstr = (N_mu - 1) + Eigen::Index(5) * NStepsTimings + N_mu;
  assert(A_out.rows() == N_cstr && A_out.cols() == NVariables);
  assert(b_out.rows() == A_out.rows());

  // Resize and reset A_out and b_out so they have the expected size
  A_out.setZero();
  b_out.setZero();

  // All mu must be decreasing in the horizon
  A_out.block(0, 1, N_mu - 1, N_mu - 1).diagonal().setConstant(1.0);
  A_out.block(0, 0, N_mu - 1, N_mu - 1).diagonal().setConstant(-1.0);
  A_out(0, 0) *= exp(-eta_ * delta_);

  Eigen::Index cstrRow = N_mu - 1;
  for(Eigen::Index i = 0; i < NStepsTimings; ++i)
  {
    // Step must be bounded
    // Sg supp must be bounded
    auto A_min = A_out.middleRows(cstrRow, 3);
    auto b_min = b_out.segment(cstrRow, 3);
    auto A_max = A_out.middleRows(cstrRow + 3, 2);
    auto b_max = b_out.segment(cstrRow + 3, 2);
    cstrRow += 5;
    A_min.block(0, (N_ds_ + 1) * i + N_ds_, 1, 2) = Eigen::RowVector2d{-1, 1};

    // Sg supp time cstr
    A_min(0, (N_ds_ + 1) * i + N_ds_) *= mu_ss_min;

    // Step time cstr
    A_min(1, (i + 1) * (N_ds_ + 1)) = 1;

    // Ds cstr
    A_min(2, (N_ds_ + 1) * i) = -mu_ds_min;
    A_min(2, (N_ds_ + 1) * i + N_ds_) = 1;

    // Ts <= ts_max
    A_max(0, (i + 1) * (N_ds_ + 1)) = -1;

    // Ts - Tds <= tss_max
    A_max(1, (i + 1) * (N_ds_ + 1)) = -1;
    A_max(1, (N_ds_ + 1) * i + N_ds_) = mu_ss_max;

    if(i != 0)
    {
      // Step time cstr
      A_min(1, i * (N_ds_ + 1)) = -mu_s_min;

      A_max(0, i * (N_ds_ + 1)) = mu_s_max;
    }

    else
    {
      if(!doubleSupport_)
      {
        // Sg supp time cstr
        // t_step >= t_ds + min_ss
        A_min(0, N_ds_) = 0;
        b_min(0) = mu_ss_min * exp(-eta_ * tLift_);

        // Ds cstr
        A_min.row(2).setZero();
      }
      else
      {
        A_min(2, 0) = 0;
        b_min(2) = mu_ds_min;
      }

      // t_step_0 >= t_s_min
      b_min(1) = mu_s_min;

      // Ts <= ts_max
      b_max(0) = -mu_s_max;
    }
  }
  // mu must be positive
  A_out.bottomRows(N_mu).diagonal().setConstant(-1.0);
  b_out.tail(N_mu).setConstant(-exp(-eta_ * 10));
  Eigen::internal::set_is_malloc_allowed(true);
}

void feasibility_solver::build_time_feasibility_matrix(Eigen::Ref<Eigen::MatrixXd> A_f,
                                                       Eigen::Ref<Eigen::VectorXd> b_f,
                                                       const sva::PTransformd & X_0_supp,
                                                       const sva::PTransformd & X_0_swg)
{
  Eigen::internal::set_is_malloc_allowed(false);
  const int NStepsTimings = N_steps;
  const int N_slack = static_cast<int>(N_.rows());
  const int N_variables = NStepsTimings * (N_ds_ + 1) + N_tdsLast + N_slack;

  const Eigen::Vector2d & P_supportFoot_0 = X_0_supp.translation().segment(0, 2);
  const Eigen::Vector2d & P_swingFoot_0 = X_0_swg.translation().segment(0, 2);

  // Variables are organised as such  : timings then slack variables
  // step i occur at indx (N_ds + 1) * i
  // sg suport i start at indx (N_ds + 1) * i + N_ds

  // DCM must remain inside the feasibility region
  assert(A_f.rows() == N_.rows() && A_f.cols() == N_variables);
  assert(b_f.rows() == A_f.rows());
  A_f.setZero();
  b_f.setZero();
  // A_f * x + b + slck >= N * P_u

  // i = 0
  //  const int j_start = doubleSupport_ ? 0 : N_ds_;
  for(int j = 0; j <= N_ds_; j++)
  {
    double alpha_j = static_cast<double>(j) / static_cast<double>(N_ds_);

    Eigen::Vector4d Oij = offsetCstrZMP_ + N_ * P_supportFoot_0;
    if(doubleSupport_ && j < N_ds_) { Oij = offsetCstrZMPDblInit_ + N_ * 0.5 * (P_supportFoot_0 + P_swingFoot_0); }

    A_f.block(0, j, 4, 1) += Oij;
    A_f.block(0, j + 1, 4, 1) -= Oij;
  }

  // Remainings
  for(int i = 1; i <= NStepsTimings; i++)
  {

    for(int j = 0; j <= (i != NStepsTimings ? N_ds_ : N_tdsLast - 1); j++)
    {
      const double alpha_j = static_cast<double>(j + 1) / static_cast<double>(N_ds_ + 1);

      Eigen::Vector2d xStepI = xStep_.segment(2 * (i - 1), 2);
      Eigen::Vector4d O_ij = offsetCstrZMP_ + N_ * (alpha_j * xStep_.segment<2>(2 * (i - 1)));

      if(i > 1) { O_ij += N_ * ((1 - alpha_j) * xStep_.segment<2>(2 * (i - 2))); }
      else { O_ij += N_ * ((1 - alpha_j) * P_supportFoot_0); }

      A_f.block(0, i * (N_ds_ + 1) + j, 4, 1) += O_ij;
      if(!(i == NStepsTimings && j == N_tdsLast - 1)) { A_f.block(0, i * (N_ds_ + 1) + j + 1, 4, 1) -= O_ij; }
    }
  }
  b_f *= kappa_;
  b_f -= N_ * gamma_ * exp(-eta_ * t_);
  A_f *= kappa_;
  Eigen::internal::set_is_malloc_allowed(true);
}

bool feasibility_solver::solve_timings(const std::vector<double> & refTimings, const double & refTds)
{
  const Eigen::DenseIndex NStepsTimings = N_steps;
  const Eigen::DenseIndex N_slack = static_cast<int>(N_.rows());
  const Eigen::DenseIndex N_mu = NStepsTimings * (N_ds_ + 1) + N_tdsLast;
  const Eigen::DenseIndex N_variables = N_mu + N_slack;

  // Variables are organised as such  : timings then slack variables
  // step i occur at indx (N_ds + 1) * i
  // sg suport i start at indx (N_ds + 1) * i + N_ds

  const Eigen::Index N_dcm_cstr = N_.rows();
  const Eigen::Index N_timings_cstr = (N_mu - 1) + Eigen::Index(5) * NStepsTimings + N_mu;
  const Eigen::Index N_cstr = N_dcm_cstr + N_timings_cstr;
  if(A_ineq_buffer_.size() < N_cstr * N_variables) { A_ineq_buffer_.resize(N_cstr * N_variables); }
  auto A_ineq_ = Eigen::Map<Eigen::MatrixXd>(A_ineq_buffer_.data(), N_cstr, N_variables);
  if(b_ineq_buffer_.size() < A_ineq_.rows()) { b_ineq_buffer_.resize(A_ineq_.rows()); }
  auto b_ineq_ = Eigen::Map<Eigen::VectorXd>(b_ineq_buffer_.data(), A_ineq_.rows());

  // DCM must remain inside the feasibility region
  A_f_.resize(N_.rows(), N_variables);
  b_f_.resize(A_f_.rows());
  // A_f * x + b + slck >= N * P_u

  build_time_feasibility_matrix(A_f_, b_f_, X_0_SupportFoot_, X_0_SwingFoot_);
  A_ineq_.topRows(N_dcm_cstr).noalias() = -A_f_ * exp(eta_ * t_);
  b_ineq_.head(A_f_.rows()).noalias() = b_f_ * exp(eta_ * t_) - N_ * dcm_;

  // Steps timings constraints
  auto A_Tsteps = A_ineq_.bottomRows(N_timings_cstr);
  auto b_Tsteps = b_ineq_.tail(A_Tsteps.rows());
  timings_constraints(A_Tsteps, b_Tsteps, NStepsTimings);

  const int NineqCstr = static_cast<int>(A_ineq_.rows());

  // Slack Variables
  if(Niter_ != 0) { A_ineq_.block(0, N_variables - N_slack, N_slack, N_slack).setIdentity(); }

  if(A_eq_buffer_.size() < 1 * N_variables) { A_eq_buffer_.resize(1 * N_variables); }
  auto A_eq_ = Eigen::Map<Eigen::MatrixXd>(A_eq_buffer_.data(), 1, N_variables);
  A_eq_.setZero();
  A_eq_(0, 0) = 1;
  if(b_eq_buffer_.size() < A_eq_.rows()) { b_eq_buffer_.resize(A_eq_.rows()); }
  auto b_eq_ = Eigen::Map<Eigen::VectorXd>(b_eq_buffer_.data(), A_eq_.rows());
  b_eq_(0) = exp(-eta_ * t_);
  const int NeqCstr = static_cast<int>(A_eq_.rows());

  // Cost function
  Eigen::MatrixXd M_deltaTs = Eigen::MatrixXd::Zero(2 * N_timings, N_variables);
  Eigen::MatrixXd b_deltaTs = Eigen::VectorXd::Zero(M_deltaTs.rows());

  for(int i = 1; i < N_timings; i++)
  {
    M_deltaTs.block(2 * i, (N_ds_ + 1) * i + N_ds_, 1, 2) = Eigen::RowVector2d{-1, 1};

    if(i != 0)
    {
      M_deltaTs(2 * i, (N_ds_ + 1) * i + N_ds_) *= exp(-eta_ * (refTimings[i] - refTimings[i - 1] - refTds_));
      M_deltaTs(2 * i + 1, (N_ds_ + 1) * i + N_ds_) = 1;
      M_deltaTs(2 * i + 1, (N_ds_ + 1) * i) = -exp(-eta_ * (refTds_));
    }
    else
    {
      M_deltaTs(2 * i, N_ds_) *= exp(-eta_ * (refTimings[0] - refTds_));
      if(!doubleSupport_)
      {
        M_deltaTs(2 * i, (N_ds_ + 1) * i + N_ds_) = 0;
        b_deltaTs(2 * i) = exp(-eta_ * (tLift_ + refTimings[0] - refTds_));
      }
    }
  }

  Eigen::MatrixXd M_timings = Eigen::MatrixXd::Identity(N_variables - N_slack, N_variables);
  Eigen::VectorXd b_timings = Eigen::VectorXd::Ones(M_timings.rows());

  Eigen::MatrixXd M_slack = Eigen::MatrixXd::Zero(N_slack, N_variables);
  M_slack.block(0, N_variables - N_slack, N_slack, N_slack) = Eigen::MatrixXd::Identity(N_slack, N_slack);

  // std::cout << "n timings " << N_timings << std::endl;
  double t_im1 = 0;
  for(int j = 0; j <= N_ds_; j++)
  {

    double alpha_j = static_cast<double>(j) / static_cast<double>(N_ds_);
    if(doubleSupport_) { b_timings(j) = exp(-eta_ * (t_ + alpha_j * (refTds - t_))); }
    else { b_timings(j) = exp(-eta_ * t_); }
  }

  for(int i = 1; i <= NStepsTimings; i++)
  {
    t_im1 = refTimings[i - 1];
    for(int j = 0; j <= (i != N_steps ? N_ds_ : N_tdsLast - 1); j++)
    {

      double alpha_j = static_cast<double>(j) / static_cast<double>(N_ds_);
      b_timings((N_ds_ + 1) * i + j) = exp(-eta_ * (t_im1 + alpha_j * (refTds)));
    }
  }

  // Keeping slack only on broken cstr
  Eigen::VectorXd x_init = Eigen::VectorXd::Zero(N_variables);
  // x_init.segment(0,N_variables - N_slack) = b_timings;
  x_init.segment(0, N_variables - N_slack) = xTimings_.segment(0, N_variables - N_slack);

  Eigen::Vector4d feasibilityOffsetInit = exp(eta_ * t_) * (A_f_ * x_init + b_f_);
  Eigen::Vector4d dcm_pose = N_ * dcm_;

  // std::cout << "[Pendulum feasibility solver][Timing solver] xStep " << std::endl << xStep_ << std::endl;
  // std::cout << "[Pendulum feasibility solver][Timing solver] init offset " << std::endl << feasibilityOffsetInit <<
  // std::endl; bool ok = true;
  for(int i = 0; i < 4; i++)
  {
    if(feasibilityOffsetInit(i) < dcm_pose(i) - 1e-5)
    {
      std::cout << "[Pendulum feasibility solver][Timing  solver] "
                << "[iter : " << Niter_ << "] broken cstr on " << i << std::endl;
      std::cout << "offset delta " << feasibilityOffsetInit(i) - dcm_pose(i) << std::endl;
      // ok = false;
    }
  }
  // if(ok && Niter_ == 0)
  // {
  //     return true;
  // }

  if(Q_cost_buffer_.size() < N_variables * N_variables)
  {
    Q_cost_buffer_.resize(N_variables * N_variables);
    c_cost_buffer_.resize(N_variables);
  }
  auto Q_cost_ = Eigen::Map<Eigen::MatrixXd>(Q_cost_buffer_.data(), N_variables, N_variables);
  auto c_cost_ = Eigen::Map<Eigen::VectorXd>(c_cost_buffer_.data(), N_variables);

  Q_cost_ = betaTsteps * M_timings.transpose() * M_timings;
  Q_cost_ += 5e0 * M_deltaTs.transpose() * M_deltaTs;
  Q_cost_ += 1e7 * M_slack.transpose() * M_slack;
  c_cost_ = betaTsteps * (-M_timings.transpose() * b_timings);
  c_cost_ += 5e0 * -M_deltaTs.transpose() * b_deltaTs;

  if(Niter_ <= 2)
  {
    for(int i = 0; i < N_timings; i++)
    {

      Eigen::MatrixXd M_plan = Eigen::MatrixXd::Zero(2, N_variables);
      Eigen::VectorXd b_plan = Eigen::VectorXd::Zero(M_plan.rows());
      M_plan(0, (i) * (N_ds_ + 1)) = 1;
      M_plan(0, (i + 1) * (N_ds_ + 1)) = -1;
      M_plan(1, (N_ds_ + 1) * i + N_ds_) = -10;
      M_plan(1, (N_ds_ + 1) * i) = 10;
      const double steps_error = (optimalSteps_[i].translation() - refSteps_[i].translation()).norm();
      Q_cost_ += 5e0 * steps_error * M_plan.transpose() * M_plan;
      c_cost_ += 5e0 * steps_error * -M_plan.transpose() * b_plan;
    }
  }

  qp_solver_.problem(N_variables, NeqCstr, NineqCstr);
  bool QPsuccess = qp_solver_.solve(Q_cost_, c_cost_, A_eq_, b_eq_, A_ineq_, b_ineq_);

  // std::cout << "init " << Niter_ << std::endl;
  // for (int i = 0 ; i < N_variables - N_slack ; i++)
  // {
  //     std::cout << -log(b_timings(i))/eta_ << std::endl;
  // }
  // std::cout << "dcm" << std::endl << dcm_ << std::endl;

  if(!QPsuccess)
  {
    std::cout << "[Pendulum feasibility solver][Timing solver] "
              << "[iter : " << Niter_ << "] QP Failed" << std::endl;
    return true;
  }

  xTimings_ = qp_solver_.result().segment(0, N_variables - N_slack);
  Eigen::Vector4d feasibilityOffset = exp(eta_ * t_) * (A_f_ * qp_solver_.result() + b_f_);
  // std::cout << "[Pendulum feasibility solver][Timing solver] output offset " << std::endl << feasibilityOffset <<
  // std::endl;
  Polygon feasibilityPolygon = Polygon(N_, feasibilityOffset);
  feasibilityRegion_ = feasibilityPolygon.Get_Polygone_Corners();

  for(int i = 0; i < 4; i++)
  {
    if(feasibilityOffset(i) < dcm_pose(i) - 1e-5)
    {
      std::cout << "[Pendulum feasibility solver][Timing  solver] "
                << "[iter : " << Niter_ << "] QP output has broken cstr on " << i << std::endl;
      // std::cout << "offset delta " << feasibilityOffset(i) - dcm_pose(i) << std::endl;
      std::cout << "slack: " << qp_solver_.result().segment(N_variables - N_slack + i, 1) << std::endl;

      // ok = false;
    }
  }

  // std::cout << "sol " << Niter_ << std::endl;
  // for (int i = 0 ; i < N_variables - N_slack ; i++)
  // {
  //     std::cout << -log(xTimings_(i))/eta_ << std::endl;
  // }
  // std::cout << "Slack Timings: " << QP.result().segment(N_variables - N_slack ,N_slack) << std::endl;

  optimalStepsTimings_.clear();
  optimalDoubleSupportDuration_.clear();
  double tds_i = -log(xTimings_(N_ds_)) / eta_;
  double ts_i = -log(xTimings_((N_ds_ + 1))) / eta_;
  optimalStepsTimings_.push_back(ts_i);
  optimalDoubleSupportDuration_.push_back(tds_i);

  for(int i = 1; i < NStepsTimings; i++)
  {
    tds_i = (-log(xTimings_((N_ds_ + 1) * i + N_ds_)) / eta_) - ts_i;
    ts_i = -log(xTimings_((N_ds_ + 1) * (i + 1))) / eta_;

    optimalStepsTimings_.push_back(ts_i);
    optimalDoubleSupportDuration_.push_back(tds_i);
  }
  if(!doubleSupport_)
  {
    if(optimalDoubleSupportDuration_.size() >= 2)
    {
      optimalDoubleSupportDuration_[0] = optimalDoubleSupportDuration_[1];
    }
    else { optimalDoubleSupportDuration_[0] = refTds; }
  }

  // std::cout << "solution Ts" << std::endl;
  // for (int i = 0 ; i < optimalStepsTimings_.size() ; i++)
  // {
  //     std::cout << optimalStepsTimings_[i] << std::endl;
  // }
  // std::cout << "solution Tds" << std::endl;
  // for (int i = 0 ; i < optimalDoubleSupportDuration_.size() ; i++)
  // {
  //     std::cout << optimalDoubleSupportDuration_[i] << std::endl;
  // }
  return true;
}
