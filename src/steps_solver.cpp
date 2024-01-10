#define EIGEN_RUNTIME_NO_MALLOC
#define eigen_assert(x) \
  if(!(x)) { throw std::runtime_error("Wrong assertion: " #x); }

#include "../include/pendulum_feasibility_solver/feasibility_solver.h"
#include "../include/pendulum_feasibility_solver/polygons.h"

void feasibility_solver::kinematics_constraints(Eigen::Ref<Eigen::MatrixXd> A_out,
                                                Eigen::Ref<Eigen::VectorXd> b_out,
                                                const std::vector<sva::PTransformd> & refSteps)
{
  // Eigen::internal::set_is_malloc_allowed(false);
  const int N_variables = static_cast<int>(A_out.cols());
  const int N_slack = static_cast<int>(N_.rows());
  const int n_steps = static_cast<int>(static_cast<double>(N_variables - N_slack) / 2);
  // const int n_steps = 1;

  assert(n_steps > 0);
  std::vector<Eigen::VectorXd> b_kin_cstr_vec;
  std::vector<Eigen::MatrixX2d> kin_cstr_normals_vec;
  std::vector<Eigen::MatrixX2d> step_cstr_normals_vec;
  std::vector<Eigen::VectorXd> b_step_cstr_vec;
  Eigen::MatrixXd Delta = Eigen::MatrixXd::Identity(2 * n_steps, 2 * n_steps); // Matrix to differentiate two footsteps

  double l = 1;
  if(supportFoot_ == "LeftFoot") { l *= -1; }
  int N_footsteps_kin_cstr = 0;
  for(int i = 0; i < n_steps; i++)
  {
    const double theta_i = rpyFromMat(refSteps[i].rotation()).z();
    sva::PTransformd X_0_step_im1 = X_0_SupportFoot_;
    if(i != 0) { X_0_step_im1 = refSteps[i - 1]; }
    const Eigen::Matrix3d R_Theta_i_0 = X_0_step_im1.rotation().transpose();
    const Eigen::Vector3d offset = R_Theta_i_0 * Eigen::Vector3d{0, l * (feetDistance_ + 0 * stepCstrSize_.y() / 2), 0};

    Rectangle Kinematic_Rectangle = Rectangle(theta_i, stepCstrSize_, offset);

    if(i > 0) { Delta.block(2 * i, 2 * (i - 1), 2, 2) = -Eigen::Matrix2d::Identity(); }
    if(i == 0) { Kinematic_Rectangle = Rectangle(X_0_SupportFoot_, stepCstrSize_, offset); }
    Polygon Kinematic_Poly = Polygon(Kinematic_Rectangle);
    b_kin_cstr_vec.push_back(Kinematic_Poly.offsets());
    kin_cstr_normals_vec.push_back(Kinematic_Poly.normals());

    N_footsteps_kin_cstr += static_cast<int>(kin_cstr_normals_vec.back().rows());
    l *= -1;
  }

  Eigen::MatrixXd foosteps_kin_cstr = Eigen::MatrixXd::Zero(N_footsteps_kin_cstr, 2 * n_steps);
  Eigen::VectorXd b_kin_cstr(N_footsteps_kin_cstr);

  int step = 0;
  int cstr_index = 0;
  for(Eigen::Index i_ineq = 0; i_ineq < static_cast<Eigen::Index>(kin_cstr_normals_vec.size()); i_ineq++)
  {

    Eigen::MatrixX2d ineq = kin_cstr_normals_vec[i_ineq];

    foosteps_kin_cstr.block(cstr_index, step, ineq.rows(), 2) = ineq.block(0, 0, ineq.rows(), 2);
    b_kin_cstr.segment(cstr_index, ineq.rows()) = b_kin_cstr_vec[i_ineq].segment(0, ineq.rows());

    step += 2;
    cstr_index += static_cast<int>(ineq.rows());
  }

  eigen_assert(A_out.rows() == N_footsteps_kin_cstr);
  eigen_assert(A_out.cols() == N_variables);
  eigen_assert(b_out.rows() == A_out.rows());
  A_out.setZero();
  b_out.setZero();
  A_out.block(0, 0, N_footsteps_kin_cstr, Eigen::Index(2) * n_steps) = foosteps_kin_cstr * Delta;
  b_out.segment(0, N_footsteps_kin_cstr) = b_kin_cstr;
  Eigen::internal::set_is_malloc_allowed(true);
}

void feasibility_solver::build_steps_feasibility_matrix(Eigen::Ref<Eigen::MatrixXd> A_out,
                                                        Eigen::Ref<Eigen::VectorXd> b_out,
                                                        const sva::PTransformd & X_0_supp,
                                                        const sva::PTransformd & X_0_swg)
{
  Eigen::internal::set_is_malloc_allowed(false);
  const int N_slack = static_cast<int>(N_.rows());
  const int N_variables = 2 * N_steps + N_slack;

  const Eigen::Vector2d & P_supportFoot_0 = X_0_supp.translation().segment(0, 2);
  const Eigen::Vector2d & P_swingFoot_0 = X_0_swg.translation().segment(0, 2);

  // DCM must remain inside the feasibility region
  eigen_assert(A_out.rows() == N_.rows() && A_out.cols() == N_variables);
  eigen_assert(b_out.rows() == A_out.rows());
  A_out.setZero();
  b_out.setZero();
  // A_f * x + b = Of
  // A_f * x + b + slack >= N * P_u

  // We generate the cstr for each vertice of the rectangle

  double tds = optimalDoubleSupportDuration_[0];
  // i = 0

  for(int j = 0; j <= N_ds_; j++)
  {

    double alpha = static_cast<double>(j) / static_cast<double>(N_ds_);
    // if(!doubleSupport_){alpha = 1;}
    const double mu_0j = xTimings_[j];
    const double mu_0jp1 = xTimings_[j + 1];

    if(doubleSupport_ && j < N_ds_)
    {
      b_out += (offsetCstrZMPDblInit_ + N_ * (P_supportFoot_0 + P_swingFoot_0) * 0.5) * (mu_0j - mu_0jp1);
    }
    else { b_out += (offsetCstrZMP_ + N_ * P_supportFoot_0) * (mu_0j - mu_0jp1); }
  }

  // Remainings
  for(int i = 1; i <= N_steps; i++)
  {
    const int step_indx_im1 = 2 * (i - 1);
    const int step_indx_im2 = 2 * (i - 2);
    for(int j = 0; j <= (i != N_steps ? N_ds_ : N_tdsLast - 1); j++)
    {
      const double mu_ij = xTimings_[i * (N_ds_ + 1) + j];
      const double alpha = static_cast<double>(j + 1) / static_cast<double>(N_ds_ + 1);
      double mu_ijp1 = 0;
      if(!(i == N_steps && j == N_tdsLast - 1)) { mu_ijp1 = xTimings_[i * (N_ds_ + 1) + j + 1]; }

      A_out.block(0, step_indx_im1, 4, 2) += N_ * (mu_ij - mu_ijp1) * alpha;
      if(i > 1) { A_out.block(0, step_indx_im2, 4, 2) += N_ * (mu_ij - mu_ijp1) * (1 - alpha); }
      else { b_out += N_ * P_supportFoot_0 * (mu_ij - mu_ijp1) * (1 - alpha); }
      b_out += offsetCstrZMP_ * (mu_ij - mu_ijp1);
    }
  }
  b_out *= kappa_;
  b_out -= N_ * gamma_ * exp(-eta_ * t_);
  A_out *= kappa_;
  Eigen::internal::set_is_malloc_allowed(true);
}

bool feasibility_solver::solve_steps(const std::vector<sva::PTransformd> & refSteps)
{

  xTimings_ = Eigen::VectorXd::Zero(N_timings * (N_ds_ + 1) + N_tdsLast);
  double t_im1 = 0;
  for(int j = 0; j <= N_ds_; j++)
  {

    double alpha_j = static_cast<double>(j) / static_cast<double>(N_ds_);
    if(doubleSupport_) { xTimings_(j) = exp(-eta_ * (t_ + alpha_j * (optimalDoubleSupportDuration_[0] - t_))); }
    else { xTimings_(j) = exp(-eta_ * t_); }
  }
  for(int i = 1; i <= N_steps; i++)
  {

    t_im1 = i == 1 ? optimalStepsTimings_[0] : t_im1 + (optimalStepsTimings_[i - 1] - optimalStepsTimings_[i - 2]);

    for(int j = 0; j <= (i != N_steps ? N_ds_ : N_tdsLast - 1); j++)
    {

      double alpha_j = static_cast<double>(j) / static_cast<double>(N_ds_);
      xTimings_((N_ds_ + 1) * i + j) = exp(-eta_ * (t_im1 + alpha_j * (refTds_)));
    }
  }

  const int N_slack = static_cast<int>(N_.rows());
  size_t tstep_indx = 0;
  while(1.5 + t_ > optimalStepsTimings_[tstep_indx])
  {
    tstep_indx += 1;

    if(tstep_indx > optimalStepsTimings_.size())
    {
      tstep_indx -= 1;
      break;
    }
  }

  const int n_steps = tstep_indx;
  const int N_variables = 2 * n_steps + N_slack;
  if(N_variables == N_slack)
  {
    std::cout << "[Pendulum feasibility solver][Steps solver] "
              << "[iter : " << Niter_ << "] Next step too far in horizon" << std::endl;
    return true;
  }

  // DCM must remain inside the feasibility region
  Eigen::MatrixXd A_f = Eigen::MatrixXd::Zero(4, N_variables);
  Eigen::MatrixXd A_f_gen = Eigen::MatrixXd::Zero(N_.rows(), 2 * N_steps + N_.rows());
  Eigen::VectorXd b_f = Eigen::VectorXd::Zero(A_f_gen.rows());

  build_steps_feasibility_matrix(A_f_gen, b_f, X_0_SupportFoot_, X_0_SwingFoot_);
  if(n_steps < N_steps)
  {
    b_f += A_f_gen.block(0, 2 * n_steps, A_f_gen.rows(), 2 * (N_steps - n_steps))
           * xStep_.segment(2 * n_steps, 2 * (N_steps - n_steps));

    A_f.block(0, 0, A_f_gen.rows(), 2 * n_steps) = A_f_gen.block(0, 0, A_f_gen.rows(), 2 * n_steps);
  }

  // Kinematics Constraints
  Eigen::MatrixXd A_kin = Eigen::MatrixXd::Zero(4 * n_steps, N_variables);
  Eigen::VectorXd b_kin = Eigen::VectorXd::Zero(A_kin.rows());
  kinematics_constraints(A_kin, b_kin, refSteps_);

  if(A_ineq_buffer_.size() < (A_kin.rows() + A_f.rows()) * N_variables)
  {
    A_ineq_buffer_.resize((A_kin.rows() + A_f.rows()) * N_variables);
  }
  auto A_ineq_ = Eigen::Map<Eigen::MatrixXd>(A_ineq_buffer_.data(), A_kin.rows() + A_f.rows(), N_variables);
  if(b_ineq_buffer_.size() < A_ineq_.rows()) { b_ineq_buffer_.resize(A_ineq_.rows()); }
  auto b_ineq_ = Eigen::Map<Eigen::VectorXd>(b_ineq_buffer_.data(), A_ineq_.rows());
  A_ineq_ << -A_f * exp(eta_ * t_), A_kin;
  b_ineq_ << (b_f * exp(eta_ * t_) - (N_ * dcm_)), b_kin;

  // Slack Variables
  A_ineq_.block(0, 2 * n_steps, N_slack, N_slack) = Eigen::Matrix4d::Identity();

  const int NineqCstr = static_cast<int>(A_ineq_.rows());

  const int NeqCstr = 0;

  // Cost function

  Eigen::VectorXd x_init = Eigen::VectorXd::Zero(N_variables);
  // x_init.segment(0,2*n_steps) = b_steps;
  x_init.segment(0, 2 * n_steps) = xStep_.segment(0, 2 * n_steps);

  Eigen::Vector4d feasibilityOffsetInit = exp(eta_ * t_) * (A_f * x_init + b_f);
  Eigen::Vector4d dcm_pose = N_ * dcm_;

  // std::cout << "[Pendulum feasibility solver][Steps solver] init offset " << std::endl << feasibilityOffsetInit <<
  // std::endl;

  for(int i = 0; i < 4; i++)
  {
    if(feasibilityOffsetInit(i) < dcm_pose(i) - 1e5)
    {
      std::cout << "[Pendulum feasibility solver][Steps solver] "
                << "[iter : " << Niter_ << "] broken cstr on " << i << std::endl;
      std::cout << "offset delta " << feasibilityOffsetInit(i) - dcm_pose(i) << std::endl;
    }
  }

  if(Q_cost_buffer_.size() < N_variables * N_variables)
  {
    Q_cost_buffer_.resize(N_variables * N_variables);
    c_cost_buffer_.resize(N_variables);
  }
  auto Q_cost_ = Eigen::Map<Eigen::MatrixXd>(Q_cost_buffer_.data(), N_variables, N_variables);
  auto c_cost_ = Eigen::Map<Eigen::VectorXd>(c_cost_buffer_.data(), N_variables);
  Q_cost_.setZero();
  Q_cost_.topLeftCorner(2 * n_steps, 2 * n_steps).diagonal().setConstant(betaSteps);
  Q_cost_.bottomRightCorner(N_slack, N_slack).diagonal().setConstant(1e7);
  for(Eigen::Index i = 0; i < n_steps; i++)
  {
    c_cost_.segment(2 * i, 2) = -betaSteps * refSteps[i].translation().segment(0, 2);
  }
  c_cost_.tail(N_slack).setZero();

  qp_solver_.problem(N_variables, NeqCstr, NineqCstr);
  bool QPsuccess = qp_solver_.solve(Q_cost_, c_cost_, Eigen::MatrixXd{}, Eigen::VectorXd{}, A_ineq_, b_ineq_);

  // if(!QPsuccess)
  // {
  //     std::cout << "[Pendulum feasibility solver][Steps solver] " << "[iter : " << Niter_ <<"] QP Failed, lowering
  //     slack" << std::endl; M_slack.block(0,2 * n_steps,N_slack,N_slack) = 1e0 *
  //     Eigen::MatrixXd::Identity(N_slack,N_slack); Q_cost = betaSteps * ( M_steps.transpose() * M_steps) + (
  //     M_slack.transpose() * M_slack) ; c_cost = betaSteps * (-M_steps.transpose() * b_steps) + (-M_slack.transpose()
  //     * b_slack) ; QPsuccess = QP.solve(Q_cost, c_cost, A_eq, b_eq, A_ineq, b_ineq);
  // }

  if(!QPsuccess)
  {
    std::cout << "[Pendulum feasibility solver][Steps solver] "
              << "[iter : " << Niter_ << "] QP Failed" << std::endl;
    return false;
  }

  solution_ = qp_solver_.result();

  Eigen::Vector4d feasibilityOffset =
      exp(eta_ * t_) * (A_f.block(0, 0, 4, 2 * n_steps) * solution_.segment(0, 2 * n_steps) + b_f);
  // std::cout << "[Pendulum feasibility solver][Steps solver] output offset " << std::endl << feasibilityOffset <<
  // std::endl;

  for(int i = 0; i < 4; i++)
  {
    if(feasibilityOffset(i) < dcm_pose(i) - 1e-5)
    {
      std::cout << "[Pendulum feasibility solver][Steps solver] "
                << "[iter : " << Niter_ << "] solution broken cstr on " << i << std::endl;
      std::cout << "slack : " << solution_.segment(2 * n_steps + i, 1) << std::endl;
    }
  }
  Polygon feasibilityPolygon = Polygon(N_, feasibilityOffset);
  feasibilityRegion_ = feasibilityPolygon.Get_Polygone_Corners();

  // std::cout << "Slack Steps: " << solution_.segment(2 * n_steps,N_slack) << std::endl;

  // optimalSteps_.clear();
  for(int i = 0; i < n_steps; i++)
  {
    optimalSteps_[i] = (sva::PTransformd(refSteps[i].rotation(), Eigen::Vector3d{solution_(2 * i), solution_(2 * i + 1),
                                                                                 refSteps[i].translation().z()}));
    xStep_.segment(2 * i, 2) = solution_.segment(2 * i, 2);
  }
  return true;
}
