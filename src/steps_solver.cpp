#include "../include/pendulum_feasibility_solver/feasibility_solver.h"
#include "../include/pendulum_feasibility_solver/polygons.h"

void feasibility_solver::kinematics_contraints(Eigen::MatrixXd & A_out, Eigen::VectorXd & b_out)
{
    const int N_variables = static_cast<int>(A_out.cols());
    const int N_steps = static_cast<int>(refSteps_.size());
    assert(N_steps > 0);
    std::vector<Eigen::VectorXd> b_kin_cstr_vec;
    std::vector<Eigen::MatrixX2d> kin_cstr_normals_vec;
    std::vector<Eigen::MatrixX2d> step_cstr_normals_vec;
    std::vector<Eigen::VectorXd> b_step_cstr_vec;
    Eigen::MatrixXd Delta = Eigen::MatrixXd::Identity(2 * N_steps, 2 * N_steps); // Matrix to differentiate two footsteps

    double l = feetDistance_;
    if(supportFoot_ == "LeftFoot"){l*=-1;}
    int N_footsteps_kin_cstr = 0;
    for(int i = 0; i < N_steps; i++)
    {
        const double theta_i = rpyFromMat(refSteps_[i].rotation()).z();
        sva::PTransformd X_0_step_im1 = X_0_SupportFoot_;
        if(i != 0)
        {
            X_0_step_im1 = refSteps_[i - 1];
        }
        const Eigen::Matrix3d R_Theta_i_0 = X_0_step_im1.rotation().transpose();
        const Eigen::Vector3d offset = R_Theta_i_0 * Eigen::Vector3d{0,l + (l/std::abs(l)) * stepCstrSize_.y()/2,0};

        Rectangle Kinematic_Rectangle = Rectangle(theta_i, stepCstrSize_, offset);

        if(i > 0)
        {
            Delta.block(2*i,2*(i-1),2,2) = -Eigen::Matrix2d::Identity();
        }
        if(i == 0)
        {
            Kinematic_Rectangle = Rectangle(X_0_SupportFoot_, stepCstrSize_, offset);
        }
        Polygon Kinematic_Poly = Polygon(Kinematic_Rectangle);
        b_kin_cstr_vec.push_back(Kinematic_Poly.offsets());
        kin_cstr_normals_vec.push_back(Kinematic_Poly.normals());

        N_footsteps_kin_cstr += static_cast<int>(kin_cstr_normals_vec.back().rows());
        l*=-1;
    }

    Eigen::MatrixXd foosteps_kin_cstr = Eigen::MatrixXd::Zero(N_footsteps_kin_cstr, 2 * N_steps);
    Eigen::VectorXd b_kin_cstr(N_footsteps_kin_cstr);
    A_out.resize(N_footsteps_kin_cstr , N_variables);
    A_out.setZero();
    b_out.resize(N_footsteps_kin_cstr);
    b_out.setZero();

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

    A_out.block(0, 0, N_footsteps_kin_cstr, 2 * N_steps) = foosteps_kin_cstr * Delta;
    b_out.segment(0, N_footsteps_kin_cstr) = b_kin_cstr;
}

bool feasibility_solver::solve_steps()
{
    assert (refSteps_.size() == refTimings_.size());
    const int N_tds = static_cast<int>(refTds_ / delta_);
    const int N_tdsLast = static_cast<int>(refTds_ / (2*delta_));
    const int N_steps =  static_cast<int>(refSteps_.size());

    const int N_variables = 2 * N_steps;

    const Eigen::Vector2d & P_supportFoot_0 = X_0_SupportFoot_.translation().segment(0,2);
    const Eigen::Vector2d & P_swingFoot_0 = X_0_SwingFoot_.translation().segment(0,2);
    
    //cst part of the offsets
    Eigen::Matrix<double,4,1> offset_cstr_zmp; 
    offset_cstr_zmp << zmpRange_.x()/2,  zmpRange_.x()/2, zmpRange_.y()/2,  zmpRange_.y()/2;

    //DCM must remain inside the feasibility region
    Eigen::MatrixXd A_f = Eigen::MatrixXd::Zero(4,N_variables);
    Eigen::VectorXd b_f = Eigen::VectorXd::Zero(A_f.rows());
    //A_f * x + b >= N * P_u

    //We generate the cstr for each vertice of the rectangle
    
    //i = 0 
    const int j_start = doubleSupport_ ? 0 : N_tds;
    for (int j = j_start ; j < N_tds ; j++)
    {
        const double alpha_j = static_cast<double>(j) / static_cast<double>(N_tds);
        const double e_alpha_j_Tds = exp(-eta_ * alpha_j * refDoubleSupportTimings_[0]); 
        const double alpha_jp1 = static_cast<double>(j+1) / static_cast<double>(N_tds);
        const double e_alpha_jp1_Tds = exp(-eta_ * alpha_jp1 * refDoubleSupportTimings_[0]); 
        b_f += (offset_cstr_zmp
                + N_ * (alpha_j * P_supportFoot_0 + (1 - alpha_j) * zmp_)) * 
                (e_alpha_j_Tds - e_alpha_jp1_Tds);
    }
    b_f += (offset_cstr_zmp + N_ * P_supportFoot_0 ) * (exp(-eta_ * refTds_) - exp(-eta_ * refTimings_[0]));

    //Remainings
    for (int i = 1 ; i <= N_steps; i++)
    {
        const int step_indx_im1 = 2 * (i-1);
        const int step_indx_im2 = 2 * (i-2);
        const double mu_im1 = exp(-eta_ * refTimings_[i - 1]);
        for (int j = 0 ; j < (i != N_steps ? N_tds : N_tdsLast )  ; j++)
        {
            const double alpha_j = static_cast<double>(j) / static_cast<double>(N_tds);
            const double e_alpha_j_Tds = exp(-eta_ * alpha_j * refDoubleSupportTimings_[i]); 
            const double alpha_jp1 = static_cast<double>(j+1) / static_cast<double>(N_tds);
            double e_alpha_jp1_Tds = exp(-eta_ * alpha_jp1 * refDoubleSupportTimings_[i]); 
            if( i == N_steps && j == N_tdsLast - 1)
            {
                e_alpha_jp1_Tds = 0;
            }

            
            A_f.block(0,step_indx_im1,4,2) += N_ * (e_alpha_j_Tds - e_alpha_jp1_Tds) * alpha_j * mu_im1;
            if(i > 1)
            {
                A_f.block(0,step_indx_im2,4,2) += N_ * (e_alpha_j_Tds - e_alpha_jp1_Tds) * (1  - alpha_j) * mu_im1;
            }
            else
            {
                b_f += N_ * P_supportFoot_0 * (e_alpha_j_Tds - e_alpha_jp1_Tds) * (1  - alpha_j) * mu_im1;
            }
            b_f += offset_cstr_zmp * (e_alpha_j_Tds - e_alpha_jp1_Tds) * mu_im1;
            
        }
        
        if(i != N_steps)
        {
            const double mu_i =  exp(-eta_ * refTimings_[i]);
            A_f.block(0,step_indx_im1,4,2) += N_ * ( mu_im1 * exp(-eta_ * refTds_) - mu_i);
            b_f += offset_cstr_zmp * ( mu_im1 * exp(-eta_ * refTds_) - mu_i);
            
        }


    }
    
    //Kinematics Constraints
    Eigen::MatrixXd A_kin = Eigen::MatrixXd::Zero(0,N_variables);
    Eigen::VectorXd b_kin = Eigen::VectorXd::Zero(0);
    kinematics_contraints(A_kin,b_kin);

    Eigen::MatrixXd A_ineq = Eigen::MatrixXd::Zero(A_kin.rows() + A_f.rows() , N_variables);
    Eigen::VectorXd b_ineq = Eigen::VectorXd::Zero(A_ineq.rows());
    A_ineq <<       - A_f         , A_kin;
    b_ineq << (b_f - (N_ * dcm_)) , b_kin;
    const int NineqCstr = static_cast<int>(A_ineq.rows()); 

    Eigen::MatrixXd A_eq = Eigen::MatrixXd::Zero(0 , N_variables);
    Eigen::VectorXd b_eq = Eigen::VectorXd::Zero(A_eq.rows());
    const int NeqCstr = static_cast<int>(A_eq.rows()); 
    
    //Cost function
    
    Eigen::MatrixXd M_steps = Eigen::MatrixXd::Identity(N_variables,N_variables);
    Eigen::VectorXd b_steps = Eigen::VectorXd::Zero(N_variables);
    for(int i = 0; i < N_steps; i++)
    {
        b_steps.segment(2 * i, 2) = refSteps_[i].translation().segment(0, 2);
    }
    

    const Eigen::MatrixXd Q_cost = betaSteps * M_steps.transpose() * M_steps 
                            + 1e-12 * Eigen::MatrixXd::Identity(N_variables,N_variables);
    const Eigen::VectorXd c_cost = betaSteps * (-M_steps.transpose() * b_steps) ;

    Eigen::QuadProgDense QP;
    QP.problem(N_variables, NeqCstr,NineqCstr);
    bool QPsuccess = QP.solve(Q_cost, c_cost, A_eq, b_eq, A_ineq, b_ineq);



    if(!QPsuccess)
    {
        return false;
    }
    
    Eigen::VectorXd solution_ = QP.result();
    // std::cout << solution_ << std::endl;
    

    optimalSteps_.clear();

    for (int i = 0 ; i < N_steps ; i++)
    {
        optimalSteps_.push_back(sva::PTransformd(
                                    refSteps_[i].rotation(),
                                    Eigen::Vector3d{solution_(2 * i),solution_(2 * i + 1),refSteps_[i].translation().z()})
                                );

    } 
    return true;


}