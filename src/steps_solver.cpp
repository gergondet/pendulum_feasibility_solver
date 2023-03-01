#include "../include/pendulum_feasibility_solver/feasibility_solver.h"
#include "../include/pendulum_feasibility_solver/polygons.h"

void feasibility_solver::kinematics_contraints(Eigen::MatrixXd & A_out, Eigen::VectorXd & b_out, const std::vector<sva::PTransformd> & refSteps)
{
    const int N_variables = static_cast<int>(A_out.cols());
    const int N_steps = static_cast<int>(refSteps.size());
    // const int N_steps = 1;

    assert(N_steps > 0);
    std::vector<Eigen::VectorXd> b_kin_cstr_vec;
    std::vector<Eigen::MatrixX2d> kin_cstr_normals_vec;
    std::vector<Eigen::MatrixX2d> step_cstr_normals_vec;
    std::vector<Eigen::VectorXd> b_step_cstr_vec;
    Eigen::MatrixXd Delta = Eigen::MatrixXd::Identity(2 * N_steps, 2 * N_steps); // Matrix to differentiate two footsteps

    double l = 1;
    if(supportFoot_ == "LeftFoot"){l*=-1;}
    int N_footsteps_kin_cstr = 0;
    for(int i = 0; i < N_steps; i++)
    {
        const double theta_i = rpyFromMat(refSteps[i].rotation()).z();
        sva::PTransformd X_0_step_im1 = X_0_SupportFoot_;
        if(i != 0)
        {
            X_0_step_im1 = refSteps[i - 1];
        }
        const Eigen::Matrix3d R_Theta_i_0 = X_0_step_im1.rotation().transpose();
        const Eigen::Vector3d offset = R_Theta_i_0 * Eigen::Vector3d{0,l*(feetDistance_ + stepCstrSize_.y()/2),0};

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

    A_out = Eigen::MatrixXd::Zero(N_footsteps_kin_cstr , N_variables);
    b_out = Eigen::VectorXd::Zero(N_footsteps_kin_cstr);
    A_out.block(0, 0, N_footsteps_kin_cstr, 2 * N_steps) = foosteps_kin_cstr * Delta;
    b_out.segment(0, N_footsteps_kin_cstr) = b_kin_cstr;


}

bool feasibility_solver::solve_steps(const std::vector<sva::PTransformd> & refSteps)
{

    const int N_slack = static_cast<int>(N_.rows());
    const int N_variables = 2 * N_steps + N_slack;

    const Eigen::Vector2d & P_supportFoot_0 = X_0_SupportFoot_.translation().segment(0,2);
    const Eigen::Vector2d & P_swingFoot_0 = X_0_SwingFoot_.translation().segment(0,2);
    

    //DCM must remain inside the feasibility region
    Eigen::MatrixXd A_f = Eigen::MatrixXd::Zero(4,N_variables);
    Eigen::VectorXd b_f = Eigen::VectorXd::Zero(A_f.rows());
    // A_f * x + b = Of + slack
    //A_f * x + b >= N * P_u

    //We generate the cstr for each vertice of the rectangle
    
    double tds = optimalDoubleSupportDuration_[0];
    const Eigen::Vector2d p_init = P_supportFoot_0 * t_/tds + P_swingFoot_0 * (1 - t_/tds );
    //i = 0 
 
        
    for (int j = 0 ; j <= N_ds_ ; j++)
    {

        double alpha = static_cast<double>(j) / static_cast<double>(N_ds_);
        if(!doubleSupport_){alpha = 1;}
        const double mu_0j = xTimings_[j] ; 
        const double mu_0jp1 = xTimings_[j+1] ; 
    
        b_f += (offsetCstrZMP_ 
                + N_ * (alpha * P_supportFoot_0 + (1 - alpha) * p_init)) * 
                (mu_0j - mu_0jp1);
    }

    //Remainings
    for (int i = 1 ; i <= N_steps; i++)
    {
        const int step_indx_im1 = 2 * (i-1);
        const int step_indx_im2 = 2 * (i-2);
        for (int j = 0 ; j <= (i != N_steps ? N_ds_ : N_tdsLast - 1)  ; j++)
        {
            const double mu_ij =  xTimings_[ i * (N_ds_ + 1) + j ]; 
            const double alpha = static_cast<double>(j + 1) / static_cast<double>(N_ds_ + 1);
            double mu_ijp1 = 0;
            if( !(i == N_steps && j == N_tdsLast - 1) )
            {
                mu_ijp1 = xTimings_[ i * (N_ds_ + 1) + j + 1 ]; 
            }

            
            A_f.block(0,step_indx_im1,4,2) += N_ * (mu_ij - mu_ijp1) * alpha;
            if(i > 1)
            {
                A_f.block(0,step_indx_im2,4,2) += N_ * (mu_ij - mu_ijp1) * (1 - alpha) ;
            }
            else
            {
                b_f += N_ * P_supportFoot_0 * (mu_ij - mu_ijp1) * (1 - alpha);
            }
            b_f += offsetCstrZMP_ * (mu_ij - mu_ijp1) ;
            
        }
        
    }


    //Kinematics Constraints
    Eigen::MatrixXd A_kin = Eigen::MatrixXd::Zero(0,N_variables);
    Eigen::VectorXd b_kin = Eigen::VectorXd::Zero(0);
    kinematics_contraints(A_kin,b_kin,refSteps_);

    Eigen::MatrixXd A_ineq = Eigen::MatrixXd::Zero(A_kin.rows() + A_f.rows() , N_variables);
    Eigen::VectorXd b_ineq = Eigen::VectorXd::Zero(A_ineq.rows());
    A_ineq <<           - A_f *  exp(eta_ * t_)        , A_kin ;
    b_ineq <<    (b_f * exp(eta_ * t_) - (N_ * dcm_) ) , b_kin ;

    // Slack Variables
    A_ineq.block(0,2 * N_steps , N_slack , N_slack ) = Eigen::Matrix4d::Identity();

    const int NineqCstr = static_cast<int>(A_ineq.rows()); 

    Eigen::MatrixXd A_eq = Eigen::MatrixXd::Zero(0 , N_variables);
    Eigen::VectorXd b_eq = Eigen::VectorXd::Zero(A_eq.rows());
    const int NeqCstr = static_cast<int>(A_eq.rows()); 
    
    //Cost function
    
    Eigen::MatrixXd M_steps = Eigen::MatrixXd::Zero(2 * N_steps,N_variables);
    M_steps.block(0,0,2*N_steps,2*N_steps) = Eigen::MatrixXd::Identity(2 * N_steps,2 * N_steps);
    // M_steps.block(2,2,2 *(N_steps - 1) , 2 * (N_steps - 1)).diagonal() *= 10;
    Eigen::VectorXd b_steps = Eigen::VectorXd::Zero(M_steps.rows());
    for(int i = 0; i < N_steps; i++)
    {
        b_steps.segment(2 * i, 2) = refSteps[i].translation().segment(0, 2);
    }
    
    Eigen::MatrixXd M_slack = Eigen::MatrixXd::Zero(N_slack,N_variables);
    M_slack.block(0,2 * N_steps,N_slack,N_slack) = 1e5 * Eigen::MatrixXd::Identity(N_slack,N_slack);
    Eigen::VectorXd b_slack = Eigen::VectorXd::Zero(M_slack.rows());

    //Keeping slack only on broken cstr
    Eigen::VectorXd x_init = Eigen::VectorXd::Zero(N_variables);
    x_init.segment(0,2*N_steps) = b_steps;

    Eigen::Vector4d feasibilityOffsetInit = exp(eta_ * t_) * ( A_f * x_init + b_f);
    Eigen::Vector4d dcm_pose = N_ * dcm_;

    // std::cout << "[Pendulum feasibility solver][Steps solver] init offset " << std::endl << feasibilityOffsetInit << std::endl;

    
    for (int i = 0 ; i < 4 ; i++)
    {
        if(feasibilityOffsetInit(i) < dcm_pose(i))
        {
            std::cout << "[Pendulum feasibility solver][Steps solver] " << "[iter : " << Niter_ <<"] broken cstr on " << i << std::endl;
        }
    }
    // std::cout << "A_ineq" << std::endl << A_ineq << std::endl;
    // std::cout << "b_ineq" << std::endl << b_ineq << std::endl;
    // std::cout << "N(i) " << std::endl << N_.row(i) << std::endl;

    const Eigen::MatrixXd Q_cost = betaSteps * ( M_steps.transpose() * M_steps) + ( M_slack.transpose() * M_slack) ;
    const Eigen::VectorXd c_cost = betaSteps * (-M_steps.transpose() * b_steps) + (-M_slack.transpose() * b_slack) ;

    Eigen::QuadProgDense QP;
    // QP.tolerance(5e-4);
    QP.problem(N_variables, NeqCstr,NineqCstr);
    bool QPsuccess = QP.solve(Q_cost, c_cost, A_eq, b_eq, A_ineq, b_ineq);



    if(!QPsuccess)
    {
        std::cout << "[Pendulum feasibility solver][Steps solver] " << "[iter : " << Niter_ <<"] QP Failed" << std::endl;
        return false;
    }
    
    solution_ = QP.result();

    Eigen::Vector4d feasibilityOffset = exp(eta_ * t_) * ( A_f * solution_ + b_f);
    // std::cout << "[Pendulum feasibility solver][Steps solver] output offset " << std::endl << feasibilityOffset << std::endl;


    for (int i = 0 ; i < 4 ; i++)
    {
        if(feasibilityOffset(i) < dcm_pose(i))
        {
            std::cout << "[Pendulum feasibility solver][Steps solver] " << "[iter : " << Niter_ <<"] solution broken cstr on " << i << std::endl;
        }
    }
    // Polygon feasibilityPolygon = Polygon(N_,feasibilityOffset);
    // feasibilityRegion_ = feasibilityPolygon.Get_Polygone_Corners();

    // std::cout << "Slack Steps: " << solution_.segment(2 * N_steps,N_slack) << std::endl;
    

    optimalSteps_.clear();
    for (int i = 0 ; i < N_steps ; i++)
    {
        optimalSteps_.push_back(sva::PTransformd(
                                    refSteps[i].rotation(),
                                    Eigen::Vector3d{solution_(2 * i),solution_(2 * i + 1),refSteps[i].translation().z()})
                                );
        xStep_.segment(2 * i , 2) = solution_.segment(2 * i,2);


    } 
    return true;


}