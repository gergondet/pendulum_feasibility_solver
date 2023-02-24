#include "../include/pendulum_feasibility_solver/feasibility_solver.h"
#include "../include/pendulum_feasibility_solver/polygons.h"

void feasibility_solver::timings_constraints(Eigen::MatrixXd & A_out, Eigen::VectorXd & b_out, const int NStepsTimings)
{
    const int NVariables = static_cast<int>(A_out.cols()); 
    const int N_mu = NVariables - static_cast<int>(N_.rows()); //N var minus the slack variables

    const double mu_ss_max = exp(-eta_ * t_ss_range_.y() );
    const double mu_ss_min = exp(-eta_ * t_ss_range_.x() );

    const double mu_ds_max = exp(-eta_ * t_ds_range_.y() );
    const double mu_ds_min = exp(-eta_ * t_ds_range_.x() );

    const double mu_s_max = exp(-eta_ * t_s_range_.y() );
    const double mu_s_min = exp(-eta_ * t_s_range_.x() );
    
    std::vector<Eigen::MatrixXd> A_vec;
    std::vector<Eigen::VectorXd> b_vec;
    
    int N_cstr = 0 ;
    for (int i  = 0 ; i < NStepsTimings ; i++)
    {
        Eigen::MatrixXd A_min = Eigen::MatrixXd::Zero(N_ds_+ 4,NVariables);
        //All mu must be decreasing in the horizon
        //Step must be bounded
        //Sg supp must be bounded
        //Ds supp must be bounded
        Eigen::VectorXd b_min = Eigen::VectorXd::Zero(A_min.rows());

        Eigen::MatrixXd A_max = Eigen::MatrixXd::Zero(3,NVariables);
        Eigen::VectorXd b_max = Eigen::VectorXd::Zero(A_max.rows());

        for (int j = 0 ; j <= N_ds_ ; j++)
        {
            //mu_i <= mu_i-1
            A_min.block( j , (N_ds_ + 1 ) * i + j,1,2) = Eigen::RowVector2d{-1 , 1};
        }
        
        //Step time cstr
        A_min(N_ds_ + 1, (i+1) * (N_ds_ + 1) ) = 1;

        //Ts <= ts_max
        A_max(0, (i+1) * (N_ds_ + 1)) = -1;

        //Ts - Tds <= tss_max
        A_max(1, (i+1) * (N_ds_ + 1)) = -1;
        A_max(1,  (N_ds_ + 1) * i + N_ds_ ) = mu_ss_max;

        
        

        if(i != 0)
        {
            //Step time cstr
            A_min(N_ds_ + 1,     i * (N_ds_ + 1) ) = -mu_s_min;
            
            //Sg supp time cstr
            A_min(N_ds_ , (N_ds_ + 1) * i + N_ds_) *= mu_ss_min;
            
            A_max(0, i * (N_ds_ + 1) ) = mu_s_max;

            // A_max(2,  (N_ds_ + 1) * i + N_ds_ ) = -1;
            // A_max(2,  (N_ds_ + 1) * i  ) = mu_ds_max;
        
        }
        
        else
        {
            if(!doubleSupport_)
            {
                //Sg supp time cstr
                //t_step >= t_ds + min_ss
                A_min(N_ds_ ,  N_ds_) *= 0;
                b_min(N_ds_) = mu_ss_min * exp(-eta_ * tLift_);

                // //Tstep >= t_ + > 0.1 
                // A_min(N_ds_ + 3 , N_ds_ + 1) = 1;
                // b_min(N_ds_ + 3) = exp(-eta_ * ( t_ + 0.1 ));

            }
            else
            {
                //Tds <= tds_max
                A_max(2,  (N_ds_ + 1) * i + N_ds_ ) = -1;
                b_max(2) = -mu_ds_max;
            }

            //t_step_0 >= t_s_min
            b_min(N_ds_ + 1) = mu_s_min;

            //mu_0 >=  t_
            A_min(N_ds_ + 2 , 0) = 1;
            b_min(N_ds_ + 2) = exp(-eta_ * t_);
            
            //Ts <= ts_max
            b_max(0) = -mu_s_max;
            
        }

        A_vec.push_back(A_min);
        b_vec.push_back(b_min);
        N_cstr += static_cast<int>(A_vec.back().rows());
        A_vec.push_back(A_max);
        b_vec.push_back(b_max);
        N_cstr += static_cast<int>(A_vec.back().rows());
    }

    A_out = Eigen::MatrixXd::Zero(N_cstr + 2 * N_mu,NVariables);
    b_out = Eigen::VectorXd::Zero(A_out.rows());
    
    Eigen::Index indx = 0;
    for (size_t i = 0 ; i < A_vec.size() ; i++)
    {
        const Eigen::MatrixXd & A = A_vec[i];
        const Eigen::VectorXd & b = b_vec[i];
        A_out.block(indx,0,A.rows(),NVariables) = A;

        b_out.segment(indx,b.rows()) = b;

        indx+= A.rows();
    }

    //mu must be positive and below 1
    A_out.block(N_cstr       ,0,N_mu,N_mu) = -Eigen::MatrixXd::Identity(N_mu,N_mu);
    A_out.block(N_cstr + N_mu,0,N_mu,N_mu) = Eigen::MatrixXd::Identity(N_mu,N_mu);
    b_out.segment(N_cstr + N_mu,N_mu) = Eigen::VectorXd::Ones(N_mu);

    
}

bool feasibility_solver::solve_timings(const std::vector<sva::PTransformd> & refSteps,const std::vector<double> & refTimings, const double & refTds)
{
    assert (refSteps.size() == refTimings.size());
    const int NStepsTimings = N_steps;
    const int N_slack = static_cast<int>(N_.rows());
    const int N_variables =  NStepsTimings * (N_ds_ + 1) + N_tdsLast + N_slack;

    

    const Eigen::Vector2d & P_supportFoot_0 = X_0_SupportFoot_.translation().segment(0,2);
    const Eigen::Vector2d & P_swingFoot_0 = X_0_SwingFoot_.translation().segment(0,2);

    const Eigen::Vector2d p_init = P_supportFoot_0 * t_/refTds_ + P_swingFoot_0 * (1 - t_/refTds_ );


    //Variables are organised as such  : timings then slack variables
    //step i occur at indx (N_ds + 1) * i
    //sg suport i start at indx (N_ds + 1) * i + N_ds


    //DCM must remain inside the feasibility region
    Eigen::MatrixXd A_f = Eigen::MatrixXd::Zero(N_.rows(),N_variables);
    Eigen::VectorXd b_f = Eigen::VectorXd::Zero(A_f.rows());
    //A_f * x + b >= N * P_u
    
    //i = 0 
    // const int j_start = doubleSupport_ ? 0 : N_ds_; 
    for (int j = 0 ; j <= N_ds_ ; j++)
    {
        const double alpha_j = static_cast<double>(j) / static_cast<double>(N_ds_);
        
        Eigen::Vector4d Oij = (offsetCstrZMP_ + N_ * (alpha_j * P_supportFoot_0 + (1 - alpha_j) * zmp_));
        if(!doubleSupport_)
        {
            Oij = offsetCstrZMP_ + N_ * alpha_j * P_supportFoot_0;
        }
        A_f.block(0, j , 4, 1)     += Oij;
        A_f.block(0, j + 1 , 4, 1) -= Oij;
    }
 
    //Remainings
    for (int i = 1 ; i <= N_steps; i++)
    {
        
        for (int j = 0 ; j <= (i != N_steps ? N_ds_ : N_tdsLast - 1 )  ; j++)
        {
            const double alpha_j = static_cast<double>(j) / static_cast<double>(N_ds_);
        
            Eigen::Vector4d O_ij = offsetCstrZMP_ + N_ * (alpha_j * refSteps[i-1].translation().segment(0,2));
            
            if( i > 1)
            {
                O_ij +=  N_ * ((1 - alpha_j) * refSteps[i-2].translation().segment(0,2));
            }
            else
            {   
                O_ij +=  N_ * ((1 - alpha_j) * P_supportFoot_0);
            }
            
            A_f.block(0,i * (N_ds_ + 1) + j     ,4,1) += O_ij;
            if(!(i == N_steps && j == N_tdsLast - 1))
            {
                A_f.block(0,i * (N_ds_ + 1) + j + 1 ,4,1) -= O_ij;
            }
            
            
        }

    }
    
    //Steps timings constraints
    Eigen::MatrixXd A_Tsteps = Eigen::MatrixXd::Zero(0,N_variables);
    Eigen::VectorXd b_Tsteps = Eigen::VectorXd::Zero(0);
    timings_constraints(A_Tsteps,b_Tsteps,NStepsTimings);

    //Slack variables cstr
    Eigen::MatrixXd Aineq_slack = Eigen::MatrixXd::Zero(2 * N_.rows(),N_variables);
    Aineq_slack.block(0,2 * N_steps,N_slack,N_slack) = Eigen::MatrixXd::Identity(N_slack,N_slack);
    Aineq_slack.block(N_slack,2 * N_steps,N_slack,N_slack) = -Eigen::MatrixXd::Identity(N_slack,N_slack);
    Eigen::VectorXd bineq_slack = Eigen::VectorXd::Ones(Aineq_slack.rows()) * 0.3;

    Eigen::MatrixXd A_ineq = Eigen::MatrixXd::Zero(A_Tsteps.rows() + A_f.rows() + Aineq_slack.rows(), N_variables);
    Eigen::VectorXd b_ineq = Eigen::VectorXd::Zero(A_ineq.rows());
    A_ineq <<       - A_f * exp(eta_ * t_)        , A_Tsteps, 0*Aineq_slack;
    b_ineq << b_f * exp(eta_ * t_) - (N_ * dcm_ ) , b_Tsteps, 0*bineq_slack;
    const int NineqCstr = static_cast<int>(A_ineq.rows()); 

    //Slack Variables
    A_ineq.block(0,N_variables - N_slack , N_slack , N_slack ) = Eigen::Matrix4d::Identity();


    Eigen::MatrixXd A_eq = Eigen::MatrixXd::Zero(0 , N_variables);
    Eigen::VectorXd b_eq = Eigen::VectorXd::Zero(A_eq.rows());
    const int NeqCstr = static_cast<int>(A_eq.rows()); 
    
    //Cost function
    Eigen::MatrixXd M_timings = Eigen::MatrixXd::Identity(N_variables - N_slack,N_variables);
    Eigen::VectorXd b_timings = Eigen::VectorXd::Ones(M_timings.rows());
    
    double t_im1 =  0;
    for (int j = 0 ; j <= N_ds_  ; j ++)
    {
        
        double alpha_j = static_cast<double>(j)/static_cast<double>(N_ds_);
        if(doubleSupport_)
        {
            b_timings(j) = exp(-eta_ * ( t_ + alpha_j * (refTds - t_)) );
        }
        else
        {
            b_timings(j) = exp(-eta_ * t_ );
        }

    }

    for(int i =  1; i <= NStepsTimings; i++)
    {
        t_im1 = refTimings[i-1];
        for (int j = 0 ; j <= (i != N_steps ? N_ds_ : N_tdsLast - 1 )  ; j ++)
        {
            
            double alpha_j = static_cast<double>(j)/static_cast<double>(N_ds_);
            b_timings( (N_ds_ + 1) * i + j) = exp(-eta_ * ( t_im1 + alpha_j * (refTds)) );

        }
        
    }

    Eigen::MatrixXd M_slack = Eigen::MatrixXd::Zero(N_slack,N_variables);
    M_slack.block(0,N_variables - N_slack,N_slack,N_slack) = 10 * Eigen::MatrixXd::Identity(N_slack,N_slack);
    Eigen::VectorXd b_slack = Eigen::VectorXd::Zero(M_slack.rows());

    //Keeping slack only on broken cstr
    // Eigen::VectorXd x_init = Eigen::VectorXd::Zero(N_variables);
    // x_init.segment(0,N_variables - N_slack) = b_timings;

    // Eigen::Vector4d feasibilityOffsetInit = exp(eta_ * t_) * ( A_f * x_init + b_f);
    // Eigen::Vector4d dcm_pose = N_ * dcm_;
    
    // for (int i = 0 ; i < 4 ; i++)
    // {
    //     if(feasibilityOffsetInit(i) < dcm_pose(i))
    //     {
    //         std::cout << "[Pendulum feasibility solver][Timing   solver] broken cstr on " << i << std::endl;
    //     }
    //     else
    //     {
    //         A_ineq.row(i).setZero();
    //         b_ineq(i) = 0;
    //     }
    // }

    
    Eigen::MatrixXd Q_cost = betaTsteps * M_timings.transpose() * M_timings + ( M_slack.transpose() * M_slack) ;
    Eigen::VectorXd c_cost = betaTsteps * (-M_timings.transpose() * b_timings) + (-M_slack.transpose() * b_slack) ;

    Eigen::QuadProgDense QP;
    // std::cout << A_ineq << std::endl;
    QP.problem(N_variables, NeqCstr,NineqCstr);
    bool QPsuccess = QP.solve(Q_cost, c_cost, A_eq, b_eq, A_ineq, b_ineq);



    if(!QPsuccess)
    {
        std::cout << "[Pendulum feasibility solver][Timing solver] QP Failed" << std::endl;
        Eigen::VectorXd Tds = Eigen::VectorXd::Ones(refTimings.size()) * (t_s_range_.x() - t_ss_range_.x());
        Eigen::VectorXd Ts = Eigen::VectorXd::Ones(refTimings.size()) * (t_s_range_.x());
        optimalDoubleSupportDuration_ = std::vector<double>(Tds.data(), Tds.data() + Tds.rows() );
        optimalStepsTimings_ = std::vector<double>(Ts.data(), Ts.data() + Ts.rows() );
        return false;
    }
    
    Eigen::VectorXd solution_ = QP.result();
    Eigen::Vector4d feasibilityOffset = exp(eta_ * t_) * ( A_f * solution_ + b_f);
    Polygon feasibilityPolygon = Polygon(N_,feasibilityOffset);
    feasibilityRegion_ = feasibilityPolygon.Get_Polygone_Corners();
    
    optimalStepsTimings_.clear();
    optimalDoubleSupportDuration_.clear();
    double tds_i = 0;
    for (int i = 0 ; i < NStepsTimings ; i++)
    {
        tds_i = -log(solution_( (N_ds_+1) * i + N_ds_))/eta_  - -log(solution_( (N_ds_ + 1 ) * i ))/eta_;

        if(i == 0){tds_i += t_;}
        
        optimalStepsTimings_.push_back(-log(solution_( (N_ds_ + 1) * (i+1)))/eta_);
        if(doubleSupport_)
        {
            optimalDoubleSupportDuration_.push_back(tds_i);
        }
        else
        {
            optimalDoubleSupportDuration_.push_back(refTds_);
        }

    } 
    return true;


}