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
    //All mu must be decreasing in the horizon
    Eigen::MatrixXd A_decrease = Eigen::MatrixXd::Zero(N_mu,NVariables);
    A_decrease.block(0,0,N_mu,N_mu) = Eigen::MatrixXd::Identity(N_mu,N_mu);
    A_decrease.block(1,0,N_mu - 1,N_mu - 1).diagonal() = -Eigen::VectorXd::Ones(N_mu - 1) ;
    A_decrease(1,0) *= exp(-eta_ * delta_);
    A_vec.push_back(A_decrease.block(1,0,N_mu - 1 , NVariables));
    b_vec.push_back(Eigen::VectorXd::Zero(A_vec.back().rows()));
    N_cstr += static_cast<int>(A_vec.back().rows());
    

    for (int i  = 0 ; i < NStepsTimings ; i++)
    {
        Eigen::MatrixXd A_min = Eigen::MatrixXd::Zero(3,NVariables);
        //Step must be bounded
        //Sg supp must be bounded
        Eigen::VectorXd b_min = Eigen::VectorXd::Zero(A_min.rows());

        Eigen::MatrixXd A_max = Eigen::MatrixXd::Zero(2,NVariables);
        Eigen::VectorXd b_max = Eigen::VectorXd::Zero(A_max.rows());

        A_min.block(0, (N_ds_ + 1) * i + N_ds_,1,2) = Eigen::RowVector2d{-1 , 1};

        //Sg supp time cstr
        A_min(0 , (N_ds_ + 1) * i + N_ds_) *= mu_ss_min;

        //Step time cstr
        A_min(1, (i+1) * (N_ds_ + 1) ) = 1;

        //Ds cstr
        A_min(2,(N_ds_ + 1) * i) = -mu_ds_min;
        A_min(2,(N_ds_ + 1) * i + N_ds_) = 1;

        //Ts <= ts_max
        A_max(0, (i+1) * (N_ds_ + 1)) = -1;

        //Ts - Tds <= tss_max
        A_max(1, (i+1) * (N_ds_ + 1)) = -1;
        A_max(1,  (N_ds_ + 1) * i + N_ds_ ) = mu_ss_max;

        
        

        if(i != 0)
        {
            //Step time cstr
            A_min( 1,     i * (N_ds_ + 1) ) = -mu_s_min;
            
            A_max(0, i * (N_ds_ + 1) ) = mu_s_max;

        
        }
        
        else
        {
            if(!doubleSupport_)
            {
                //Sg supp time cstr
                //t_step >= t_ds + min_ss
                A_min(0 ,  N_ds_) = 0;
                b_min(0) = mu_ss_min * exp(-eta_ * tLift_);

                //Ds cstr
                A_min.row(2).setZero();

            }
            else
            {
                A_min(2,0) = 0;
                b_min(2) = mu_ds_min;
            
            }
            

            //t_step_0 >= t_s_min
            b_min(1) = mu_s_min;
            
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

    A_out.resize(N_cstr + N_mu,NVariables);
    A_out.setZero();
    b_out.resize(A_out.rows());
    b_out.setZero();
    
    Eigen::Index indx = 0;
    for (size_t i = 0 ; i < A_vec.size() ; i++)
    {
        const Eigen::MatrixXd & A = A_vec[i];
        const Eigen::VectorXd & b = b_vec[i];
        A_out.block(indx,0,A.rows(),NVariables) = A;

        b_out.segment(indx,b.rows()) = b;

        indx+= A.rows();
    }

    //mu must be positive 
    A_out.block(N_cstr ,0,N_mu,N_mu) = -Eigen::MatrixXd::Identity(N_mu,N_mu);
    b_out.segment(N_cstr,N_mu) = -Eigen::VectorXd::Ones(N_mu) * exp(-eta_ * 10);
    
}

void feasibility_solver::build_time_feasibility_matrix(Eigen::MatrixXd & A_f, Eigen::VectorXd & b_f, const sva::PTransformd & X_0_supp, const sva::PTransformd & X_0_swg)
{
    const int NStepsTimings = N_steps;
    const int N_slack = static_cast<int>(N_.rows());
    const int N_variables =  NStepsTimings * (N_ds_ + 1) + N_tdsLast + N_slack;

    const Eigen::Vector2d & P_supportFoot_0 = X_0_supp.translation().segment(0,2);
    const Eigen::Vector2d & P_swingFoot_0 = X_0_swg.translation().segment(0,2);

    //Variables are organised as such  : timings then slack variables
    //step i occur at indx (N_ds + 1) * i
    //sg suport i start at indx (N_ds + 1) * i + N_ds


    //DCM must remain inside the feasibility region
    A_f.resize(N_.rows(),N_variables);
    b_f.resize(A_f.rows());
    A_f.setZero();
    b_f.setZero();
    //A_f * x + b + slck >= N * P_u
    
    //i = 0 
    // const int j_start = doubleSupport_ ? 0 : N_ds_; 
    for (int j = 0 ; j <= N_ds_ ; j++)
    {
        double alpha_j = static_cast<double>(j) / static_cast<double>(N_ds_);

        Eigen::Vector4d Oij = offsetCstrZMP_  + N_ * P_supportFoot_0 ;
        if(doubleSupport_ && j < N_ds_)
        {
            Oij = offsetCstrZMPDblInit_ + N_ * 0.5 * ( P_supportFoot_0 + P_swingFoot_0);
        }

        
        A_f.block(0, j , 4, 1)     += Oij;
        A_f.block(0, j + 1 , 4, 1) -= Oij;
    }
 
    //Remainings
    for (int i = 1 ; i <= NStepsTimings; i++)
    {
        
        for (int j = 0 ; j <= (i != NStepsTimings ? N_ds_ : N_tdsLast - 1 )  ; j++)
        {
            const double alpha_j = static_cast<double>(j + 1) / static_cast<double>(N_ds_ + 1);
        
            Eigen::Vector4d O_ij = offsetCstrZMP_ + N_ * (alpha_j * xStep_.segment(2 * (i - 1),2));
            
            if( i > 1)
            {
                O_ij +=  N_ * ((1 - alpha_j) * xStep_.segment(2 * (i - 2),2));
            }
            else
            {   
                O_ij +=  N_ * ((1 - alpha_j) * P_supportFoot_0);
            }
            
            A_f.block(0,i * (N_ds_ + 1) + j     ,4,1) += O_ij;
            if(!(i == NStepsTimings && j == N_tdsLast - 1))
            {
                A_f.block(0,i * (N_ds_ + 1) + j + 1 ,4,1) -= O_ij;
            }   
        }
    }
}

bool feasibility_solver::solve_timings(const std::vector<double> & refTimings, const double & refTds)
{
    const int NStepsTimings = N_steps;
    const int N_slack = static_cast<int>(N_.rows());
    const int N_variables =  NStepsTimings * (N_ds_ + 1) + N_tdsLast + N_slack;

    //Variables are organised as such  : timings then slack variables
    //step i occur at indx (N_ds + 1) * i
    //sg suport i start at indx (N_ds + 1) * i + N_ds


    //DCM must remain inside the feasibility region
    Eigen::MatrixXd A_f = Eigen::MatrixXd::Zero(N_.rows(),N_variables);
    Eigen::VectorXd b_f = Eigen::VectorXd::Zero(A_f.rows());
    //A_f * x + b + slck >= N * P_u
    
    build_time_feasibility_matrix(A_f,b_f,X_0_SupportFoot_,X_0_SwingFoot_);
        
    //Steps timings constraints
    Eigen::MatrixXd A_Tsteps = Eigen::MatrixXd::Zero(0,N_variables);
    Eigen::VectorXd b_Tsteps = Eigen::VectorXd::Zero(0);
    timings_constraints(A_Tsteps,b_Tsteps,NStepsTimings);

    Eigen::MatrixXd A_ineq = Eigen::MatrixXd::Zero(A_Tsteps.rows() + A_f.rows(), N_variables);
    Eigen::VectorXd b_ineq = Eigen::VectorXd::Zero(A_ineq.rows());
    A_ineq <<       - A_f * exp(eta_ * t_)        , A_Tsteps;
    b_ineq << b_f * exp(eta_ * t_) - (N_ * dcm_ ) , b_Tsteps;
    const int NineqCstr = static_cast<int>(A_ineq.rows()); 

    //Slack Variables
    A_ineq.block(0,N_variables - N_slack , N_slack , N_slack ) = Eigen::Matrix4d::Identity();

    // std::cout << A_f << std::endl;

    Eigen::MatrixXd A_eq = Eigen::MatrixXd::Zero(1 , N_variables);
    A_eq(0,0) = 1;
    Eigen::VectorXd b_eq = Eigen::VectorXd::Zero(A_eq.rows());
    b_eq(0) = exp(-eta_ * t_);
    const int NeqCstr = static_cast<int>(A_eq.rows()); 
    
    //Cost function
    Eigen::MatrixXd M_deltaTs = Eigen::MatrixXd::Zero(2 * N_timings,N_variables);
    Eigen::MatrixXd b_deltaTs = Eigen::VectorXd::Zero(M_deltaTs.rows());

    for(int i = 1 ; i < N_timings ; i ++)
    {
        M_deltaTs.block(2 * i,(N_ds_ + 1) * i + N_ds_,1,2) = Eigen::RowVector2d{-1,1};

        if(i != 0)
        {
            M_deltaTs(2 * i     , (N_ds_ + 1) * i + N_ds_) *= exp(-eta_ * ( refTimings[i] - refTimings[i-1] - refTds_));
            M_deltaTs(2 * i + 1 , (N_ds_ + 1) * i + N_ds_) = 1;
            M_deltaTs(2 * i + 1 , (N_ds_ + 1) * i) = -exp(-eta_ * (refTds_));
        }
        else
        {
            M_deltaTs(2 * i     , N_ds_ ) *= exp(-eta_ * ( refTimings[0] - refTds_));
            if(!doubleSupport_)
            {
                M_deltaTs(2 * i     , (N_ds_ + 1) * i + N_ds_) = 0;
                b_deltaTs(2 * i) = exp(-eta_ * (tLift_ + refTimings[0] - refTds_) );
            }
        }

    }

    Eigen::MatrixXd M_timings = Eigen::MatrixXd::Identity(N_variables - N_slack,N_variables);
    Eigen::VectorXd b_timings = Eigen::VectorXd::Ones(M_timings.rows());

    Eigen::MatrixXd M_slack = Eigen::MatrixXd::Zero(N_slack,N_variables);
    M_slack.block(0,N_variables - N_slack,N_slack,N_slack) = Eigen::MatrixXd::Identity(N_slack,N_slack);
    
    // std::cout << "n timings " << N_timings << std::endl;
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


    //Keeping slack only on broken cstr
    Eigen::VectorXd x_init = Eigen::VectorXd::Zero(N_variables);
    x_init.segment(0,N_variables - N_slack) = b_timings;

    Eigen::Vector4d feasibilityOffsetInit = exp(eta_ * t_) * ( A_f * x_init + b_f);
    Eigen::Vector4d dcm_pose  = N_ * dcm_;

    // std::cout << "[Pendulum feasibility solver][Timing solver] xStep " << std::endl << xStep_ << std::endl;
    // std::cout << "[Pendulum feasibility solver][Timing solver] init offset " << std::endl << feasibilityOffsetInit << std::endl;
    // bool ok = true;
    for (int i = 0 ; i < 4 ; i++)
    {
        if(feasibilityOffsetInit(i) < dcm_pose(i))
        {
            std::cout << "[Pendulum feasibility solver][Timing  solver] " << "[iter : " << Niter_ <<"] broken cstr on " << i << std::endl;
            // ok = false;
        }
    }
    // if(ok && Niter_ == 0)
    // {
    //     return true;
    // }

    
    Eigen::MatrixXd Q_cost = betaTsteps * M_timings.transpose() * M_timings;
    Q_cost += 5e0 * M_deltaTs.transpose() * M_deltaTs;
    Q_cost += 1e5 * M_slack.transpose() * M_slack;
    Eigen::VectorXd c_cost = betaTsteps * (-M_timings.transpose() * b_timings) ;
    c_cost += 5e0 * -M_deltaTs.transpose() * b_deltaTs;

    for (int i = 0 ; i < N_timings ; i++)
    {

        Eigen::MatrixXd M_plan = Eigen::MatrixXd::Zero(2,N_variables);
        Eigen::VectorXd b_plan = Eigen::VectorXd::Zero(M_plan.rows());
        M_plan(0,(i) * (N_ds_ + 1)) = 1;
        M_plan(0,(i + 1) * (N_ds_ + 1)) = -1;
        M_plan(1,(N_ds_ + 1) * i + N_ds_) = -10;
        M_plan(1,(N_ds_ + 1) * i) = 10;
        const double steps_error = (optimalSteps_[i].translation() - refSteps_[i].translation()).norm();  
        Q_cost += 5e0 * steps_error * M_plan .transpose() * M_plan;
        c_cost += 5e0 * steps_error * -M_plan.transpose() * b_plan; 
        
    }


    Eigen::QuadProgDense QP;
    // std::cout << A_ineq << std::endl;
    // std::cout << b_ineq << std::endl;
    QP.problem(N_variables, NeqCstr,NineqCstr);
    bool QPsuccess = QP.solve(Q_cost, c_cost, A_eq, b_eq, A_ineq, b_ineq);

    // std::cout << "init " << Niter_ << std::endl;
    // for (int i = 0 ; i < N_variables - N_slack ; i++)
    // {
    //     std::cout << -log(b_timings(i))/eta_ << std::endl;
    // }
    // std::cout << "dcm" << std::endl << dcm_ << std::endl;


    if(!QPsuccess)
    {
        std::cout << "[Pendulum feasibility solver][Timing solver] " << "[iter : " << Niter_ <<"] QP Failed" << std::endl;
        return true;
    }
    
    xTimings_ = QP.result().segment(0,N_variables - N_slack);
    Eigen::Vector4d feasibilityOffset = exp(eta_ * t_) * ( A_f * QP.result() + b_f);
    // std::cout << "[Pendulum feasibility solver][Timing solver] output offset " << std::endl << feasibilityOffset << std::endl;
    Polygon feasibilityPolygon = Polygon(N_,feasibilityOffset);
    feasibilityRegion_ = feasibilityPolygon.Get_Polygone_Corners();

    // std::cout << "sol " << Niter_ << std::endl;
    // for (int i = 0 ; i < N_variables - N_slack ; i++)
    // {
    //     std::cout << -log(xTimings_(i))/eta_ << std::endl;
    // }
    // std::cout << "Slack Timings: " << QP.result().segment(N_variables - N_slack ,N_slack) << std::endl;


    
    optimalStepsTimings_.clear();
    optimalDoubleSupportDuration_.clear();
    double tds_i = -log(xTimings_( N_ds_))/eta_ ;
    double ts_i = -log(xTimings_( (N_ds_ + 1)))/eta_;
    optimalStepsTimings_.push_back(ts_i);
    optimalDoubleSupportDuration_.push_back(tds_i);

    for (int i = 1 ; i < NStepsTimings ; i++)
    {
        tds_i = (-log(xTimings_( (N_ds_+1) * i + N_ds_))/eta_)  - ts_i;
        ts_i = -log(xTimings_( (N_ds_ + 1) * (i+1)))/eta_;

        
        optimalStepsTimings_.push_back(ts_i);
        optimalDoubleSupportDuration_.push_back(tds_i);

    } 
    if(!doubleSupport_)
    {
        if(optimalDoubleSupportDuration_.size() >= 2)
        {
            optimalDoubleSupportDuration_[0] = optimalDoubleSupportDuration_[1];
        }
        else
        {
            optimalDoubleSupportDuration_[0] = refTds;
        }
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