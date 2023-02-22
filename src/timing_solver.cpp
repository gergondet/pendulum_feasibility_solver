#include "../include/pendulum_feasibility_solver/feasibility_solver.h"
#include "../include/pendulum_feasibility_solver/polygons.h"

void feasibility_solver::timings_constraints(Eigen::MatrixXd & A_out, Eigen::VectorXd & b_out)
{
    const int NVariables = static_cast<int>(A_out.cols());
    const int NStepsTimings = static_cast<int>(refTimings_.size());

    const double mu_ss_max = exp(-eta_ * t_ss_range_.y() );
    const double mu_ss_min = exp(-eta_ * t_ss_range_.x() );

    const double mu_s_max = exp(-eta_ * t_s_range_.y() );
    const double mu_s_min = exp(-eta_ * t_s_range_.x() );
    
    std::vector<Eigen::MatrixXd> A_vec;
    std::vector<Eigen::VectorXd> b_vec;
    
    int N_cstr = 0 ;
    for (int i  = 0 ; i < NStepsTimings ; i++)
    {
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N_ds_+ 3,NVariables);
        Eigen::VectorXd b = Eigen::VectorXd::Zero(A.rows());

        for (int j = 0 ; j <= N_ds_ ; j++)
        {

            A.block( j , (N_ds_ + 1 ) * i + j,1,2) = Eigen::RowVector2d{-1 , 1};
        }
        

        
        A(N_ds_ + 1, (i+1) * (N_ds_ + 1) ) = 1;
        if(i != 0)
        {
            A(N_ds_ + 1,     i * (N_ds_ + 1) ) = -mu_s_min;
            A(N_ds_ , (N_ds_ + 1) * i + N_ds_) *= mu_ss_min;
        }
        
        else
        {
            if(!doubleSupport_)
            {
                A(N_ds_ ,  N_ds_) = 0;
                b(N_ds_) = mu_ss_min * exp(-eta_ * (-tLift_));
            }
            b(N_ds_ + 1) = mu_s_min;
            A(N_ds_ + 2 , 0) = 1;
            b(N_ds_ + 2) = exp(-eta_ * t_);
        }

        A_vec.push_back(A);
        b_vec.push_back(b);
        N_cstr += static_cast<int>(A_vec.back().rows());
    }

    A_out = Eigen::MatrixXd::Zero(N_cstr + NVariables,NVariables);
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

    //mu must be positive
    A_out.block(N_cstr,0,NVariables,NVariables) = -Eigen::MatrixXd::Identity(NVariables,NVariables);
    
}

bool feasibility_solver::solve_timings()
{
    assert (refSteps_.size() == refTimings_.size());
    const int N_tdsLast = static_cast<int>(static_cast<double>(N_ds_) / 2);
    const int N_steps =  static_cast<int>(refSteps_.size());
    const int NStepsTimings = N_steps;

    const int N_variables =  NStepsTimings * (N_ds_ + 1) + N_tdsLast;

    const Eigen::Vector2d & P_supportFoot_0 = X_0_SupportFoot_.translation().segment(0,2);
    const Eigen::Vector2d & P_swingFoot_0 = X_0_SwingFoot_.translation().segment(0,2);


    //Variables are organised as such  : timings then footsteps
    //footstep i indx : NStepsTimings + 2 * i 



    
    //cst part of the offsets
    Eigen::Matrix<double,4,1> offset_cstr_zmp; 
    offset_cstr_zmp << zmpRange_.x()/2,  zmpRange_.x()/2, zmpRange_.y()/2,  zmpRange_.y()/2;

    //DCM must remain inside the feasibility region
    Eigen::MatrixXd A_f = Eigen::MatrixXd::Zero(N_.rows(),N_variables);
    Eigen::VectorXd b_f = Eigen::VectorXd::Zero(A_f.rows());
    //A_f * x + b >= N * P_u
    
    //i = 0 
    const int j_start = doubleSupport_ ? 0 : N_ds_; 
    for (int j = j_start ; j <= N_ds_ ; j++)
    {
        const double alpha_j = static_cast<double>(j) / static_cast<double>(N_ds_);
        
        Eigen::Vector4d Oij = (offset_cstr_zmp + N_ * (alpha_j * P_supportFoot_0 + (1 - alpha_j) * zmp_));
        A_f.block(0, j , 4, 1) += Oij;
        A_f.block(0, j + 1 , 4, 1) -= Oij;
    }
 
    //Remainings
    for (int i = 1 ; i <= N_steps; i++)
    {
        
        for (int j = 0 ; j <= (i != N_steps ? N_ds_ : N_tdsLast - 1 )  ; j++)
        {
            const double alpha_j = static_cast<double>(j) / static_cast<double>(N_ds_);
        
            Eigen::Vector4d O_ij = offset_cstr_zmp + N_ * (alpha_j * refSteps_[i-1].translation().segment(0,2));
            
            if( i > 1)
            {
                O_ij +=  N_ * ((1 - alpha_j) * refSteps_[i-2].translation().segment(0,2));
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
    timings_constraints(A_Tsteps,b_Tsteps);

    Eigen::MatrixXd A_ineq = Eigen::MatrixXd::Zero(A_Tsteps.rows() + A_f.rows() , N_variables);
    Eigen::VectorXd b_ineq = Eigen::VectorXd::Zero(A_ineq.rows());
    A_ineq <<       - A_f * exp(eta_ * t_)       , A_Tsteps;
    b_ineq << b_f - (N_ * dcm_ ) , b_Tsteps;
    const int NineqCstr = static_cast<int>(A_ineq.rows()); 


    Eigen::MatrixXd A_eq = Eigen::MatrixXd::Zero(0 , N_variables);
    Eigen::VectorXd b_eq = Eigen::VectorXd::Zero(A_eq.rows());
    const int NeqCstr = static_cast<int>(A_eq.rows()); 
    
    //Cost function
    Eigen::MatrixXd M_timings = Eigen::MatrixXd::Identity(N_variables,N_variables);
    Eigen::VectorXd b_timings = Eigen::VectorXd::Ones(M_timings.rows());
    
    double t_im1 = doubleSupport_ ? 0 : refTimings_[0] ;
    double j0 = j_start;
    for(int i = doubleSupport_ ? 0 : 1; i <= NStepsTimings; i++)
    {
        for (int j = 0 ; j <= (i != N_steps ? N_ds_ : N_tdsLast - 1 )  ; j ++)
        {
            double alpha_j = static_cast<double>(j)/static_cast<double>(N_ds_);
            b_timings( (N_ds_ + 1) * i + j) = exp(-eta_ * ( t_im1 + alpha_j * (refTds_)) );
  
        }
        t_im1 = refTimings_[i];
        j0 = 0;
    }

    
    Eigen::MatrixXd Q_cost = betaTsteps * M_timings.transpose() * M_timings 
                            + 1e-12 * Eigen::MatrixXd::Identity(N_variables,N_variables);
    Eigen::VectorXd c_cost = betaTsteps * (-M_timings.transpose() * b_timings) ;

    Eigen::QuadProgDense QP;
    QP.problem(N_variables, NeqCstr,NineqCstr);
    bool QPsuccess = QP.solve(Q_cost, c_cost, A_eq, b_eq, A_ineq, b_ineq);



    if(!QPsuccess)
    {
        Eigen::VectorXd Tds = Eigen::VectorXd::Ones(refTimings_.size()) * (t_s_range_.x() - t_ss_range_.x());
        Eigen::VectorXd Ts = Eigen::VectorXd::Ones(refTimings_.size()) * (t_s_range_.x());
        optimalDoubleSupportTimings_ = std::vector<double>(Tds.data(), Tds.data() + Tds.rows() );
        optimalStepsTimings_ = std::vector<double>(Ts.data(), Ts.data() + Ts.rows() );
        return false;
    }
    
    Eigen::VectorXd solution_ = QP.result();
    
    optimalStepsTimings_.clear();
    optimalDoubleSupportTimings_.clear();
    double tds_i = 0;
    for (int i = 0 ; i < NStepsTimings ; i++)
    {
        tds_i = -log(solution_( (N_ds_+1) * i + N_ds_))/eta_  - -log(solution_( (N_ds_ + 1 ) * i ))/eta_;

        if(i == 0){tds_i += t_;}
        
        optimalStepsTimings_.push_back(-log(solution_( (N_ds_ + 1) * (i+1)))/eta_);
        optimalDoubleSupportTimings_.push_back(tds_i);

    } 
    return true;


}