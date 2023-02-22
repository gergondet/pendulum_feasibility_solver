#include "../include/pendulum_feasibility_solver/feasibility_solver.h"
#include "../include/pendulum_feasibility_solver/polygons.h"
#include <ScsEigen/Solver.h>



bool feasibility_solver::solve(double t,double t_lift,bool dbl_supp,const Eigen::Vector2d & dcm, const Eigen::Vector2d & zmp , const std::string & supportFoot, const sva::PTransformd & X_0_supportFoot , const sva::PTransformd & X_0_swingFoot,
                   double tds_ref , std::vector<sva::PTransformd> & steps_ref,
                   std::vector<double> & timings_refs)
{

    assert (steps_ref.size() == timings_refs.size());

    doubleSupport_ = dbl_supp;
    tLift_ = t_lift; 
    refSteps_ = steps_ref;
    refTimings_ = timings_refs;
    optimalSteps_ = steps_ref;
    optimalStepsTimings_ = timings_refs;
    supportFoot_ = supportFoot;
    X_0_SupportFoot_ = X_0_supportFoot;
    X_0_SwingFoot_ = X_0_swingFoot;
    dcm_ = dcm;
    zmp_ = zmp;
    refTds_ = tds_ref;
      
    Eigen::VectorXd Tds = Eigen::VectorXd::Ones(refTimings_.size()) * refTds_;
    refDoubleSupportTimings_ = std::vector<double>(Tds.data(), Tds.data() + Tds.rows() );
    t_ = t;
    const Eigen::Matrix2d & R_supportFoot_0 = X_0_SupportFoot_.rotation().transpose().block(0,0,2,2);
    N_ <<  1 ,  0,
         -1 ,  0,
          0 ,  1,
          0 , -1;
    
    N_*= R_supportFoot_0;

    bool feasible_ = true;
    update_x();
    for (int k = 0 ; k < 4 ; k++)
    {
        Eigen::MatrixXd Q;
        Eigen::VectorXd c;
        double r;
        feasibility_constraint(Q,c,r,k);
        feasible_ = feasible_ && (N_ * dcm_)(k) < x_.transpose() * Q * x_ + c.dot(x_) + r;
    } 
    bool ret = true;
    ret = ret && solve_timings();

    // refTimings_ = optimalStepsTimings_; 
    // refDoubleSupportTimings_ = optimalDoubleSupportTimings_;
    // ret = ret && solve_steps();
    // refSteps_ = optimalSteps_;
    // ret = ret && solve_timings();

    return ret; 

}
void feasibility_solver::feasibility_constraint(Eigen::MatrixXd & Q_out, Eigen::VectorXd & c_out, double & r_out,int k)
{

    const int N_tds = static_cast<int>(refTds_ / delta_);
    const int N_tdsLast = static_cast<int>(refTds_ / (2*delta_));
    const int N_steps =  static_cast<int>(refSteps_.size());
    const int N_timings = N_steps;
    const int N_variables = 2 * N_steps + N_timings;

    const Eigen::Vector2d & P_supportFoot_0 = X_0_SupportFoot_.translation().segment(0,2);
    const Eigen::Vector2d & P_swingFoot_0 = X_0_SwingFoot_.translation().segment(0,2);

    ScsEigen::Solver solver;
    solver.mathematicalProgram().setNumberOfVariables(N_variables);

    //Variables are organised as such  : timings then footsteps
    //footstep i indx : N_timings + 2 * i 

    
    //cst part of the offsets
    Eigen::Matrix<double,4,1> offset_cstr_zmp; 
    offset_cstr_zmp << zmpRange_.x()/2,  zmpRange_.x()/2, zmpRange_.y()/2,  zmpRange_.y()/2;

    //DCM must remain inside the feasibility region

    //We generate the cost function for each vertice of the rectangle
    Q_out = Eigen::MatrixXd::Zero(N_variables,N_variables);
    c_out = Eigen::VectorXd::Zero(N_variables);
    r_out = 0;
    const Eigen::Matrix<double,1,2> & n_k = N_.block(k,0,1,2); 
    
    //i = 0 
    const int j_start = static_cast<int>(t_/refTds_)*N_tds;
    for (int j = j_start ; j < N_tds ; j++)
    {
        const double alpha_j = static_cast<double>(j) / static_cast<double>(N_tds);
        const double e_alpha_j_Tds = exp(-eta_ * alpha_j * refTds_); 
        const double alpha_jp1 = static_cast<double>(j+1) / static_cast<double>(N_tds);
        const double e_alpha_jp1_Tds = exp(-eta_ * alpha_jp1 * refTds_); 
        r_out += (offset_cstr_zmp(k) 
                + n_k * (alpha_j * P_supportFoot_0 + (1 - alpha_j) * P_swingFoot_0)) * 
                (e_alpha_j_Tds - e_alpha_jp1_Tds);
    }
    r_out += (offset_cstr_zmp(k) + n_k * P_supportFoot_0 ) * exp(-eta_ * refTds_);
    c_out(0) += - (offset_cstr_zmp(k) + n_k * P_supportFoot_0 );

    //Remainings
    for (int i = 1 ; i <= N_steps; i++)
    {
        
        for (int j = 0 ; j < (i != N_steps ? N_tds : N_tdsLast )  ; j++)
        {
            const double alpha_j = static_cast<double>(j) / static_cast<double>(N_tds);
            const double e_alpha_j_Tds = exp(-eta_ * alpha_j * refTds_); 
            const double alpha_jp1 = static_cast<double>(j+1) / static_cast<double>(N_tds);
            double e_alpha_jp1_Tds = exp(-eta_ * alpha_jp1 * refTds_); 
            if( i == N_steps && j == N_tdsLast - 1)
            {
                e_alpha_jp1_Tds = 0;
            }


            Q_out.block( i - 1 , N_timings + 2 * (i-1),1,2) +=  n_k * alpha_j * (e_alpha_j_Tds - e_alpha_jp1_Tds);

            c_out(i - 1) += offset_cstr_zmp(k) * (e_alpha_j_Tds - e_alpha_jp1_Tds);
            
            if( i > 1)
            {
                Q_out.block( i - 1 , N_timings + 2 * (i-2),1,2) += n_k * (1 - alpha_j) * (e_alpha_j_Tds - e_alpha_jp1_Tds);
            }
            else
            {   
                const double n_P = n_k * P_supportFoot_0;
                c_out( i - 1 ) += n_P * ((1 - alpha_j) * (e_alpha_j_Tds - e_alpha_jp1_Tds));
            }
            
        }
        if(i != N_steps)
        {
            Q_out.block(i     , N_timings + 2 * (i-1),1,2) += - n_k  ;
            Q_out.block(i - 1 , N_timings + 2 * (i-1),1,2) += n_k * exp(-eta_ * refTds_);
            
            c_out(i)     += -offset_cstr_zmp(k);
            c_out(i - 1) += offset_cstr_zmp(k) * exp(-eta_ * refTds_);

            
        }


    }

  
    
}

void feasibility_solver::update_x()
{
    x_ = Eigen::VectorXd::Zero(2*optimalSteps_.size() + optimalStepsTimings_.size());
    x_.segment(0,optimalStepsTimings_.size()) = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(optimalStepsTimings_.data(), optimalStepsTimings_.size());
    for (size_t i = 0 ; i < optimalSteps_.size() ; i++)
    {
        x_.segment(2 * i , 2) = optimalSteps_[i].translation().segment(0,2);
    }
}