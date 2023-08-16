#include "../include/pendulum_feasibility_solver/feasibility_solver.h"
#include "../include/pendulum_feasibility_solver/polygons.h"
// #include <ScsEigen/Solver.h>



bool feasibility_solver::solve(double t,double t_lift,
                    bool dbl_supp,
                    const Eigen::Vector2d & dcm, const Eigen::Vector2d & zmp , 
                    const std::string & supportFoot, const sva::PTransformd & X_0_supportFoot , const sva::PTransformd & X_0_swingFoot,
                    double tds_ref , std::vector<sva::PTransformd> steps_ref,
                    std::vector<double> timings_refs)
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

    N_tdsLast = static_cast<int>(static_cast<double>(N_ds_) / 2);
    N_steps =  static_cast<int>(refSteps_.size());
    N_timings = N_steps;
      
    Eigen::VectorXd Tds = Eigen::VectorXd::Ones(refTimings_.size()) * refTds_;
    Eigen::VectorXd Tds_min = Eigen::VectorXd::Ones(refTimings_.size()) * (t_s_range_ - t_ss_range_).x();
    refDoubleSupportDuration_ = std::vector<double>(Tds.data(), Tds.data() + Tds.rows() );
    // optimalDoubleSupportDuration_ = std::vector<double>(Tds_min.data(), Tds_min.data() + Tds_min.rows() );
    optimalDoubleSupportDuration_ = std::vector<double>(Tds.data(), Tds.data() + Tds.rows() );
    // optimalStepsTimings_.clear();
    // for(int i = 0 ; i < N_timings ; i++)
    // {
    //    optimalStepsTimings_.push_back( (i+1) * t_s_range_.x() );
    // }
    t_ = std::max(0. , t);
    Eigen::Matrix2d R_0_rect = X_0_SupportFoot_.rotation().block(0,0,2,2);
    if(doubleSupport_)
    {
        Eigen::Vector2d y1 = (X_0_supportFoot.translation() - X_0_swingFoot.translation()).segment(0,2);
        y1.normalize();
        Eigen::Vector2d x1 = -(Eigen::Vector3d{0,0,1}.cross(Eigen::Vector3d{y1(0),y1(1),0})).segment(0,2);
        R_0_rect.block(0,0,1,2) = x1.transpose();
        R_0_rect.block(1,0,1,2) = y1.transpose();
    }
    N_ <<  1.  ,  0.,
           0.  , -1.,
          -1.  ,  0.,
           0.  ,  1.;
          
    N_*= R_0_rect;

    xStep_ = Eigen::VectorXd::Zero(2 * N_steps);
    xTimings_ = Eigen::VectorXd::Zero(N_timings * (N_ds_ + 1) + N_tdsLast);
    double t_im1 =  0;
    for (int j = 0 ; j <= N_ds_  ; j ++)
    {
        
        double alpha_j = static_cast<double>(j)/static_cast<double>(N_ds_);
        if(doubleSupport_)
        {
            xTimings_(j) = exp(-eta_ * ( t_ + alpha_j * (tds_ref - t_)) );
        }
        else
        {
            xTimings_(j) = exp(-eta_ * t_ );
        }

    }
    for(int i = 0 ; i <= N_steps ; i++)
    {
        if(i!=N_steps)
        {
            xStep_.segment(2*i,2) = optimalSteps_[i].translation().segment(0,2);
        }
        if(i != 0)
        {
            t_im1 = optimalStepsTimings_[i-1];
            for (int j = 0 ; j <= (i != N_steps ? N_ds_ : N_tdsLast - 1 )  ; j ++)
            {
                
                double alpha_j = static_cast<double>(j)/static_cast<double>(N_ds_);
                xTimings_( (N_ds_ + 1) * i + j) = exp(-eta_ * ( t_im1 + alpha_j * tds_ref) );

            }
        }
    }

    //cst part of the offsets
    offsetCstrZMP_ << zmpRange_.x()/2,  zmpRange_.y()/2, zmpRange_.x()/2,  zmpRange_.y()/2;
    const double l = (X_0_supportFoot.translation() - X_0_SwingFoot_.translation()).norm();

    offsetCstrZMPDblInit_ << zmpRange_.x()/2,  l/2, zmpRange_.x()/2,  l/2;

    bool ret = true;
    Niter_ = 0;
    ret = ret && solve_steps(refSteps_);
    Niter_ += 1;
    ret = ret && solve_timings(refTimings_,refTds_);
    Niter_ += 1;
    ret = ret && solve_steps(refSteps_);
    Niter_ += 1;
    ret = ret && solve_timings(refTimings_,refTds_);

    // if(!ret)
    // {
    //     ret = true;
    //     Niter_ += 1;
    //     ret = ret && solve_steps(refSteps_);
    //     Niter_ += 1;
    //     ret = ret && solve_timings(refTimings_,refTds_);
    // }

    

    bool feasible_ = true;

    // update_x();
    // Eigen::Vector4d feasibilityOffset = Eigen::Vector4d::Zero();
    // for (int k = 0 ; k < 4 ; k++)
    // {
    //     Eigen::MatrixXd Q;
    //     Eigen::VectorXd c;
    //     double r;
    //     feasibility_constraint(Q,c,r,k);
    //     feasibilityOffset(k) = (x_.transpose() * Q * x_ + c.dot(x_) + r) * exp(eta_ * t_);
    //     feasible_ = feasible_ && (N_ * dcm_)(k) < feasibilityOffset(k) ;
    // } 
    // Polygon feasibilityPolygon = Polygon(N_,feasibilityOffset);
    // feasibilityRegion_ = feasibilityPolygon.Get_Polygone_Corners();

    return ret; 

}
void feasibility_solver::feasibility_constraint(Eigen::MatrixXd & Q_out, Eigen::VectorXd & c_out, double & r_out,int k)
{


    const int N_variables = 2 * N_steps + N_timings * (N_ds_ + 1) + N_tdsLast;

    const Eigen::Vector2d & P_supportFoot_0 = X_0_SupportFoot_.translation().segment(0,2);
    const Eigen::Vector2d & P_swingFoot_0 = X_0_SwingFoot_.translation().segment(0,2);

    // ScsEigen::Solver solver;
    // solver.mathematicalProgram().setNumberOfVariables(N_variables);

    //Variables are organised as such  : timings then footsteps
    //footstep i indx : N_timings + 2 * i 

    //DCM must remain inside the feasibility region

    //We generate the cost function for each vertice of the rectangle
    Q_out = Eigen::MatrixXd::Zero(N_variables,N_variables);
    c_out = Eigen::VectorXd::Zero(N_variables);
    r_out = 0;
    const Eigen::Matrix<double,1,2> & n_k = N_.block(k,0,1,2); 
    
    //i = 0 
    const int j_start = doubleSupport_ ? 0 : N_ds_;    
    for (int j = 0 ; j <= N_ds_ ; j++)
    {
        const double alpha_j = static_cast<double>(j) / static_cast<double>(N_ds_);
        double Ojk = (offsetCstrZMP_(k) + n_k * (alpha_j * P_supportFoot_0 + (1 - alpha_j) * zmp_));
        c_out(j) += Ojk;
        c_out(j+1) -= Ojk;
    }

    //Remainings
    for (int i = 1 ; i <= N_steps; i++)
    {
        const int step_indx_im1 = N_timings * (N_ds_ + 1) + N_tdsLast + 2 * (i-1);
        const int step_indx_im2 = N_timings * (N_ds_ + 1) + N_tdsLast + 2 * (i-2);
        for (int j = 0 ; j < (i != N_steps ? N_ds_ : N_tdsLast )  ; j++)
        {
            const double alpha_j = static_cast<double>(j) / static_cast<double>(N_ds_);
            const int mu_indx_j = (N_ds_ + 1) * i + j;
            const int mu_indx_jp1 = (N_ds_ + 1) * i + j + 1;


            Q_out.block( mu_indx_j   , step_indx_im1 ,1,2) +=  n_k * alpha_j;
            c_out(mu_indx_j) += offsetCstrZMP_(k);
            if( !(i == N_steps && j == N_tdsLast - 1))
            {
                Q_out.block( mu_indx_jp1 , step_indx_im1 ,1,2) -=  n_k * alpha_j;
                c_out(mu_indx_jp1) -= offsetCstrZMP_(k);
            }

            
            if( i > 1)
            {
                Q_out.block( mu_indx_j   , step_indx_im2 ,1,2) +=   n_k * (1 - alpha_j);
                Q_out.block( mu_indx_jp1 , step_indx_im2 ,1,2) -=   n_k * (1 - alpha_j);
            }
            else
            {   
                const double n_P = n_k * P_supportFoot_0;
                c_out( mu_indx_j)   +=   n_P * (1 - alpha_j);
                c_out( mu_indx_jp1) -=   n_P * (1 - alpha_j);
            }
            
        }
  
        
    }
    Eigen::MatrixXd Q_outT = Q_out.transpose();
    Q_out += Q_outT;


  
    
}

void feasibility_solver::update_x()
{
    const int N_variables = 2 * N_steps + N_timings * (N_ds_ + 1) + N_tdsLast;
    x_ = Eigen::VectorXd::Zero(N_variables);
    
    
    double tds_i = optimalDoubleSupportDuration_[0];
    for(size_t j =  0  ; j <=  static_cast<size_t>(N_ds_) ; j++)
    {        
        const double alpha_j = static_cast<double>(j) / static_cast<double>(N_ds_);
        double t_ij = t_ + (tds_i - t_) * alpha_j;
        if(!doubleSupport_)
        {
            t_ij = t_;
        }

        x_( j ) = exp(-eta_ * t_ij);
    }
    
    for (size_t i = 1 ; i <= static_cast<size_t>(N_timings) ; i++)
    {
        double t_i = optimalStepsTimings_[i - 1];
        tds_i = i < static_cast<size_t>(N_timings) ? optimalDoubleSupportDuration_[i] : tds_i;
        for(size_t j = 0 ; j <= (i != static_cast<size_t>(N_timings) ? static_cast<size_t>(N_ds_) : static_cast<size_t>(N_tdsLast - 1)) ; j++)
        {
            const double alpha_j = static_cast<double>(j) / static_cast<double>(N_ds_);
            double t_ij = t_i + tds_i * alpha_j;
            x_( (N_ds_ + 1) * i + j ) = exp(-eta_ * t_ij);
        }
       
        
    }
    for (size_t i = 0 ; i < static_cast<size_t>(N_steps); i++)
    {
        x_.segment( static_cast<Eigen::Index>( N_timings * (N_ds_ + 1) + N_tdsLast + 2 * i ) , 2) = optimalSteps_[i].translation().segment(0,2);
    }
}

std::vector<Eigen::Vector3d> feasibility_solver::get_feasibility_region(const sva::PTransformd & X_0_supportFoot , const sva::PTransformd & X_0_swingFoot)
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
    
    build_time_feasibility_matrix(A_f,b_f,X_0_supportFoot,X_0_swingFoot);
    
    Eigen::Vector4d feasibilityOffset = exp(eta_ * t_) * ( A_f.block(0,0,N_.rows(),N_variables - N_slack) * xTimings_ + b_f);
    // std::cout << "[Pendulum feasibility solver][Timing solver] output offset " << std::endl << feasibilityOffset << std::endl;
    Polygon feasibilityPolygon = Polygon(N_,feasibilityOffset);
    return feasibilityPolygon.Get_Polygone_Corners();

}