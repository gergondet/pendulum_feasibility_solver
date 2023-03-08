#pragma once
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <assert.h>
#include <eigen-quadprog/QuadProg.h>
#include <eigen-quadprog/eigen_quadprog_api.h>
#include <SpaceVecAlg/SpaceVecAlg>

class feasibility_solver
{
    public:

        feasibility_solver() = default;
        ~feasibility_solver(){};

        void configure(const double eta,const double delta , const Eigen::Vector2d t_ds_range ,const Eigen::Vector2d t_ss_range, const Eigen::Vector2d t_s_range , const Eigen::Vector2d stepCstrSize, 
                      const Eigen::Vector2d zmpRange, const double feetDistance, const int N_ds)
        {
            delta_ = delta;
            eta_ = eta;
            t_ds_range_ = t_ds_range;
            t_ss_range_ = t_ss_range;
            t_s_range_ = t_s_range;
            stepCstrSize_ = stepCstrSize;
            feetDistance_ = feetDistance;
            zmpRange_ = zmpRange;
            N_ds_ = N_ds;
        }
        
        /**
         * @brief Compute the steps timing and position to remain in feasibility
         * 
         * @param t actual timing, must be between 0 and the next steps timings (tds0 + tss0)
         * @param dcm current pendulum dcm
         * @param X_0_supportFoot 
         * @param tds_ref 
         * @param steps_ref 
         * @param timings_refs timing of the future steps
         * @return true the solution has been found and the dcm is inside the found feasibility region
         * @return false no solution has been found or the dcm is not inside the found feasibility region
         */
        bool solve(double t,double t_lift,bool dbl_supp,const Eigen::Vector2d & dcm, const Eigen::Vector2d & zmp , const std::string & supportFoot, const sva::PTransformd & X_0_supportFoot , const sva::PTransformd & X_0_swingFoot,
                   double tds_ref , std::vector<sva::PTransformd> steps_ref,
                   std::vector<double> timings_refs);
        
        bool solve_timings(const std::vector<double> & refTimings, const double & refTds);

        bool solve_steps(const std::vector<sva::PTransformd> & refSteps);

        const std::vector<sva::PTransformd> & get_optimal_footsteps()
        {
            return optimalSteps_;
        }

        const std::vector<Eigen::Vector3d> & get_feasibility_region()
        {
            return feasibilityRegion_;
        }

        const std::vector<double> & get_optimal_steps_timings()
        {
            return optimalStepsTimings_;
        }

        const std::vector<double> & get_optimal_steps_ds_duration()
        {
            return optimalDoubleSupportDuration_;
        }
    
    private:

        void timings_constraints(Eigen::MatrixXd & A_out, Eigen::VectorXd & b_out, const int NStepsTimings);

        void kinematics_contraints(Eigen::MatrixXd & A_out, Eigen::VectorXd & b_out, const std::vector<sva::PTransformd> & refSteps);

        /**
         * @brief Compute quadrtic cstr such as
         * 
         * (N Pu) (k) <= xT Q x + cT x + r
         * 
         * @param Q_out 
         * @param c_out 
         * @param r_out 
         * @param k
         */
        void feasibility_constraint(Eigen::MatrixXd & Q_out, Eigen::VectorXd & c_out, double & r_out,int k);

        void update_x();
    
        int N_tdsLast;
        int N_steps;
        int N_timings;
        int Niter_ = 0;
        
        double delta_ = 5e-3;
        double eta_ = 3.4;
        double feetDistance_ = 0.2;
        Eigen::Vector2d t_ds_range_ = Eigen::Vector2d::Zero(); //min_max
        Eigen::Vector2d t_ss_range_ = Eigen::Vector2d::Zero(); //min_max
        Eigen::Vector2d t_s_range_ = Eigen::Vector2d::Zero(); //min_max
        Eigen::Vector2d stepCstrSize_ = Eigen::Vector2d::Zero(); //size of the kinematics cstr rectangle
        Eigen::Vector2d zmpRange_ = Eigen::Vector2d::Zero(); //size of the zmp cstr rectangle
        double betaSteps = 1;
        double betaTsteps = 1;

        double refTds_ = 0.2;
        bool doubleSupport_ = true;

        int N_ds_ = 3;

        sva::PTransformd X_0_SupportFoot_ = sva::PTransformd::Identity();
        sva::PTransformd X_0_SwingFoot_ = sva::PTransformd::Identity();

        double t_ = 0;
        double tLift_ = 0; //time when contact has been released
        Eigen::Matrix<double,4,2> N_;
        Eigen::Matrix<double,4,1> offsetCstrZMP_;
        Eigen::Matrix<double,4,1> offsetCstrZMPDblInit_; 
        Eigen::Vector2d dcm_;
        Eigen::Vector2d zmp_;
        std::string supportFoot_;
        Eigen::VectorXd x_; //desicion variables : steps timing then footsteps pose

        std::vector<sva::PTransformd> refSteps_;
        std::vector<double> refTimings_;
        std::vector<double> refDoubleSupportDuration_;

        std::vector<Eigen::Vector3d> feasibilityRegion_;

        Eigen::VectorXd solution_; 

        std::vector<sva::PTransformd> optimalSteps_;
        std::vector<double> optimalStepsTimings_;
        std::vector<double> optimalDoubleSupportDuration_;

        Eigen::VectorXd xStep_; //current decision variables for step QP
        Eigen::VectorXd xTimings_; //current decision variables for timing QP



};