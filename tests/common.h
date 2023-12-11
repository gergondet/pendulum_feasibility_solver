#include <pendulum_feasibility_solver/feasibility_solver.h>

#include <mc_rtc/Configuration.h>

struct ConfigureParams
{
  double eta;
  double delta;
  Eigen::Vector2d t_ds_range;
  Eigen::Vector2d t_ss_range;
  Eigen::Vector2d t_s_range;
  Eigen::Vector2d stepCstrSize;
  Eigen::Vector2d zmpRange;
  double feetDistance;
  int N_ds;

  static ConfigureParams fromConfiguration(const mc_rtc::Configuration & cfg)
  {
    return {cfg("eta"),          cfg("delta"),    cfg("t_ds_range"),   cfg("t_ss_range"), cfg("t_s_range"),
            cfg("stepCstrSize"), cfg("zmpRange"), cfg("feetDistance"), cfg("N_ds")};
  }
};

struct SolveParams
{
  double t;
  double t_lift;
  bool dbl_supp;
  Eigen::Vector2d dcm;
  Eigen::Vector2d zmp;
  std::string supportFoot;
  sva::PTransformd X_0_supportFoot;
  sva::PTransformd X_0_swingFoot;
  double tds_ref;
  std::vector<sva::PTransformd> steps_ref;
  std::vector<double> timings_refs;
  Eigen::Vector2d gamma;
  double kappa;

  static SolveParams fromConfiguration(const mc_rtc::Configuration & cfg)
  {
    return {cfg("t"),       cfg("t_lift"),      cfg("dbl_supp"),        cfg("dcm"),
            cfg("zmp"),     cfg("supportFoot"), cfg("X_0_supportFoot"), cfg("X_0_swingFoot"),
            cfg("tds_ref"), cfg("steps_ref"),   cfg("timings_refs"),    cfg("gamma"),
            cfg("kappa")};
  }
};

struct Input
{
  ConfigureParams config;
  SolveParams solve;

  static Input fromConfiguration(const mc_rtc::Configuration & cfg) { return {cfg("configure"), cfg("solve")}; }
};

inline std::vector<Input> loadInputs(const std::string & path)
{
  return mc_rtc::Configuration(path).operator std::vector<Input>();
}
