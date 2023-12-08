#include <pendulum_feasibility_solver/feasibility_solver.h>

#include <mc_rtc/Configuration.h>
#include <mc_rtc/clock.h>

#ifdef WITH_VALGRIND_SUPPORT
#  include <valgrind/callgrind.h>
#else
#  define CALLGRIND_START_INSTRUMENTATION
#  define CALLGRIND_STOP_INSTRUMENTATION
#endif

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

int main(int argc, char * argv[])
{
  if(argc < 2)
  {
    mc_rtc::log::error("[usage] {} [inputs.yaml]");
    return 1;
  }
  auto inputs = mc_rtc::Configuration(argv[1]).operator std::vector<Input>();
  feasibility_solver solver;
  auto now = mc_rtc::clock::now();
  CALLGRIND_START_INSTRUMENTATION;
  for(const auto & [cfg, solve] : inputs)
  {
    solver.configure(cfg.eta, cfg.delta, cfg.t_ds_range, cfg.t_ss_range, cfg.t_s_range, cfg.stepCstrSize, cfg.zmpRange,
                     cfg.feetDistance, cfg.N_ds);
    solver.solve(solve.t, solve.t_lift, solve.dbl_supp, solve.dcm, solve.zmp, solve.supportFoot, solve.X_0_supportFoot,
                 solve.X_0_swingFoot, solve.tds_ref, solve.steps_ref, solve.timings_refs, solve.gamma, solve.kappa);
  }
  CALLGRIND_STOP_INSTRUMENTATION;
  auto end = mc_rtc::clock::now();
  mc_rtc::log::info("Took {} milliseconds to run {} times", mc_rtc::duration_ms(end - now).count(), inputs.size());
  return 0;
}
