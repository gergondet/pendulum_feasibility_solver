#include "common.h"

#include <mc_rtc/clock.h>

#ifdef WITH_VALGRIND_SUPPORT
#  include <valgrind/callgrind.h>
#else
#  define CALLGRIND_START_INSTRUMENTATION
#  define CALLGRIND_STOP_INSTRUMENTATION
#endif

int main(int argc, char * argv[])
{
  if(argc < 2)
  {
    mc_rtc::log::error("[usage] {} [inputs.yaml]");
    return 1;
  }
  auto inputs = loadInputs(argv[1]);
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
