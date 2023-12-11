#include "common.h"

#include <benchmark/benchmark.h>

extern std::string bench_inputs();

static auto inputs = loadInputs(bench_inputs());

static void BM_FeasibilitySolver(benchmark::State & state)
{
  feasibility_solver solver;
  size_t idx = 0;
  for(auto _ : state)
  {
    const auto & [cfg, solve] = inputs[idx];
    solver.configure(cfg.eta, cfg.delta, cfg.t_ds_range, cfg.t_ss_range, cfg.t_s_range, cfg.stepCstrSize, cfg.zmpRange,
                     cfg.feetDistance, cfg.N_ds);
    solver.solve(solve.t, solve.t_lift, solve.dbl_supp, solve.dcm, solve.zmp, solve.supportFoot, solve.X_0_supportFoot,
                 solve.X_0_swingFoot, solve.tds_ref, solve.steps_ref, solve.timings_refs, solve.gamma, solve.kappa);
    idx++;
    if(idx == inputs.size()) { idx = 0; }
  }
}

BENCHMARK(BM_FeasibilitySolver)->Unit(benchmark::kMicrosecond)->Repetitions(10);

BENCHMARK_MAIN();
