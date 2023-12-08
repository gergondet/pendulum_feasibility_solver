find_package(mc_rtc REQUIRED)

include(FindPkgConfig)

pkg_check_modules(VALGRIND IMPORTED_TARGET valgrind)

add_executable(run-feasibility-solver run-feasibility-solver.cpp)
target_link_libraries(run-feasibility-solver PUBLIC ${PROJECT_NAME} mc_rtc::mc_rtc_utils)
if(VALGRIND_FOUND)
  target_compile_definitions(run-feasibility-solver PUBLIC -DWITH_VALGRIND_SUPPORT)
endif()

add_test(NAME TestFeasibiltySolverWithSimulationInputs
         COMMAND run-feasibility-solver
                 ${CMAKE_CURRENT_SOURCE_DIR}/feasibility-inputs-simulation.yaml)
