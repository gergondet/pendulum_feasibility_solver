# pendulum_feasibility_solver

Linear Inverted Pendulum based footsteps location and steps timing adptation for Bipedal Locomotion.

Using the current walking plan and the robot state, this will generate the closest location and timings from the plan that will guarantee a feasible LIPM-based centroidal trajectory

# Dependencies

[SpaceVecAlg](https://github.com/jrl-umi3218/SpaceVecAlg)

[eigen-quadprog](https://github.com/jrl-umi3218/eigen-quadprog)

## Installation
In repo directory
```bash
mkdir build && cd build
cmake ..
make
sudo make install
