# Self-balanced differential robot with sampling-based Model Predictive Control(MPC)

#### Project for the course of Advanced Control Methods at Skoltech, Prof. Pavel Oseninko.

#### Team: Ali Alridha Abdulkarim + Dimitry Nalbersky

### Introduction





### /src
- `sigway_imu.ino` contains the functionalities  of the imu unit measurements, without MPPI. These measurements are used for state estimation of the robot, with the states being mainly: `"angle of inclination"` and `"angular velocity"`, other measurements are also calculated by explicit integration being `"position"` and `"linear velocity"`.

- `sigway_eigen.cpp` bare implementation of the `MPPI (Model Predictive Path Integral)`, just the algorithm on its own.

- `sigway_MPPI.ino` the `"Arduino"` version of `MPPI`, without `Eigen` using only the data structures of `arduino`, and without the `IMU` measurements.

- `sigway_PID.ino` contains straightforward implementation of the `PID` controller with state estimation recieved from the `IMU`.

- `sigway_MPPI_IMU.ino` contains the fully functional `MPPI` code, with state updates recieved from the IMU unit, and being used in the update iterations of the `stochastic MPC`.

### /lib
Contains the required library header files, required for the motros and the `IMU`.
### docs
Contains the research papers and materials we have used and learned from.

Requirements:
- Arduino Uno
- Commercially available arduino motor shield
- IMU (MPU6050)

