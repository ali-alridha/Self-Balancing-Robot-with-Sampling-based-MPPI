# Self-balanced differential robot with Sampling-Based Model Predictive Control (SBMPC)

#### Project for the course of Advanced Control Methods at Skoltech, Prof. Pavel Oseninko.

#### Team: Ali Alridha Abdulkarim + Dimitry Nalbersky

### Introduction
## Controller : MPPI Controller

### Note
The following MPPI implementation follows Algorithms 1 and 2 of the reference paper. 

### Reference
1. G. Williams et al. "Information-Theoretic Model Predictive Control: Theory and Applications to Autonomous Driving" 
    - URL : https://ieeexplore.ieee.org/document/8558663
    - PDF : https://arxiv.org/pdf/1707.02342.pdf

### Brief overview of MPPI algorithm
General process flow to calculate optimal input with mppi algorithm goes as:

**[Step 1]** ramdomly sample input sequence

Mean input sequence $U$ and ramdomly sampled input sequence $V$ are defined as follows.  
Usually, optimal input sequence on the previous step is used as $U$.
```math
$$
    \begin{align}
        & (\mathbf{u}_0, \mathbf{u}_1, ... \mathbf{u}_{T-1}) = U \in \mathbb{R}^{m \times T}, \nonumber \\
        & (\mathbf{v}_0, \mathbf{v}_1, ... \mathbf{v}_{T-1}) = V \in \mathbb{R}^{m \times T}, \nonumber \\
        & \mathbf{v}_t = \mathbf{u}_t + \epsilon_t, \nonumber \\
        & \epsilon_t \sim \mathcal{N}(0, \Sigma).\nonumber 
    \end{align}
$$
```

**[Step 2]** predict future states and evaluate cost for each sample

We assume a discrete time, continuous state-action dynamical system as a control target.  
$\mathbf{x}$ is a system state, and $\mathbf{v}$ is a sampled control input.
```math
$$
\begin{align}
\mathbf{x}_t  &\in \mathbb{R}^{n}, \nonumber \\
\mathbf{x}_{t+1} &= \mathbf{F}(\mathbf{x}_t, \mathbf{v}_t).\nonumber 
\end{align}
$$
```

Then costs (i.e. penalties to be minimized) for sampled sequences $S(V; \mathbf{x}_0)$ can be evaluated with following formulations.
```math
$$
    \begin{align}
        & S(V; \mathbf{x}_0) = C(\mathcal{H}(V; \mathbf{x}_0)), \nonumber \\
        & C(\mathbf{x}_0, \mathbf{x}_1, ... \mathbf{x}_T) = \phi(\mathbf{x}_T) + \sum_{t=0}^{T-1}c(\mathbf{x}_t), \nonumber \\
        & \mathcal{H}(V; \mathbf{x}_0) = \left( \mathbf{x}_0, \mathbf{F}(\mathbf{x}_0, \mathbf{v}_0), \mathbf{F}(\mathbf{F}(\mathbf{x}_0, \mathbf{v}_0), \mathbf{v}_1), ... \right).\nonumber 
    \end{align}
$$
```

**[Step 3]** calculate weight for each sample sequence

Weight for a each sample sequence is derived on the basis of information theory.  
There are K sample sequences in total, represented with an index k.  
Good control sequence with small cost value get more weight, and vice versa.  
```math
$$
\begin{align}
& w(V) = \frac{1}{\eta} \exp
\left( 
    -\frac{1}{\lambda}
    \left(
        S(V) + \lambda(1-\alpha) \sum^{T-1}_{t=0} \mathbf{u}_t^T \Sigma^{-1} (\mathbf{u}_t + \epsilon_t) - \rho
    \right)
\right) \nonumber \\
& \eta = 
\sum_{k=1}^K \exp
\left( 
    -\frac{1}{\lambda}
    \left(
        S(U + \mathcal{E}_k) + \lambda(1-\alpha) \sum^{T-1}_{t=0} \mathbf{u}_t^T \Sigma^{-1} (\mathbf{u}_t + \epsilon_t^k) - \rho
    \right)
\right)\nonumber \\
& \rho = 
\min_k 
\left( S(V_k) + \lambda(1-\alpha) \sum^{T-1}_{t=0} \mathbf{u}_t^T \Sigma^{-1} (\mathbf{u}_t + \epsilon_t^k) \right)\nonumber
\end{align}
$$
```
Note that $\rho$ is inserted into the formulation to avoid overflow errors during implementation.

**[Step 4]** get optimal control input sequence

Finally, optimal input trajectory for the next ($i+1$) step is given adding weighted sample sequences to the previous solution.
```math
$$
\begin{align}
    \mathbf{u}_t^{i+1} % &= \mathbb{E}_{\mathbb{Q}_{\hat{U}, \Sigma}}[w(V)\mathbf{v}_t]
                 = u_t^i + \sum_{k=1}^K w(V_k) \epsilon_t^k \nonumber 
\end{align}
$$
```



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

