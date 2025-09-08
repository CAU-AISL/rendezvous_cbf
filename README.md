### Cartesian coordinate frames
ECI frame
``` math
\mathcal{F}_i \triangleq \left\{ O_i, \hat{i}_i, \hat{j}_i, \hat{k}_i\right\}
```

Body fixed frame of *chaser*
``` math
\mathcal{F}_c \triangleq \left\{ O_c, \hat{i}_c, \hat{j}_c, \hat{k}_c\right\}
```

Body fixed frame of *target* 
``` math
\mathcal{F}_t \triangleq \left\{ O_t, \hat{i}_t, \hat{j}_t, \hat{k}_t\right\}
```

## Relative motion dynamics
### Relative attitude kinematics and dynamics
**Modified Rodrigues Parameters (MRP)** $\boldsymbol{\sigma} = [\sigma_x, \sigma_y, \sigma_z]^\mathsf{T} \in \mathbb{R}^3$ describes the relative attitude between the chaser and the target. 

``` math
\begin{aligned}
&\dot{\boldsymbol{\sigma}} = G(\boldsymbol{\sigma}) \boldsymbol{\omega} \\
&G(\boldsymbol{\sigma}) = \frac{1}{4}\left[(1 - \boldsymbol{\sigma}^\mathsf{T}\boldsymbol{\sigma} )I_3 + 2\Omega(\boldsymbol{\sigma}) + 2\boldsymbol{\sigma} \boldsymbol{\sigma}^\mathsf{T}\right] 
\end{aligned}
```

$\boldsymbol{\omega} \in \mathbb{R}^3$ denotes **the relative angular velocity between the two spacecraft**. And $R_t^c \in \text{SO}(3)$ represents the rotation matrix from $\mathcal{F}_t$ to $\mathcal{F}_c$.
``` math
\begin{aligned}
&\boldsymbol{\omega} = \boldsymbol{\omega}_{i,c}^c - R_t^c \boldsymbol{\omega}_{i, t}^t \\
&R_t^c = I_3 - \frac{4 (1 - \boldsymbol{\sigma}^\mathsf{T}\boldsymbol{\sigma})}{(1 + \boldsymbol{\sigma}^\mathsf{T} \boldsymbol{\sigma})^2} \Omega (\boldsymbol{\sigma}) + \frac{8 \Omega(\boldsymbol{\sigma})^2}{(1 + \boldsymbol{\sigma}^\mathsf{T}\boldsymbol{\sigma})^2}
\end{aligned}
```

$\boldsymbol{\omega}_{i, c}^c$ : The angular velocity of $\mathcal{F}_c$ relative to $\mathcal{F}_i$, referenced in $\mathcal{F}_c$.

$\boldsymbol{\omega}_{i, t}^t$ : The angular velocity of $\mathcal{F}_t$ relative to $\mathcal{F}_i$, referenced in $\mathcal{F}_t$.

$\Omega(x)$ : Skew symmetric matrix from the vector $x$.

$J_c \in \mathbb{R}^{3\times3}$ is the inertia matrix of the chaser spacecraft, $\tau \in \mathbb{R}^3$ denotes the control torque, and $\tau_d \in \mathbb{R}^3$ denotes the external disturbance torque.

``` math
J_c \dot{\boldsymbol{\omega}} = C_1 (\boldsymbol{\omega}) \boldsymbol{\omega} + D_1(\boldsymbol{\omega}) + \tau + \tau_d
```

where:
``` math
\begin{aligned}
&C_1(\boldsymbol{\omega}) = -J_c \Omega(R_t^c \boldsymbol{\omega}_{i,t}^t) - \Omega(R_t^c \boldsymbol{\omega}_{i,t}^t) J_c + \Omega(J_c(\boldsymbol{\omega} + R_t^c \boldsymbol{\omega}_{i,t}^t)) \\
&D_1(\boldsymbol{\omega}) = -\Omega(R_t^c \boldsymbol{\omega}_{i,t}^t) J_c R_t^c \boldsymbol{\omega}_{i,t}^t - J_c R_t^c \dot{\boldsymbol{\omega}}_{i,t}^t 
\end{aligned}
```

### Relative translation kinematics and dynamics
The relative position vector between the chaser spacecraft and the target represented in the frame $\mathcal{F}_c$:

``` math
\boldsymbol{\rho} = \boldsymbol{r}_c - R_t^c \boldsymbol{r}_t
```

where $\boldsymbol{r}_c \in \mathbb{R}^3$ and $\boldsymbol{r}_t \in \mathbb{R}^3$ denotes the inertial position vectors of the chaser spacecraft and the target, respectively. According to **the fundamental equation of the two-body problem**, the relative position dynamics can be expressed as:
``` math
\begin{aligned}
&\dot{\boldsymbol{\rho}} = \boldsymbol{v} - \Omega(\boldsymbol{\omega}_{i,c}^c)\boldsymbol{\rho}\\
&m_c \dot{\boldsymbol{v}} = m_c C_2\boldsymbol{v} + m_c D_2 + \boldsymbol{f} + \boldsymbol{f}_d 
\end{aligned}
```

where $\boldsymbol{v} \in \mathbb{R}^3$ is the relative velocity between the two spacecraft represented in the $\mathcal{F}_c$.
``` math
\begin{aligned}
&C_2 = -\Omega(\boldsymbol{\omega}_{i,c}^c)\\
&D_2 = -\frac{\mu}{\|\boldsymbol{r}_c\|^3} \boldsymbol{r}_c - R_t^c \dot{\boldsymbol{v}}_t
\end{aligned}
```
