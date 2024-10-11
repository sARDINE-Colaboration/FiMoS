# Methods on fish movement

The movements of individual fishes is modeled by a model with inertia introduced for 2D-movement in [Klamser & Gómez-Nava et al. (2021)](https://doi.org/10.3389/fphy.2021.715996). The velocity vector $\mathbf{v} _i$ of individual $i$ changes over time according to

$$
\frac{\text{d}\mathbf{v} _i}{\text{d}t} =  \frac{1}{m}
\left( \mathbf{F} _i - \alpha
    \left(
        \frac{\text{d} \varphi _i}{\text{d}t} \hat{\mathbf{e}} _{\varphi, i} +
        \frac{\text{d} \theta _i}{\text{d}t} \hat{\mathbf{e}} _{\theta, i}
    \right)
\right)
$$

where $\mathbf{F} _i$ is the force acting on fish $i$, $m$ is its mass, $\alpha$ is the turning friction parameter, $\varphi _i$ is its heading angle in xy-plane, $\theta _i$ is its depth angle,
$\hat{\mathbf{e}} _{\varphi, i}$ is the unit vector pointing perpendicular to the heading direction in the xy-plane and 
$\hat{\mathbf{e}} _{\theta, i}$ is the unit vector pointing perpendicular to the heading direction in the depth-plane.

To understand the second term in the equation, the turning friction, we split (see ) the change in velocity into a (i) speed $v _i$ and (ii) heading angle $\varphi _i$ change:

$$
\begin{split}
\frac{\text{d} v _i}{\text{d} t} &=  \frac{\mathbf{F} _{i}}{m}  \cdot \hat{\mathbf{e}} _{v, i} \\
\frac{\text{d} \varphi _i}{\text{d} t} &= \frac{\mathbf{F} _{i}}{m v _i + \alpha} \cdot \hat{\mathbf{e}} _{\varphi, i} \\
\frac{\text{d} \theta _i}{\text{d} t} &= \frac{\mathbf{F} _{i}}{m v _i + \alpha} \cdot \hat{\mathbf{e}} _{\theta, i} \ .
\end{split}
$$

We see that the turning friction has no effect on the change in speed, but alters the change in heading and depth angle such that it does not allow instantaneous turning at $v_i=0$.

The force $\mathbf{F} _i$ acting on individual $i$ has a self-propulsion, an environmental and a social term:

$$
\mathbf{F} _i =
      \mathbf{F} _{i,\text{sp}}
    + \mathbf{F} _{i, \text{env}}
    + \mathbf{F} _{i, \text{soc}}
$$

The self-propulsion force regulates the individual's current speed to the preferred speed $v _0$, and adds fluctuations to speed $v$ and angular speeds $\dot{\varphi}$, $\dot{\theta}$:

$$
\mathbf{F} _{i, \text{sp}} = 
\left(
    \beta (v _0 - v _i) + \sqrt{2 D _v} \xi _{v}
\right) \hat{\mathbf{e}} _{i, v}+ 
\left(
    \sqrt{2 D _{\varphi}} \xi _{\varphi} 
\right) \hat{\mathbf{e}} _{i, \varphi}+ 
\left(
    \sqrt{2 D _{\theta}} \xi _{\theta} 
\right) \hat{\mathbf{e}} _{i, \theta}
$$

with $\beta$ as the speed relaxation coefficient, $\xi _v$, $\xi _{\varphi}$, $\xi _{\theta}$ as independent Gaussian white noise processes and $D _v$, $D _{\varphi}$ and $D _{\theta}$ as diffusion coefficients which set the noise intensity.

The environmental forces modulate the individual's response to the water surface, to the waters' ground and its shore line: 

$$
\mathbf{F} _{i, \text{env}} = \sum _{k \in [\text{surface, ground, shore}]} \mathbf{F} _{i, k} 
$$

with a force depending on the collision time $t _{i, \text{col}}$ with the respective border

$$
\mathbf{F} _{i, k} =
    \begin{cases}
        \text{sgn}(\hat{\mathbf{e}} _{i, x} \cdot \hat{\mathbf{n}} _{i, \text{B}}) \mu _k \hat{\mathbf{e}} _{i, x} &\text{ if } t _{i, \text{col}} \leq t _c \\
        0 &\text{ if } t _{i, \text{col}} > t _c
    \end{cases}
$$

where the signum $\text{sgn}(\dots)$ function ensures that the force points away from the border, $\hat{\mathbf{n}} _{i, \text{B}}$ is the normal vector at the
collision point from the border and the $x$ in $\hat{\mathbf{e}} _{i, x}$ is $\varphi$ for shore-repulsion and $\theta$ otherwise.
In order to prevent a forth and back switching of the signum function due to the individual approaching a corner, the signum function is locked at first collision detection with a border.
Of course, the lock is released as soon as no collision is detected anymore $t _{i, \text{col}} > t _c$.


The social forces !!!!!!!!!!!!!!!!!!!!MISSING!!!!!!!!!!!!!!!!!!!!