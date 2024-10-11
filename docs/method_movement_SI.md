## SI: derivation of speed and angular change terms  

Here we derive from the change in velocity the speed and angular change. The change in velocity in 3 dimensions is

$$
\begin{split}
\frac{\text{d}\mathbf{v} _i}{\text{d}t} 
    &= 
    \frac{\text{d}}{\text{d}t} ( v_i\  \hat{\mathbf{e}} _{v, i})
    =
    \frac{\text{d}}{\text{d}t} v_i \left(
        \begin{array}{c}
            \cos(\varphi) \sin(\theta) \\
            \sin(\varphi) \sin(\theta) \\
            \cos(\theta)
        \end{array} \right)\\
    &= 
    \frac{\text{d} v _i}{\text{d}t} \hat{\mathbf{e}} _{v, i}
        + v _i \frac{\text{d} \varphi}{\text{d}t}  \frac{\text{d} \hat{\mathbf{e}} _{v, i}}{\text{d} \varphi}
        + v _i \frac{\text{d} \theta}{\text{d}t}  \frac{\text{d} \hat{\mathbf{e}} _{v, i}}{\text{d} \theta} \\
    &=
    \frac{\text{d} v _i}{\text{d}t} \hat{\mathbf{e}} _{v, i}
        + v _i \frac{\text{d} \varphi}{\text{d}t}  \hat{\mathbf{e}} _{\varphi, i}
        + v _i \frac{\text{d} \theta}{\text{d}t}  \hat{\mathbf{e}} _{\theta, i}
\end{split}
$$

Since the unit vectors [$\hat{\mathbf{e}} _{v, i}$, $\hat{\mathbf{e}} _{\varphi, i}$, $\hat{\mathbf{e}} _{\theta, i}$] are a orthonormal basis, we can compute the change in speed, heading angle and depth angle by:

$$
\begin{split}
\frac{\text{d} v _i}{\text{d}t} 
    &= 
    \frac{\text{d}\mathbf{v} _i}{\text{d}t} \cdot
    \hat{\mathbf{e}} _{v, i} = \\
\frac{\text{d} \varphi _i}{\text{d}t} 
    &= 
    \frac{1}{v_i}
    \frac{\text{d}\mathbf{v} _i}{\text{d}t} \cdot
    \hat{\mathbf{e}} _{\varphi, i}  \\
\frac{\text{d} \theta _i}{\text{d}t} 
    &= 
    \frac{1}{v_i}
    \frac{\text{d}\mathbf{v} _i}{\text{d}t} \cdot
    \hat{\mathbf{e}} _{\theta, i}  \\
\end{split}
$$

In the following we compute $\frac{\text{d} \varphi}{\text{d} t}$ which can be repeated analogously for $\theta$:


$$
\begin{split}
&\frac{\text{d} \varphi _i}{\text{d}t} 
    = 
    \frac{1}{v_i}
    \frac{\text{d}\mathbf{v} _i}{\text{d}t} \cdot
    \hat{\mathbf{e}} _{\varphi, i}
    =
    \frac{ 1 }{m v_i}
    \left( \mathbf{F} _i \cdot \hat{\mathbf{e}} _{\varphi, i} - \alpha
            \frac{\text{d} \varphi _i}{\text{d}t} 
    \right) \\
\rightarrow\ 
&\frac{\text{d} \varphi _i}{\text{d}t}
    \underbrace{
    \left( 
        1 + \frac{\alpha}{m v _i}
    \right)
    } _{\frac{m v _i - \alpha}{m v _i}}
    =
    \frac{ \mathbf{F} _i }{m v_i}
     \cdot \hat{\mathbf{e}} _{\varphi, i}  \\
\rightarrow\ 
&\frac{\text{d} \varphi _i}{\text{d}t}
    =
    \frac{ \mathbf{F} _i }{m v _i + \alpha}
     \cdot \hat{\mathbf{e}} _{\varphi, i} \ .
\end{split}
$$


## SI: Notes on the rotational friction coefficient

Note that we assume here the overdamped case, i.e. we neglect the moment of inertia due to dominating turning friction.

At first glance it seems odd, that the turn friction still exists, even if the mass is set to zero.
However, the introduced turn friction relaxes the point-like particle assumption and is meant to represent the friction of the surrounding fluid with an object of non-zero volume.
That means the friction parameter rather relates to the surface of the object, than to its mass.
Thus, an extremely light object with $m \rightarrow 0$ that still has an surface area $A \gg 0$ does still experience finite turning speed.

Of course, the surface area for objects of similar shape and density, is a function of the mass, and for spherical particle, the relation is

$$ 
\begin{split}
    m &= \hat{\rho} V = \rho \frac{4}{3} \pi r ^3 \\
    A &= 4 \pi r ^2 \\
    \rightarrow A &\propto m ^{2/3} \ .
\end{split}
$$