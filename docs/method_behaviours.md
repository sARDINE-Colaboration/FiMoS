## Behavioural states

### State-specific movement

We assume the fish to have a repertoire of $K$ behavioral states. The movement parameters vary between states which are defined by the following state-specific parameters

$$
\begin{itemize}
   
    \item speed relaxation coefficient $ \beta $,
     \item preferred speed $ v _0 $,
    \item diffusion coefficients $ D _{\varphi} $, $ D _{\theta} $, and $ D _v $, 
    \item patch attraction strength $ \mu _{\mathrm{patchatt}} $,
    \item social force parameters $ \mu _{\mathrm{att}} $ and $ \mu _{\mathrm{ali}} $, representing attraction and alignment strengths,
\end{itemize}
$$

resulting in a parameter vector $\Theta^{(k)}$ for each state $ k \in \{1, \dots, K\} $:

$$
\Theta^{(k)} = \left( 
\beta ^{(k)}, \; 
v _0 ^{(k)}, \; 
D _{\varphi} ^{(k)}, \; 
D _{\theta} ^{(k)}, \; 
D _v ^{(k)}, \; 
\mu _{\mathrm{patchatt}} ^{(k)}, \; 
\mu _{\mathrm{att}} ^{(k)}, \; 
\mu _{\mathrm{ali}}^{(k)}
\right)
$$




### State-switching dynamics

A simple Markov chain is used for simulating the state dynamics.
At each second, a draw from the transition probability matrix $\Gamma$ determines what the next state is going to be. 
For each state $ k $, a sojourn time $ \tau _k $  is specified, indicating the average time (in minutes) spent in this state before switching. If a fish switches state, it is equally likely to switch to any of the other activated states. This means that the sojourn times $ \bm{\tau} = (\tau _1, \dots, \tau _K) $ solely govern the state transitions and the duration of staying in the same state follows a geometric distribution with the mean given by the respective sojourn time.

The transition probability matrix $ \Gamma $ is defined as:

$$
\Gamma = \left(
\begin{array}{ccc}
\gamma _{11} & \cdots & \gamma _{1K} \\
\vdots & \ddots & \vdots \\
\gamma _{K1} & \cdots & \gamma _{KK}
\end{array}
\right)
$$

Where each diagonal element (probability of staying in a state) is given by:

$$
\gamma _{kk} = \frac{1}{\tau _k \cdot 60} 
$$

And the off-diagonal elements for $k, l \in \{1, \dots, K\}$ and $k \neq l$ are:

$$
\gamma _{kl} = \left(1 - \gamma _{kk}\right) \cdot \frac{1}{K-1}
$$

### Parameter Variation

The simulator allows for variability in state-dependent movement parameters between fish with $K$ parameter vectors for each fish $i$: $\Theta^{(k, i)}$ is created by rescaling $\hat{\Theta}^{(k, i)} \sim \mathcal{N}(\tilde{\Theta}^{(k)}, \sigma^2)$ with $\tilde{\Theta}^{(k)}$ being the user-provided state-dependent movement parameters and $\sigma$ the user-provided standard deviation. The scaling ensures that the parameters in $\Theta^{(k, i)}$ lie between predefined lower and upper limits, $L$ and $U$, respectively:

$$
 L = \left( 0.1, 0.01, 0.01, 0.01, 0.01, 0, 0, 0 \right), 
 U = \left( 3, 3, 1, 1, 1, 10, 5, 5 \right)
$$

These limits are chosen to make sure the parameters stay within a biologically realistic range and are applied as follows:

$$
\Theta^{(k,i)} =
\left\{
  \begin{array}{ll}
    L, & \text{if } \hat{\Theta}^{(k,i)} < 0, \\
    U, & \text{if } \hat{\Theta}^{(k,i)} > 1, \\
    \hat{\Theta}^{(k,i)} \cdot (\bm{U} - \bm{L}) + \bm{L}, & else.
  \end{array}
\right.
$$


Additional to the variation between fish, variation of state-specific parameters within fish is introduced by the idea of substates.
Each fish has a repertoire of $M _k$ substates within each behavioral state $ k $, and the movement parameters vary between substates. These realizations $\hat{\Theta}^{(k, m _k, i)} \sim \mathcal{N}(\tilde{\Theta}^{(k)}, \sigma^2)$ are drawn for each substate $m _k$ and fish $i$ once at the beginning of a simulation and redrawn every time any parameter in the GUI is changed. 
Each time a fish switches states, one of the substates is randomly chosen that defines the movement parameters of the fish. 


### Interaction in the GUI

By default, only the resting state is activated when the application starts.
For the predefined states, namely resting, foraging, and active, there are defaults for the vectors of state-dependent unscaled movement parameter means $\tilde{\Theta}^{(k)}$, the sojourn times $ \tau _k $, and the number of substates $M _k$, $k \in [1,2,3]$. The global standard deviation $\sigma$ also has a default.

The unscaled parameter means, the sojourn times, the number of substates and the standard deviation can be adjusted in the GUI. 
Note for the choice of $\tilde{\Theta}^{(k)}$ and $\sigma$ that draws $\hat{\Theta}^{(k,i)}$ outside of the interval $[0,1]$ will be set to according to the parameter limits.

The number of states $K$ can be changed in the GUI by activating and deactivating behavioural states. Additional states can be added manually. 
When no behavioral state is selected, the fish freeze while the time keeps running.