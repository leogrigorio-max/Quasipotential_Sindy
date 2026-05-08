
Repository with codes to implement the method for learning the drift decomposition into gradient and orthogonal components from instanton trajectories using the SINDy algorithm
according to the paper Quasi-Potential of Stochastic Dynamics via Instanton-Guided
Sparse Identification.


To run the sparse identification step, execute:

    'sindy.m'

This script applies a Stepwise Sparse Regressor (SSR) version of SINDy, adapted here to recover the gradient and orthogonal components of the drift from instanton data.

The instanton trajectories used by the code are stored in the folder:

    'trajectories'

If you would like to generate your own instanton data, go to the 'instanton_evaluation' folder and follow the steps below:

1. Run:

    'initialize.m'

This script initializes the required variables.

2. Then run:

    'repar.m'

This script computes instantons by the Chernykh–Stepanov (CS) method.

Time integration is performed with a symplectic scheme based on the midpoint rule, namely

    x - x_ = h * H((x + x_)/2, (p + p_)/2)_p
    p - p_ = h * H((x + x_)/2, (p + p_)/2)_x

Good practices

1. We recommend starting with small values of lambda (the final condition for p), and then progressively increasing lambda to reach points farther away from the limit cycle.

2. The folder

    'NR_functions'

contains the Newton–Raphson routines used by repar.m, since the symplectic discretization is implicit.
