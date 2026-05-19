# Quasi-Potential SINDy

This repository implements a sparse-identification workflow for learning the drift decomposition of stochastic dynamics from instanton trajectories. The method uses SINDy/Stepwise Sparse Regression (SSR) to recover the gradient and transverse components of the drift field associated with the quasi-potential.

The code accompanies the project:

**Quasi-Potential of Stochastic Dynamics via Instanton-Guided Sparse Identification**

## Authors

- Leo Grigorio
- Mnerh Alqahtani

## Mathematical background

For a stochastic dynamical system with drift $b(x)$, the drift is decomposed as

$$
b(x) = -\nabla V(x) + \ell(x),
$$

where $V(x)$ is the quasi-potential and $\ell(x)$ is the transverse/orthogonal component. Along instanton trajectories, the momentum variable satisfies

$$
p = 2\nabla V(x).
$$

Therefore, the sparse-regression step can be used to learn either the gradient component

$$
-\nabla V(x) = -\frac{p}{2},
$$

or the transverse component

$$
\ell(x) = \dot{x} - \frac{p}{2}.
$$

## Repository structure

```text
.
├── sindy.m                  # Original sparse-identification script
├── sindy_fct.m              # Function version of the SINDy/SSR workflow
├── demo_sindy.m             # Example script for running sindy_fct
├── trajectories/            # Precomputed instanton trajectories
├── instanton_evaluation/    # Scripts for generating instantons
├── NR_functions/            # Newton-Raphson routines
├── CITATION.cff             # Citation metadata
└── LICENSE                  # Software license
```

## Requirements

- MATLAB R2025a or newer

## Quick start

Run the demo script from the root folder of the repository:

```matlab
demo_sindy
```

The demo uses the trajectory file

```matlab
file = 'trajectories/limit_cycle_lambda_-15.000.mat';
```

and calls `sindy_fct.m` for both components:

```matlab
results_gradient = sindy_fct(file, 'gradient', true);
results_transverse = sindy_fct(file, 'transverse', true);
```

The third input controls plotting. Use `true` to show the cross-validation plot and `false` to run without plotting.

## Method summary

The sparse-identification step follows the paper notation

$$
\mathbb{Y} \approx \mathbb{X}c,
$$

where $\mathbb{Y}$ is the target data matrix, $\mathbb{X}$ is the candidate library matrix, and $c$ is the sparse coefficient matrix.

In the MATLAB implementation, the same relation is written as

$$
F \approx \Theta \Xi,
$$

where $F$ corresponds to $\mathbb{Y}$, $\Theta$ corresponds to $\mathbb{X}$, and $\Xi$ corresponds to $c$.

The SSR procedure starts from a least-squares fit and progressively removes small coefficients. Cross-validation is used to choose the sparsity level that balances accuracy and model simplicity.

## Reproducing the results

Precomputed instanton trajectories are stored in:

```text
trajectories/
```

To reproduce the sparse-identification results, run:

```matlab
demo_sindy
```

To generate new instanton trajectories, go to:

```text
instanton_evaluation/
```

and run:

```matlab
initialize
repar
```

The script `initialize.m` initializes the required variables, and `repar.m` computes instantons using the Chernykh-Stepanov method. The implicit midpoint scheme uses Newton-Raphson routines stored in `NR_functions/`.

## Practical recommendation for generating new instanton trajectories

Start with small values of the final momentum parameter $\lambda$, then progressively increase $\lambda$ to reach points farther away from the limit cycle.

## Output

For Example 4.4, the limit-cycle example, running `demo_sindy` with the `gradient` option identifies the gradient component of the drift decomposition. In this case, the target vector field is

$$
F = -p/2,
$$

where $p = 2\nabla V(x)$.

The function prints the selected sparsity level, the validation score, and the sparse coefficients for the two components of the learned vector field. For the gradient run, the nonzero terms give

$$
F_1(x,y) \approx x - x^3 - xy^2,
\qquad
F_2(x,y) \approx y - x^2y - y^3.
$$

Equivalently,

$$
F(x,y) \approx
\begin{pmatrix}
x(1-x^2-y^2) \\
y(1-x^2-y^2)
\end{pmatrix}
= -\nabla V(x,y).
$$

When plotting is enabled, the function also displays the cross-validation score as a function of the sparsity level.

## Citation

If you use this code, please cite the repository using the information in `CITATION.cff`.
