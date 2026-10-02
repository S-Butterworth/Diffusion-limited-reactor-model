# Methodology

This document walks through the key equations and modeling decisions behind the pellet-scale diffusion-reaction model and its coupling to the reactor-scale simulation. It assumes familiarity with the governing equations used in reactor modelling. The focus here is on the specific derivations and numerical techniques used in this project.

## 1. Governing Equation

For a spherical catalyst pellet with two reacting species A and B, fractional order in concentration kinetics ($rate ∝ C_A^{1/3}C_B^{2/3}$), the steady-state diffusion-reaction balance for each species is:

$$D_i\left(\frac{d^2C_i}{dr^2} + \frac{2}{r}\frac{dC_i}{dr}\right) = \nu_i \cdot Act \cdot k(T) \cdot C_A^{1/3}C_B^{2/3}$$

Where *r* is the radial coordinate, $D_i$ is the effective diffusivity of species *i*, $ν_i$ is its stoichiometric coefficient, $Act$ is catalyst activity (Section 5), and $k(T)$ follows Arrhenius temperature dependence.

## 2. Normalisation and the Thiele Modulus

### 2.1 Dimensionless concentration and radius

Define $\tilde{C}_i = C_i/C_{i,s}$ (normalised by surface concentration) and $\tilde{r} = r/R$ (normalised by pellet radius). Substituting and applying the chain rule twice:

$$\frac{d^2C_i}{dr^2} = \frac{1}{R^2}\frac{d^2\tilde{C}_i}{d\tilde{r}^2}, \qquad 
\frac{1}{r}\frac{dC_i}{dr} = \frac{1}{R\tilde{r}}\cdot\frac{1}{R}\frac{d\tilde{C}_i}{d\tilde{r}}$$

This produces the dimensionless form used throughout the solver:

$$\frac{d^2\tilde{C}_i}{d\tilde{r}^2} + \frac{2}{\tilde{r}}\frac{d\tilde{C}_i}{d\tilde{r}} = \phi_i^2 \cdot \tilde{C}_A^{1/3}\tilde{C}_B^{2/3}$$

### 2.2 The mass-action normalisation rule

For a general rate law $k\prod_i C_i^{n_i}$, normalising species *X*'s own equation (dividing by $D_X C_{X,s}$) leaves every *other* species' surface concentration at its rate-law exponent unchanged, but reduces *X*'s own exponent by exactly 1:

$$\phi_X^2 \propto \prod_i C_{i,s}^{\,n_i - [i = X]}$$

Applied here: $\phi_A^2 \propto C_{As}^{-2/3}C_{Bs}^{2/3}$ and $\phi_B^2 \propto C_{As}^{1/3}C_{Bs}^{-1/3}$ each depends on composition **only through the ratio** $C_{As}/C_{Bs}$.

## 3. Handling the Singularity at r = 0

The $2/\tilde{r}$ term is singular at the pellet center. Since $C_i'(0)=0$ by spherical symmetry, this is a *removable* singularity as the limit is finite, but naive evaluation divides by zero.

`bvp4c`'s `SingularTerm` option handles this internally: the system is written as $y' = f(x,y) + \frac{S}{x}y$, with $S$ a constant matrix whose entries identify which equations carry a $1/\tilde{r}$ term and with what coefficient (here, $S_{22}=S_{44}=-2$, all other entries zero). The solver uses a collocation scheme aware of this structure near the origin rather than evaluating the division directly.

## 4. Why Effectiveness Factor Is Independent of Conversion for a Stoichiometric Feed

Tracking species flow rates with conversion *X*: $F_A(X) = F_{A0}(1-X)$, $F_B(X) = F_{B0} - \tfrac{1}{2}F_{A0}X$. If the feed is exactly stoichiometric ($F_{A0}=2F_{B0}$), substituting gives $F_B(X) = F_{B0}(1-X)$, **both species scale by the identical factor** $(1-X)$, so their ratio $F_A/F_B$, and hence $C_{As}/C_{Bs}$, is invariant with conversion.

Since (Section 2.2) $\phi_A^2$ and $\phi_B^2$ depend on composition only through this ratio, they and therefore the effectiveness factor η are independent of conversion for this specific feed. This was confirmed numerically (identical η across widely different conversion levels at fixed temperature) before being traced to this algebraic cause. For a non-stoichiometric feed the ratio drifts with X and this independence does not hold, the (T,X) sweep code generalises to that case without modification.

## 5. Catalyst Deactivation

Activity follows the general n-th order deactivation law:

$$-\frac{da}{dt} = k_d a^n \;\Rightarrow\; a(t) = \left[a_0^{1-n} + k_d(n-1)t\right]^{1/(1-n)} \quad (n\neq 1)$$

Reducing to $a(t)=a_0 e^{-k_d t}$ for $n=1$. Activity enters the pellet model multiplicatively on the rate constant ($\phi_i^2 \propto Act\cdot k$), so a lower activity directly reduces the Thiele modulus, a deactivated catalyst is effectively less diffusion-limited, since its intrinsic reaction rate has slowed relative to a fixed diffusion environment.

## 6. The Dead-Core Bifurcation

For reaction orders less than 1, once $\phi_A$ exceeds a critical value, the true solution develops a **dead core**: an interior region $r<r_c$ where $C_A\equiv 0$. This is because the rate law's derivative, $\frac{d}{dC}C^{1/3} = \frac{1}{3}C^{-2/3}$, diverges as $C\to 0$, which is incompatible with the smooth, positive profile the standard BVP assumes.

### 6.1 Approaches used

**Regularised (soft-floor).** The concentration fed into the rate law is smoothed via $C_{\text{used}} = \tfrac{1}{2}(C+\sqrt{C^2+\delta^2})$, which matches $C$ away from zero but keeps the rate law's derivative finite everywhere. This removes the singularity at the cost of a small, controlled approximation near $C=0$.

**Exact free-boundary.** The domain is restricted to the active zone $r\in[r_c,1]$, stretched onto a fixed computational domain via $s=(r-r_c)/(1-r_c)$, with $r_c$ solved for as an unknown parameter. Boundary conditions at $s=0$: $C_A=0$ (continuity with the dead zone) and $C_A'=0$. Since $r$ never reaches zero on this domain, no singular-term handling is required.

A further simplification: inside the dead zone, $C_A\equiv0$ forces the reaction rate to zero regardless of $C_B$, so B's dead-zone equation reduces to source-free diffusion with the unique bounded solution $C_B=\text{const}$, meaning the dead zone need not be solved numerically at all; it contributes only the condition $C_B'(r_c)=0$ to the active-zone problem.

### 6.2 Why onset is numerically difficult

Near the critical $\phi$, $r_c \propto \sqrt{\phi-\phi_{crit}}$ the dead-core radius has an **infinite initial slope** at onset. Newton's method (used internally by `bvp4c`) relies on this slope to size its correction steps. Where the true slope is infinite, the linearisation becomes ill-conditioned. This was observed directly as reciprocal Jacobian condition numbers near $10^{-18}$–$10^{-28}$ when attempting to solve the free-boundary formulation within a narrow band immediately above the critical temperature.

### 6.3 Resolution: bridging via the regularised formulation

The full pellet-scale sweep uses the soft-floor formulation up to and including this narrow onset band, switching to the exact free-boundary formulation only once the estimated dead-core radius clears a minimum threshold, with the switch attempt verified before being committed to, falling back to the regularised formulation otherwise.

## 7. Reactor-Scale Coupling

Pellet-scale results are precomputed into lookup tables η(T, t) (and η(T, X) where composition independence does not hold) rather than solved inline within the reactor ODE, decoupling the expensive, bifurcation-aware pellet solve from the reactor march. The plug-flow reactor model (`Reactor.m`) integrates conversion and temperature along catalyst weight using `ode45`, querying the table via `interp2` at each step; requests outside the table's validated range raise an explicit error rather than extrapolating silently.

Two downstream applications built on this coupling:

- **Temperature compensation**: for a fixed target conversion, the required inlet temperature is solved via root-finding (`fzero`) against the table at each point in the catalyst's service life, producing a temperature-ramping schedule that maintains performance as the catalyst deactivates.
- **Catalyst sizing**: analogous root-finding on catalyst mass, bounded to a finite search range with an explicit unreachability check.