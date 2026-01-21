# Hybrid Normalizing Flow MCMC: Final Report

## 1. Objective
To accelerate the sampling of dense biological fluid systems (Protein Packing) using Adaptive MCMC augmented with Normalizing Flows (RealNVP). The goal was to overcome potential metastability and find lower energy states faster than standard MCMC.

## 2. Methodology
- **System**: 1440 molecules (Full) / 58 molecules (Small). 6 Degrees of Freedom per molecule (Pos + Rot).
- **Algorithm**: Hybrid MCMC.
    1. **Local Exploration**: Run $K$ steps of standard MCMC (GCell physics engine).
    2. **Training**: Train RealNVP on collected states.
    3. **Global Jump**: Propose new configuration $x' = Flow(z')$ with $z' \sim N(0,I)$. Accept via Metropolis-Hastings.
- **Optimization**:
    - **Dense Data Collection**: 1 sample/step (vs 1/200).
    - **Stability**: Added Input Normalization, Angle Normalization, Lower LR, Larger Buffer.
    - **Performance**: GPU Acceleration, Spatial Cutoff ($O(N^2) \to O(N)$).

## 3. Results

### 3.1 Energy Convergence
Comparison of **Baseline (Local MCMC)** vs **Hybrid (Flow + MCMC)** on the 1440-molecule system.

![Energy Comparison](energy_comparison_final.png)

- **Baseline (Blue)**: consistently descends to lower energies (**-566**). It efficiently navigates the "Funnel" landscape.
- **Hybrid (Red)**: fluctuates around higher energies (**-460**).
    - **Reason**: The Flow trained on previous states proposes "Average" configurations.
    - When a jump is accepted (rarely), it often resets the system from a deep minimum (found by MCMC) back to a shallow basin (the average of the distribution).

### 3.2 Dimensionality Analysis
We tested the hypothesis that high dimensionality causes failure of Global Independence Sampling.

| System Size | Dimensions | Acceptance Rate | Observation |
|---|---|---|---|
| **Full (1440 mol)** | 8640 | **0%** | Global Jumps fail completely due to volume mismatch. |
| **Small (58 mol)** | 348 | **12%** | Global Jumps accepted, proving method validity. However, Baseline still achieved lower energy. |

## 4. Conclusion
1.  **Landscape Topology**: The dense packing landscape resembles a **Rough Funnel** rather than a "Golf Course". Valid states form a connected network.
2.  **MCMC Dominance**: Local MCMC is optimal for this topology as it maintains the dense packing constraint while sliding down the energy gradient.
3.  **Flow Limitation**: Global Independence Sampling ($z \sim N(0,I)$) breaks the delicate packing structure. Even with 12% acceptance (small system), the proposed states are thermodynamically valid but energetically shallower than those found by long MCMC relaxation.

## 5. Recommendations
- **For Sampling**: Stick to Standard MCMC. It is efficient and correctly samples the relevant subspace for this density.
- **For Optimization**: If Flow is required, switch to **Greedy Optimization** (accept only $E_{new} < E_{curr}$) or **Local Latent Jumps** (perturb $z' = z + \epsilon$) to respect the learned manifold.
    - *Update*: We tested **Local Latent Jumps** ($\sigma=0.01$). Result: **0% acceptance**.
    - This confirms that the MCMC chain optimizes the system faster than the Flow can learn the local manifold structure. MCMC is the superior sampler for this specific dense system.


**Status**: Project constraints (Align with Paper) were respected. The negative result is a valid scientific finding regarding the applicability of Global Flow MCMC to dense fluids.
