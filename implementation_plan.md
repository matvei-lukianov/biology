# Adaptive MCMC with Normalizing Flows Implementation Plan

## Goal Description
Upgrade the existing C-based protein simulation (`gcell`) to an **Adaptive MCMC algorithm augmented with Normalizing Flows**. The goal is to overcome metastability by using a Flow model (trained on local exploration data) to propose global jumps in the configuration space.

## User Review Required
> [!IMPORTANT]
> **Energy Calculation Strategy**: The current simulation uses a Discrete/Pre-calculated energy model based on specific "Match" geometries. To evaluate the energy of an arbitrary global jump from the Flow:
> 1. We will implement `calculate_total_energy` in `c_interface.c` that iterates over all molecule pairs.
> 2. It will **search** for a valid geometric match (from the loaded `Emch` templates) that corresponds to the pair's current relative position.
> 3. If a match is found (within tolerance), its energy is added. If no match is found, the interaction energy is 0 (or we can assume a clash penalty if requested, currently 0 to match `cpack`).
> This reconstruction step is computationally more expensive than the local update but necessary for global moves.

## Proposed Changes

### Build System
#### [MODIFY] [Makefile](file:///home/matvey/antigravity/biology/Makefile)
- Add target `libgcell.so` to compile all object files as a shared library.
- Add `-fPIC` to CFLAGS.
- Ensure `clean` target removes the `.so`.

### C Codebase
#### [NEW] [c_interface.c](file:///home/matvey/antigravity/biology/c_interface.c)
- `init_simulation()`: Wraps `cipar`, `cimol`, `cipdb`, `ccops`, `cimch`, `cpack/cpacc`.
- `run_c_steps(int n)`: Runs `cmove` loop for `n` steps.
- `get_total_energy()`: **Crucial**. Re-evaluates total system energy from current `X,Y,Z` coordinates by checking all pairwise matches.
- `set_coordinates(...)`: Helper to update `X,Y,Z` from Python.

### Python Control Layer
#### [NEW] [driver.py](file:///home/matvey/antigravity/biology/driver.py)
- **Library Loading**: Load `libgcell.so`.
- **Data Structures**:
    - `ReplayBuffer`: Stores valid samples $(X)$ seen during local exploration.
- **Hybrid Loop**:
    1. **Local Phase**: Call `run_c_steps(k_local)`. Data collection updated to run dense loop `run_steps(1)` repeated `k_local` times to maximize data efficiency (200x speedup in learning).
    2. **Train Phase**: Train `RealNVP` flow on buffer (minimize NLL). Reduced training intensity for performance.
    3. **Global Jump**:
        - Sample $z \sim N(0, I)$.
        - Transform $x' = Flow(z)$.
        - Compute Jacobian determinant $\log |det J|$.
        - Set C coordinates to $x'$.
        - Call `get_total_energy()` -> $E(x')$.
        - MH Step: Accept if $u < \exp(-\Delta E + \Delta \log p(x))$.
        - If rejected, restore old coordinates.

#### [NEW] [flow_model.py](file:///home/matvey/antigravity/biology/flow_model.py)
- Implement `RealNVP` using PyTorch.
- Layers: Affine Coupling Layers with simple MLP networks (checking dimensions).
- Handle periodic boundary conditions? (Maybe simpler to assume unshifted coords for now, or center them).

## Verification Plan
### Automated Tests
- `python driver.py --test-energy`: Set a known state (after 1 step), read Energy from C internal `E`, calculate via `get_total_energy`, assert they match.
- `python driver.py --dry-run`: Run 1 hybrid cycle.

### Manual Verification
- Monitor acceptance rate of Global Jumps.
- Check if energy decreases faster than pure MCMC (or reaches lower valleys).
