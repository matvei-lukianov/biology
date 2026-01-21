# Hybrid Flow-MCMC Algorithm

**Input:**
1. Initialize variables from current state $mol(t)$:
   - Center positions $x^0_l$, atom offsets $x^i_l$, energies $e_l$, orientation matrices $U_l$ for $l=1 \dots N_{mol}$.
   - Flow Model parameters $\theta$, Replay Buffer $\mathcal{R}$, Empty Queue $Q$, Busy Set $B$.

**Algorithm Loop:**

1. **Local Exploration Phase (Parallel MCMC)**
   Run for $K$ steps. In each step:
   
   a. Randomly sample $M$ molecule indices. Store in queue $Q$.
   b. **While** $Q$ is not empty:
      i.   Wait for free thread.
      ii.  Take $L$ from head of $Q$.
      iii. **If** $L$ is a potential neighbor to any molecule in $B$:
             Move $L$ to tail of $Q$ (defer).
           **Else**:
             Add $L$ to $B$ and start thread:
             
             1. Sample random neighbor $LL$ for $L$. If none, remove $L$ from $B$, continue.
             2. Sample random match $k$ (docking configuration) between $L$ and $LL$.
             3. Form putative state $x'_L$ using match $k$.
             4. **Validate** move length: If $\sum ||x'_L - x_L|| > d_{max}$, reject.
             5. **Check Collision**: If hard sphere overlap, reject.
             6. **Metropolis Balance**:
                Calculate $p = \min(1, \exp(-\frac{E_{after} - E_{before}}{T}) \frac{N_{after}}{N_{before}})$.
                With probability $p$, accept:
                  - Detach $L$ from old neighbors ($e \leftarrow e - e_{old}$).
                  - Attach $L$ to $LL$ ($e \leftarrow e + e_{match}^k$).
                  - Update coordinates $x_L \leftarrow x'_L$, $U_L \leftarrow U'_L$.
             7. Remove $L$ from $B$.

   c. Store valid states $x = \{x^0_l, U_l\}_{l=1}^{N_{mol}}$ into Replay Buffer $\mathcal{R}$.

2. **Flow Training Phase**
   **If** $|\mathcal{R}| > N_{batch}$:
   a. Sample batch $X \sim \mathcal{R}$.
   b. Compute Loss: $\mathcal{L} = -\frac{1}{|X|} \sum_{x \in X} (\log p_Z(f_\theta(x)) + \log |\det \frac{\partial f}{\partial x}|)$.
   c. Update $\theta \leftarrow \theta - \eta \nabla_\theta \mathcal{L}$.

3. **Global Jump Phase (Flow Proposal)**
   **If** Training is sufficiently converged (or Pretraining active):
   
   a. Sample latent vector $z' \sim \mathcal{N}(0, I)$ (Global) OR $z' = f_\theta(mol(t)) + \epsilon$ (Local Latent).
   b. Generate candidate $x' = f_\theta^{-1}(z')$.
   c. Compute Jacobian correction: $\Delta \log J = \log q(mol(t)) - \log q(x')$.
   d. Compute Energy Change: $\Delta E = E(x') - E(mol(t))$.
   e. **Metropolis Acceptance**:
      Calculate $\alpha = \min(1, \exp(-\frac{\Delta E}{T} + \Delta \log J))$.
      Sample $u \sim U[0,1]$.
      **If** $u < \alpha$:
        $mol(t) \leftarrow x'$ (Accept Jump).
      **Else**:
        $mol(t) \leftarrow mol(t)$ (Reject).

**Output:**
Updated state $mol(t+1)$.
