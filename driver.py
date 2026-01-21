import ctypes
import numpy as np
import torch
import torch.optim as optim
import os
import sys
import time
import argparse
from flow_model import RealNVP

# --- C interface configuration ---
class CInterface:
    def __init__(self, lib_path="./libgcell.so"):
        self.lib = ctypes.CDLL(os.path.abspath(lib_path))
        
        # Define argtypes/restypes
        self.lib.init_simulation_interface.restype = ctypes.c_int
        self.lib.calculate_total_energy.restype = ctypes.c_float
        
        # Global pointers
        self.lib.get_ptr_mnum.restype = ctypes.POINTER(ctypes.c_int)
        self.lib.get_ptr_nat.restype = ctypes.POINTER(ctypes.c_int)
        self.lib.get_ptr_E.restype = ctypes.POINTER(ctypes.c_float)
        self.lib.get_ptr_X.restype = ctypes.POINTER(ctypes.c_float)
        self.lib.run_c_steps.argtypes = [ctypes.c_int]
        self.lib.get_nmol.restype = ctypes.c_int
        
        # Stats getter
        self.lib.get_stats_interface.argtypes = [ctypes.POINTER(ctypes.c_float)]
        
        # Rotation functions
        self.lib.get_rotations_euler.argtypes = [
            ctypes.POINTER(ctypes.c_float),  # alpha
            ctypes.POINTER(ctypes.c_float),  # beta
            ctypes.POINTER(ctypes.c_float),  # gamma
            ctypes.c_int                      # n
        ]
        self.lib.set_rotations_euler.argtypes = [
            ctypes.POINTER(ctypes.c_float),  # alpha
            ctypes.POINTER(ctypes.c_float),  # beta  
            ctypes.POINTER(ctypes.c_float),  # gamma
            ctypes.c_int                      # n
        ]
        
        # Dimensions (Must match cdim.h)
        self.DIMA = 10000 
        self.DIMN = 10
        self.DIMC = 2100
        self.DIMM = self.DIMN * self.DIMC # 21000 
        # But we can just use the pointer and rely on nmol to access valid range.
        # However, for numpy mapping we need stride.
        # X is float[DIMM+2][DIMA+1].
        # In C, this is a single block of size (DIMM+2)*(DIMA+1)*4 bytes.
        
        # We need to map this carefully.
        # Let's verify DIMA from cdim.h (I read it earlier: DIMA 10000)
        # So stride is 10001 floats (0 index included?)
        # cdim.h: float X[DIMM+2][DIMA+1]
        self.STRIDE_A = self.DIMA + 1
        
    def init(self):
        res = self.lib.init_simulation_interface()
        if res != 0:
            raise RuntimeError(f"C initialization failed with code {res}")
        self.nmol = self.lib.get_nmol()
        print(f"Initialized C simulation with {self.nmol} molecules.")
        
        # Setup numpy views
        self._setup_views()
        
    def _setup_views(self):
        # pointers
        ptr_X = self.lib.get_ptr_X()
        # ptr_Y = self.lib.get_ptr_Y() # Need to expose Y and Z too!
        # Accessing global symbols directly if not exposed via function:
        # ctypes provides in_dll for globals
        
        # Load globals
        self.c_X = ctypes.c_float.in_dll(self.lib, "X")
        self.c_Y = ctypes.c_float.in_dll(self.lib, "Y")
        self.c_Z = ctypes.c_float.in_dll(self.lib, "Z")
        self.c_nat = ctypes.c_int.in_dll(self.lib, "nat") # This is nat[DIMN+2], array of atom counts per TYPE
        self.c_mnum = ctypes.c_int.in_dll(self.lib, "mnum") # mnum[DIMM+2], molecule type ID per molecule
        
        # Create numpy wrappers
        # Current shape: (DIMM+2, DIMA+1)
        # We don't know DIMM exactly at runtime unless we trust cdim.h or read a constant.
        # But we only access up to nmol.
        # Let's map a large enough buffer.
        
        total_size = (self.nmol + 2) * self.STRIDE_A
        # We can map the whole thing
        # But ctypes multidimensional arrays are tricky.
        # Let's cast pointer to array
        
        ArrayType = ctypes.c_float * total_size
        
        self.np_X = np.ctypeslib.as_array(ArrayType.from_address(ctypes.addressof(self.c_X)))
        self.np_Y = np.ctypeslib.as_array(ArrayType.from_address(ctypes.addressof(self.c_Y)))
        self.np_Z = np.ctypeslib.as_array(ArrayType.from_address(ctypes.addressof(self.c_Z)))
        
        self.np_X = self.np_X.reshape(-1, self.STRIDE_A)
        self.np_Y = self.np_Y.reshape(-1, self.STRIDE_A)
        self.np_Z = self.np_Z.reshape(-1, self.STRIDE_A)
        
        # Map mnum and nat
        MnumType = ctypes.c_int * (self.DIMM + 2)
        self.np_mnum = np.ctypeslib.as_array(MnumType.in_dll(self.lib, "mnum"))
        
        NatType = ctypes.c_int * (self.DIMN + 2)
        self.np_nat = np.ctypeslib.as_array(NatType.in_dll(self.lib, "nat"))

    def get_stats(self):
        stats = (ctypes.c_float * 4)()
        self.lib.get_stats_interface(stats)
        return list(stats) # [Eavg, Savg, AccAvg, MSDavg]

    def get_coords(self):
        # Return state vector for Flow: (nmol, 3) Center of Mass
        # 1-based indexing in C (1..nmol)
        x = self.np_X[1:self.nmol+1, 0]
        y = self.np_Y[1:self.nmol+1, 0]
        z = self.np_Z[1:self.nmol+1, 0]
        return np.stack([x, y, z], axis=1) # (nmol, 3)

    def set_coords(self, new_coords):
        # new_coords: (nmol, 3)
        # Update CM and shift atoms
        
        current_x = self.np_X[1:self.nmol+1, 0]
        current_y = self.np_Y[1:self.nmol+1, 0]
        current_z = self.np_Z[1:self.nmol+1, 0]
        
        dx = new_coords[:, 0] - current_x
        dy = new_coords[:, 1] - current_y
        dz = new_coords[:, 2] - current_z
        
        # Update CM
        self.np_X[1:self.nmol+1, 0] = new_coords[:, 0]
        self.np_Y[1:self.nmol+1, 0] = new_coords[:, 1]
        self.np_Z[1:self.nmol+1, 0] = new_coords[:, 2]
        
        # Update atoms
        # This loop operation in numpy is vectorizable if nat is constant, but nat varies.
        # Loop over molecules (faster than C loop for < 1000 molecules, ok for python)
        for i in range(self.nmol):
            mid = i + 1
            mtype = self.np_mnum[mid]
            n_atoms = self.np_nat[mtype]
            
            # 1..n_atoms
            self.np_X[mid, 1:n_atoms+1] += dx[i]
            self.np_Y[mid, 1:n_atoms+1] += dy[i]
            self.np_Z[mid, 1:n_atoms+1] += dz[i]

    def run_steps(self, n):
        self.lib.run_c_steps(ctypes.c_int(n))
        
    def get_energy(self):
        return self.lib.calculate_total_energy()
    
    def get_rotations(self):
        """Get Euler angles (alpha, beta, gamma) for all molecules. Returns (nmol, 3)."""
        alpha = (ctypes.c_float * self.nmol)()
        beta = (ctypes.c_float * self.nmol)()
        gamma = (ctypes.c_float * self.nmol)()
        self.lib.get_rotations_euler(alpha, beta, gamma, ctypes.c_int(self.nmol))
        return np.stack([
            np.ctypeslib.as_array(alpha),
            np.ctypeslib.as_array(beta),
            np.ctypeslib.as_array(gamma)
        ], axis=1)  # (nmol, 3)
    
    def set_rotations(self, euler_angles):
        """Set Euler angles from (nmol, 3) array of [alpha, beta, gamma]."""
        alpha = euler_angles[:, 0].astype(np.float32)
        beta = euler_angles[:, 1].astype(np.float32)
        gamma = euler_angles[:, 2].astype(np.float32)
        self.lib.set_rotations_euler(
            alpha.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            beta.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            gamma.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            ctypes.c_int(self.nmol)
        )
    
    def get_state(self):
        """Get full state for Flow: (nmol, 6) = [x, y, z, alpha, beta, gamma]."""
        coords = self.get_coords()      # (nmol, 3)
        rotations = self.get_rotations() # (nmol, 3)
        return np.concatenate([coords, rotations], axis=1)  # (nmol, 6)
    
    def set_state(self, state):
        """Set full state from (nmol, 6) array."""
        coords = state[:, :3]
        rotations = state[:, 3:]
        self.set_coords(coords)
        self.set_rotations(rotations)

# --- Replay Buffer ---
class ReplayBuffer:
    def __init__(self, capacity=2000):
        self.buffer = []
        self.capacity = capacity
        
    def add(self, x):
        if len(self.buffer) >= self.capacity:
            self.buffer.pop(0)
        self.buffer.append(x)
        
    def sample(self, batch_size):
        if len(self.buffer) < batch_size:
            return torch.tensor(np.array(self.buffer), dtype=torch.float32)
        indices = np.random.choice(len(self.buffer), batch_size)
        return torch.tensor(np.array([self.buffer[i] for i in indices]), dtype=torch.float32)

# --- Main logic ---
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--steps", type=int, default=50, help="Number of hybrid steps")
    parser.add_argument("--local_steps", type=int, default=200, help="C steps per hybrid step")
    parser.add_argument("--train_steps", type=int, default=10, help="Training steps per hybrid step")
    parser.add_argument("--baseline", action="store_true", help="Run in baseline mode (no flow)")
    parser.add_argument("--output", type=str, default="simulation_results.csv", help="Output CSV file")
    parser.add_argument("--load", type=str, default=None, help="Load model weights from file")
    parser.add_argument("--save", type=str, default="flow_model.pt", help="Save model weights to file")
    parser.add_argument("--pretrain", action="store_true", help="Pretrain mode (save periodically)")
    args = parser.parse_args()

    # Init
    sim = CInterface()
    sim.init()
    
    nmol = sim.nmol
    print(f"System ready. N_mol={nmol}")
    
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")
    
    # Model: Input dim = nmol * 6 (3 positions + 3 Euler angles)
    if not args.baseline:
        flow = RealNVP(input_dim=nmol*6, hidden_dim=512).to(device)
        if args.load:
            print(f"Loading model from {args.load}")
            flow.load_state_dict(torch.load(args.load, map_location=device))
        optimizer = optim.Adam(flow.parameters(), lr=1e-4) # Reduced LR for stability
    else:
        flow = None
        optimizer = None
        
    buffer = ReplayBuffer(capacity=10000) # Increased capacity for better history
    # Normalization Helpers
    # Coords: [0, 500] -> [-1, 1]
    # Angles: [-pi, pi] -> [-1, 1]
    
    def normalize(state):
        # state: (nmol, 6)
        state_norm = state.copy()
        # XYZ
        state_norm[:, :3] = (state[:, :3] - 250.0) / 250.0
        # Angles (Assuming approx -pi to pi range from C)
        state_norm[:, 3:] = state[:, 3:] / np.pi 
        return state_norm

    def denormalize(state_norm):
        state = state_norm.copy()
        # XYZ
        state[:, :3] = state_norm[:, :3] * 250.0 + 250.0
        # Angles
        state[:, 3:] = state_norm[:, 3:] * np.pi
        return state
        
    # Initial state
    current_state = sim.get_state()
    buffer.add(normalize(current_state).flatten())
    
    start_time = time.time()
    accepted_jumps = 0
    
    print(f"{'Step':<6} | {'Loss':<8} | {'Jump':<5} | {'AvgE':<8} | {'Phys(s)':<7} | {'Trn(s)':<6} | {'Jmp(s)':<6} | {'Tot(s)':<6}")
    print("-" * 80)
    
    # Open CSV for logging
    with open(args.output, "w") as f_out:
        f_out.write("Step,FlowLoss,Accepted,AvgEnergy,AvgShift,AvgAcc,AvgMSD,T_Phys,T_Train,T_Jump,T_Total\n")
        
        for step in range(args.steps):
            t_start_step = time.time()
            t_phys = 0.0
            t_train = 0.0
            t_jump = 0.0
            # 1. Local Exploration with Dense Data Collection
            # Collect data every step to fill buffer faster (200x speedup in data generation)
            for _ in range(args.local_steps):
                sim.run_steps(1)
                state = sim.get_state()
                buffer.add(normalize(state).flatten())
                
            flat_state = normalize(state).flatten() # Last state for proposal
            
            current_E = sim.get_energy()
            t_phys = time.time() - t_start_step
            
            # 2. Train Flow (Only if not baseline)
            t_train_start = time.time()
            loss_val = 0.0
            if not args.baseline and len(buffer.buffer) > 32:
                flow.train()
                for _ in range(args.train_steps):
                    batch = buffer.sample(32).to(device)
                    optimizer.zero_grad()
                    loss = -flow.log_prob(batch).mean()
                    loss.backward()
                    optimizer.step()
                    loss_val = loss.item()
                
                # Save periodically
                if args.save and step % 10 == 0:
                    torch.save(flow.state_dict(), args.save)
            t_train = time.time() - t_train_start
            
            # 3. Global Jump Proposal (Only if not baseline)
            accept = False
            t_jump_start = time.time()
            if not args.baseline:
                # Propose
                flow.eval()
                with torch.no_grad():
                    # Current state z
                    tensor_state = torch.tensor(flat_state, dtype=torch.float32).unsqueeze(0).to(device)
                    z_curr, log_det_curr = flow.forward(tensor_state)
                    
                    # Global Independence Sampler
                    z_new = torch.randn_like(z_curr)
                    x_new, _ = flow.inverse(z_new)
                    x_new_np = x_new.cpu().detach().numpy().reshape(nmol, 6)  # (nmol, 6)
                    x_new_np = denormalize(x_new_np) # Scale back to [0, 500] for Physics
                    
                    # Calculate Flow Probabilities
                    log_q_curr = flow.log_prob(tensor_state).item()
                    log_q_new  = flow.log_prob(x_new).item()
                    
                # Calculate Energy of New State
                old_state = sim.get_state()  # Save full state
                sim.set_state(x_new_np)      # Apply new positions + rotations
                new_E = sim.get_energy()
                
                # Metropolis Criterion
                T_val = 100.0 # Matching cpar.gc
                
                delta_E = new_E - current_E
                log_alpha = -delta_E / T_val + (log_q_curr - log_q_new)
                
                if np.log(np.random.rand()) < log_alpha:
                    accept = True
                    accepted_jumps += 1
                    current_E = new_E
                    # Keep new coords (already set)
                else:
                    # Revert full state
                    sim.set_state(old_state)
                
            t_jump = time.time() - t_jump_start

            # Get Stats
            stats = sim.get_stats()
            # Stats: Eavg, Savg, AccAvg, MSDavg

            t_total = time.time() - t_start_step
            cumulative_step = step * args.local_steps
            print(f"{step:<6} | {loss_val:<8.1e} | {'ACC' if accept else 'REJ'} | {stats[0]:<8.1f} | {t_phys:<7.1f} | {t_train:<6.1f} | {t_jump:<6.1f} | {t_total:<6.1f}")
            f_out.write(f"{cumulative_step},{loss_val},{1 if accept else 0},{stats[0]},{stats[1]},{stats[2]},{stats[3]},{t_phys},{t_train},{t_jump},{t_total}\n")
            f_out.flush()
            
    print(f"Finished. Total accepted jumps: {accepted_jumps}/{args.steps}")

if __name__ == "__main__":
    main()
