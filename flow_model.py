import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np

class RealNVP(nn.Module):
    def __init__(self, input_dim, hidden_dim=512, n_layers=4):
        super().__init__()
        self.input_dim = input_dim
        self.layers = nn.ModuleList()
        self.mask = self._get_mask(input_dim)
        
        for _ in range(n_layers):
            self.layers.append(CouplingLayer(input_dim, hidden_dim, self.mask))
            self.mask = 1 - self.mask # Alternate mask

    def _get_mask(self, dim):
        mask = torch.zeros(dim)
        mask[::2] = 1 # Simple checkerboard-like mask for 1D vector
        return mask

    def forward(self, x):
        log_det_sum = 0
        z = x
        for layer in self.layers:
            z, log_det = layer(z)
            log_det_sum += log_det
        return z, log_det_sum

    def inverse(self, z):
        x = z
        log_det_sum = 0
        for layer in reversed(self.layers):
            x, log_det = layer.inverse(x)
            log_det_sum += log_det
        return x, log_det_sum

    def sample(self, n_samples):
        z = torch.randn(n_samples, self.input_dim)
        x, _ = self.inverse(z)
        return x

    def log_prob(self, x):
        z, log_det = self.forward(x)
        # log p(x) = log p(z) + log |det J|
        # p(z) is unit gaussian
        log_pz = -0.5 * (z**2 + np.log(2 * np.pi)).sum(dim=1)
        return log_pz + log_det

class CouplingLayer(nn.Module):
    def __init__(self, input_dim, hidden_dim, mask):
        super().__init__()
        self.register_buffer('mask', mask)
        
        # Simple MLP for the affine transformation parameters
        self.net = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.LeakyReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.LeakyReLU(),
            nn.Linear(hidden_dim, input_dim * 2), # Outputs s and t
        )
        
    def forward(self, x):
        # x: [batch, dim]
        masked_x = x * self.mask
        out = self.net(masked_x)
        s, t = out.chunk(2, dim=1)
        
        # Scaling factor: exp(s) * (1 - mask)
        # Translation: t * (1 - mask)
        # Apply only to masked-out elements
        
        s = s * (1 - self.mask)
        t = t * (1 - self.mask)
        
        # Tanh trick to stabilize s training
        s = torch.tanh(s) 
        
        z = x * torch.exp(s) + t # Affine transform: z = x * exp(s) + t
        # But wait, RealNVP definition:
        # y_d = x_d
        # y_D = x_D * exp(s(x_d)) + t(x_d)
        # where d are masked IN (kept) and D are masked OUT (transformed)
        
        # If mask is 1 for kept, 0 for transformed:
        # z = x_kept + (x_transformed * exp(s) + t)
        
        z = masked_x + (1 - self.mask) * (x * torch.exp(s) + t)
        
        log_det = s.sum(dim=1) # determinant is product of diagonal exp(s), log is sum s
        return z, log_det

    def inverse(self, z):
        masked_z = z * self.mask
        out = self.net(masked_z)
        s, t = out.chunk(2, dim=1)
        
        s = s * (1 - self.mask)
        t = t * (1 - self.mask)
        s = torch.tanh(s)
        
        # x = (z - t) * exp(-s)
        x = masked_z + (1 - self.mask) * ((z - t) * torch.exp(-s))
        
        log_det = -s.sum(dim=1)
        return x, log_det

if __name__ == "__main__":
    # Test
    model = RealNVP(10)
    x = torch.randn(5, 10)
    z, log_det = model(x)
    x_rec, inv_log_det = model.inverse(z)
    
    print(f"Reconstruction Error: {(x - x_rec).abs().max().item()}")
    print(f"Det Sum: {(log_det + inv_log_det).abs().max().item()}")
