"""
HYDRA (Hydrothermal Dynamics of Replicative Amplification)
Version 2.1: Enhanced Multiphysics Framework for Prebiotic Polymer Selection

Author: Seyed Mohammad Reza Hashemi (Reza Hashemi)  Intelligence-Augmented (IA)
Affiliation: Former Member of the Pasteur Institute of Iran
Email: mrhashemi2000@gmail.com
DOI: 10.5281/zenodo.20771213
License: MIT

Description:
A comprehensive multiphysics framework for simulating polymer selection in 
hydrothermal pore networks. Integrates chemical kinetics, pore hydrodynamics, 
mineral surface chemistry, and temperature-gradient transport.
"""

import numpy as np
from numba import njit, prange
import scipy.ndimage as ndi
import h5py
from dataclasses import dataclass, field
from datetime import datetime
from enum import IntEnum
from pathlib import Path
import json

# ======================================================
# 0. PHYSICAL CONSTANTS & ENUMS
# ======================================================
class PhysicalConstants:
    kB = 1.380649e-23      # Boltzmann constant [J/K]
    NA = 6.02214076e23     # Avogadro number
    R = 8.314462618        # Gas constant [J/mol.K]
    VISCOSITY_W = 0.00089  # Water viscosity at 25C [Pa.s]
    T_REF = 298.15         # Reference temperature [K]

class PolymerType(IntEnum):
    RNA = 0
    DNA = 1
    PEPTIDE = 2
    NUCLEOTIDE = 3

class MineralType(IntEnum):
    MONTMORILLONITE = 0
    PYRITE = 1
    ZIRCON = 2
    QUARTZ = 3

# ======================================================
# 1. CHEMICAL KINETICS ENGINE (Module 1)
# ======================================================
@njit
def arrhenius_ph_rate(k25, Ea, T, pH):
    """Calculates hydrolysis rate with temperature and pH correction."""
    k_t = k25 * np.exp(-Ea / PhysicalConstants.R * (1.0 / T - 1.0 / 298.15))
    ph_factor = 10.0 ** (1.2 * np.abs(pH - 7.0)) # pH dependence from paper
    return k_t * ph_factor

class KineticsEngine:
    def __init__(self):
        # Data from Table 1 of the HYDRA paper
        self.params = {
            PolymerType.RNA: {'k25': 2.2e-9, 'Ea': 121000.0, 'kcat_fe': 3.5e4, 'Kassoc_fe': 1.2e3},
            PolymerType.DNA: {'k25': 7.33e-16, 'Ea': 134000.0, 'kcat_fe': 8.7e3, 'Kassoc_fe': 2.4e2}
        }

    def get_degradation_rate(self, p_type, temp, pH, fe2_conc):
        if p_type not in self.params: return 1e-12
        p = self.params[p_type]
        base_rate = arrhenius_ph_rate(p['k25'], p['Ea'], temp, pH)
        
        # Metal Catalysis (Modified Haber-Weiss)
        assoc = p['Kassoc_fe'] * fe2_conc
        enhancement = 1.0 + (p['kcat_fe'] * assoc / (1.0 + assoc))
        return base_rate * enhancement

# ======================================================
# 2. TRANSPORT PHYSICS (Module 3)
# ======================================================
class TransportSolver:
    @staticmethod
    @njit
    def faxen_diffusion(r_poly, r_pore, T):
        """Module 3: Hindered diffusion in confinement with Faxén correction."""
        D0 = (PhysicalConstants.kB * T) / (6.0 * np.pi * PhysicalConstants.VISCOSITY_W * r_poly)
        lambda_ratio = r_poly / r_pore
        if lambda_ratio >= 1.0: return 1e-22
        # Faxen formula from paper
        correction = (1.0 - lambda_ratio)**1.5 * (1.0 - 2.104*lambda_ratio + 2.089*lambda_ratio**3 - 0.948*lambda_ratio**5)
        return D0 * correction

    @staticmethod
    @njit
    def thermophoresis_drift(D, T_grad, T, S_t=0.01):
        """v_T = -D * S_T * grad(T) / T"""
        return -D * S_t * T_grad / T

# ======================================================
# 3. PORE NETWORK GENERATOR (Module 2)
# ======================================================
class PoreNetwork:
    def __init__(self, shape=(50, 50, 100), res=1e-6, porosity=0.35):
        self.shape = shape
        self.res = res
        self.grid = self._generate_grf(porosity)
        # Distance transform for local pore radii
        self.dist_map = ndi.distance_transform_edt(self.grid) * res

    def _generate_grf(self, porosity):
        noise = np.random.randn(*self.shape)
        smoothed = ndi.gaussian_filter(noise, sigma=1.5)
        threshold = np.percentile(smoothed, 100 * (1.0 - porosity))
        return smoothed > threshold

# ======================================================
# 4. CORE SIMULATION ENGINE (HYDRA)
# ======================================================
class HYDRA:
    def __init__(self, config=None):
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        print(f"[{self.timestamp}] Initializing HYDRA v2.1...")
        
        self.net = PoreNetwork()
        self.kinetics = KineticsEngine()
        self.transport = TransportSolver()
        
        # Environmental Fields (Module 5)
        self.temp_grad = 20.0e3 # 20 K/mm gradient from paper
        self.T_bottom = 353.15 # 80C
        self.pH = 7.3
        self.fe2_conc = 0.00042 # 0.42 mM optimum
        
        self.molecules = []
        self.stats = {'time': [], 'dna': [], 'rna': []}

    def seed_population(self, rna=1000, dna=50):
        pore_indices = np.argwhere(self.net.grid)
        for i in range(rna + dna):
            m_type = PolymerType.RNA if i < rna else PolymerType.DNA
            idx = pore_indices[np.random.randint(len(pore_indices))]
            self.molecules.append({
                'type': m_type,
                'pos': idx.astype(float) * self.net.res,
                'damage': 0.0,
                'length': 20
            })
        print(f"Successfully seeded {len(self.molecules)} molecules.")

    def step(self, dt=3600.0):
        """Simulation step (1 hour default)"""
        new_mols = []
        for mol in self.molecules:
            # Get local environment
            grid_idx = (mol['pos'] / self.net.res).astype(int)
            grid_idx = np.clip(grid_idx, 0, np.array(self.net.shape)-1)
            
            # Local Temp (Z-gradient)
            z_m = mol['pos'][2]
            local_T = self.T_bottom - (self.temp_grad * z_m)
            local_pore_r = max(1e-9, self.net.dist_map[tuple(grid_idx)])
            
            # 1. Kinetics: Degradation
            k_deg = self.kinetics.get_degradation_rate(mol['type'], local_T, self.pH, self.fe2_conc)
            mol['damage'] += k_deg * dt
            
            if mol['damage'] >= 1.0: continue # Degraded
            
            # 2. Transport: Diffusion + Thermophoresis
            D = self.transport.faxen_diffusion(1e-9, local_pore_r, local_T)
            v_t = self.transport.thermophoresis_drift(D, self.temp_grad, local_T)
            
            # Random walk + Drift
            noise = np.random.normal(0, np.sqrt(2 * D * dt), 3)
            mol['pos'] += noise
            mol['pos'][2] += v_t * dt # Z-drift
            
            # Boundary conditions (reflective)
            limit = np.array(self.net.shape) * self.net.res
            mol['pos'] = np.clip(mol['pos'], 0, limit)
            
            new_mols.append(mol)
            
        self.molecules = new_mols

    def run(self, hours=168):
        """Runs simulation for the selection window (e.g. 7 days)"""
        print(f"Running for {hours} hours...")
        for h in range(hours):
            self.step(dt=3600.0)
            if h % 24 == 0:
                dna_c = sum(1 for m in self.molecules if m['type'] == PolymerType.DNA)
                rna_c = sum(1 for m in self.molecules if m['type'] == PolymerType.RNA)
                print(f"Day {h//24}: DNA={dna_c}, RNA={rna_c}")
                self.stats['dna'].append(dna_c)
                self.stats['rna'].append(rna_c)

    def export(self):
        filename = f"HYDRA_V2_1_{self.timestamp}.h5"
        with h5py.File(filename, 'w') as f:
            f.attrs['Author'] = "Reza Hashemi"
            f.attrs['DOI'] = "10.5281/zenodo.18393990"
            f.create_dataset('pore_network', data=self.net.grid)
            f.create_dataset('stats_dna', data=np.array(self.stats['dna']))
            f.create_dataset('stats_rna', data=np.array(self.stats['rna']))
        print(f"Data exported to {filename}")

# ======================================================
# EXECUTION
# ======================================================
if __name__ == "__main__":
    model = HYDRA()
    model.seed_population(rna=5000, dna=100) # Trace DNA scenario
    model.run(hours=168) # 1 week selection
    model.export()
    print("HYDRA simulation completed. Results ready for analysis.")
