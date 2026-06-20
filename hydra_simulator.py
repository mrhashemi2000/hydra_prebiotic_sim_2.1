"""
HYDRA (Hydrothermal Dynamics of Replicative Amplification)
Version 2.1: Enhanced Implementation based on HYDRA300.pdf (Zenodo DOI: 10.5281/zenodo.18393990)

Comprehensive Multiphysics Framework for Prebiotic Polymer Selection in Hydrothermal Environments.

Author: سید محمدرضا هاشمی (Reza Hashemi)
Affiliation: Former Member of the Pasteur Institute of Iran
License: MIT
"""

import numpy as np
from numba import njit
import scipy.ndimage as ndi
import h5py
from dataclasses import dataclass
from datetime import datetime
from enum import Enum, auto
import matplotlib.pyplot as plt

# ======================================================
# PHYSICAL CONSTANTS & CORE ENUMS
# ======================================================
class PhysicalConstants:
    NA = 6.02214076e23
    kB = 1.380649e-23
    R = 8.3144626
    WATER_VISCOSITY_25C = 0.000890
    WATER_DIELECTRIC = 78.4


class PolymerType(Enum):
    RNA = auto()
    DNA = auto()
    PEPTIDE = auto()
    NUCLEOTIDE = auto()


# ======================================================
# DATA STRUCTURES
# ======================================================
@dataclass
class EnvironmentalConditions:
    temperature: float
    pH: float = 7.0
    water_activity: float = 1.0
    ionic_strength: float = 0.1
    fe2_conc: float = 0.0


# ======================================================
# KINETICS ENGINE (Module 1)
# ======================================================
class KineticsEngine:
    def __init__(self):
        self.rates = {
            'RNA': {'k25': 2.2e-9, 'Ea': 121000.0},
            'DNA': {'k25': 7.3e-16, 'Ea': 134000.0}
        }

    @staticmethod
    @njit
    def calculate_rate(k25, Ea, T, pH, R=8.3144626):
        k_t = k25 * np.exp(-Ea / R * (1/T - 1/298.15))
        ph_factor = 10 ** (max(0, pH - 7.0))
        return k_t * ph_factor

    def get_degradation_rate(self, poly_type: str, cond: EnvironmentalConditions) -> float:
        if poly_type not in self.rates:
            return 0.0
        base_rate = self.calculate_rate(
            self.rates[poly_type]['k25'],
            self.rates[poly_type]['Ea'],
            cond.temperature,
            cond.pH
        )
        metal_factor = 1.0
        if cond.fe2_conc > 0:
            boost = 52.0 if poly_type == 'RNA' else 6.0
            metal_factor = 1.0 + (boost * (cond.fe2_conc / 0.0005))
        return base_rate * metal_factor


# ======================================================
# TRANSPORT PHYSICS (Module 3 - Enhanced with Thermophoresis)
# ======================================================
class TransportPhysics:
    @staticmethod
    def calculate_diffusion_faxen(r_poly, r_pore, T):
        d0 = (PhysicalConstants.kB * T) / (6 * np.pi * PhysicalConstants.WATER_VISCOSITY_25C * r_poly)
        lambda_ratio = r_poly / r_pore
        if lambda_ratio >= 1.0:
            return 1e-20
        correction = (1 - lambda_ratio)**1.5 * (1 - 2.104*lambda_ratio + 2.089*lambda_ratio**3 - 0.948*lambda_ratio**5)
        return d0 * correction

    @staticmethod
    def thermophoresis_velocity(temp, grad_T, S_T=0.01):
        """Simplified thermophoresis (Soret coefficient)."""
        return -S_T * grad_T * temp


# ======================================================
# PORE NETWORK (Module 2)
# ======================================================
class PoreNetwork:
    def __init__(self, shape=(64, 64, 64), resolution=1e-6, porosity=0.35):
        self.shape = shape
        self.res = resolution
        self.grid = self._generate_synthetic_pores(porosity)
        
    def _generate_synthetic_pores(self, porosity):
        noise = np.random.randn(*self.shape)
        smoothed = ndi.gaussian_filter(noise, sigma=2.0)
        threshold = np.percentile(smoothed, 100 * (1 - porosity))
        return smoothed > threshold


# ======================================================
# MAIN SIMULATOR
# ======================================================
class HYDRASimulator:
    def __init__(self):
        print(f"[{datetime.now().strftime('%H:%M:%S')}] Initializing HYDRA 2.1 Kernel (based on Zenodo 10.5281/zenodo.18393990)...")
        self.pore_net = PoreNetwork()
        self.kinetics = KineticsEngine()
        self.molecules = []
        self.time = 0.0
        self.temp_field = np.linspace(300, 370, self.pore_net.shape[2])
        self.fe2_conc = 0.0005
        self.history = {'time': [], 'dna': [], 'rna': []}

    def seed_molecules(self, rna_count=1000, dna_count=50):
        pore_indices = np.argwhere(self.pore_net.grid)
        for i in range(rna_count + dna_count):
            m_type = PolymerType.RNA if i < rna_count else PolymerType.DNA
            idx = pore_indices[np.random.randint(len(pore_indices))]
            self.molecules.append({
                'type': m_type,
                'pos': idx.astype(float) * self.pore_net.res,
                'damage': 0.0,
                'length': 20
            })

    def run_step(self, dt=3600):
        new_molecules = []
        for mol in self.molecules:
            z_idx = int(mol['pos'][2] / self.pore_net.res)
            z_idx = np.clip(z_idx, 0, self.pore_net.shape[2]-1)
            temp = self.temp_field[z_idx]
            cond = EnvironmentalConditions(temperature=temp, fe2_conc=self.fe2_conc)
            
            # Degradation
            rate = self.kinetics.get_degradation_rate(mol['type'].name, cond)
            mol['damage'] += rate * dt
            
            # Transport
            d_coeff = TransportPhysics.calculate_diffusion_faxen(1e-9, 5e-6, temp)
            drift = np.random.normal(0, np.sqrt(2 * d_coeff * dt), 3)
            mol['pos'] += drift
            
            if mol['damage'] < 1.0:
                new_molecules.append(mol)
        
        self.molecules = new_molecules
        self.time += dt
        
        # Record history
        dna = sum(1 for m in self.molecules if m['type'] == PolymerType.DNA)
        rna = sum(1 for m in self.molecules if m['type'] == PolymerType.RNA)
        self.history['time'].append(self.time)
        self.history['dna'].append(dna)
        self.history['rna'].append(rna)

    def get_stats(self):
        dna = sum(1 for m in self.molecules if m['type'] == PolymerType.DNA)
        rna = sum(1 for m in self.molecules if m['type'] == PolymerType.RNA)
        return {'DNA': dna, 'RNA': rna, 'Total': len(self.molecules), 'Time': self.time}

    def plot_results(self):
        plt.figure(figsize=(10, 6))
        plt.plot(self.history['time'], self.history['dna'], 'b-', label='DNA')
        plt.plot(self.history['time'], self.history['rna'], 'r-', label='RNA')
        plt.xlabel('Time (hours)')
        plt.ylabel('Surviving Polymers')
        plt.title('HYDRA 2.1 - DNA vs RNA Selection in Hydrothermal Pores')
        plt.legend()
        plt.grid(True)
        plt.savefig('hydra_results.png')
        plt.show()

    def save_to_h5(self, filename="hydra_output.h5"):
        with h5py.File(filename, 'w') as f:
            f.attrs['author'] = "Reza Hashemi"
            f.attrs['version'] = "2.1"
            f.attrs['doi'] = "10.5281/zenodo.18393990"
            f.create_dataset('pore_grid', data=self.pore_net.grid)
        print(f"✅ Results saved to {filename}")


def main():
    sim = HYDRASimulator()
    sim.seed_molecules(rna_count=1000, dna_count=50)
    
    print("🚀 Starting HYDRA 2.1 Simulation...")
    for step in range(100):
        sim.run_step(dt=3600)
        if step % 10 == 0 or step == 99:
            stats = sim.get_stats()
            print(f"Step {step:3d}: DNA={stats['DNA']:3d} | RNA={stats['RNA']:3d} | Total={stats['Total']}")
    
    sim.save_to_h5()
    sim.plot_results()
    print("🎉 HYDRA 2.1 simulation completed successfully. Results plotted.")


if __name__ == "__main__":
    main()
