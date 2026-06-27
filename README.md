[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.20771213-blue)](https://doi.org/10.5281/zenodo.20771213)

## HYDRA (Hydrothermal Dynamics of Replicative Amplification)
Version 2.1: Comprehensive Multiphysics Framework for Prebiotic Polymer Selection

https://doi.org/10.5281/zenodo.18393990

https://doi.org/10.5281/zenodo.20771213

MIT License

## Overview
HYDRA is a multiphysics simulation framework designed to model the chemical evolution and selection of prebiotic polymers (RNA, DNA, peptides) within hydrothermal pore networks. Built on experimental kinetics and 3D hydrodynamic transport, HYDRA enables researchers to test origins-of-life scenarios with high physical and chemical realism.

## Multi-Module Architecture

The framework is structured into six integrated modules:
1.  Chemical Kinetics Engine: Incorporates metal catalysis (Fe²⁺), base hydrolysis, and Michaelis-Menten polymerization.
2.  Pore Network Generator: Generates synthetic or tomography-derived 3D environments.
3.  Transport Solver: Models diffusion (Faxén-corrected), thermophoresis (Soret effect), and electrostatics.
4.  Mineral Surface Chemistry: Site-specific adsorption and surface-mediated catalysis.
5.  Environmental Field Generator: Resolves temperature gradients, pH fluctuations, and metal speciation.
6.  Experimental Validator: Built-in metrics for clinical and laboratory data comparison.

Key Features in v2.1
- Faxén Diffusion Correction: Accurate modeling of molecular movement in micro-confined spaces.
- Enhanced Metal Catalysis: Dynamic Haber-Weiss kinetics for Fe²⁺-mediated degradation.
- Optimized Performance: Numba-accelerated kernels for high-speed simulations.
- Zenodo Integration: Ready for academic citation.

## Citation
If you use HYDRA in your research, please cite it as:
Hashemi, R. (2026).HYDRA: A Comprehensive Multiphysics Framework for Simulating Prebiotic Polymer Selection in Hydrothermal Environments. https://doi.org/10.5281/zenodo.18393990

https://doi.org/10.5281/zenodo.20771213

## Author

## Seyed Mohammad Reza Hashemi(Reza Hashemi)  

Former Member of the Pasteur Institute of Iran 

Email: mrhashemi2000@gmail.com

ORCID : 0009-0002-0645-5180

This project is licensed under the MIT License.

## Versions:
*Main Branch (v2.1): The current high-fidelity multiphysics framework. Focuses on rigorous chemical kinetics, thermophoresis, and pore-geometry effects. (Recommended for research).

*Tag 1.0.0: Initial version focused on visualization and basic polymer decay/selection metrics.

