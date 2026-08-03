# 🔬 Quantum Sensing

Educational material and Python simulations on quantum sensing based on nitrogen-vacancy (NV) centers in diamond.

![Python](https://img.shields.io/badge/Python-3.11-blue)
![License](https://img.shields.io/badge/License-MIT-green)
![Version](https://img.shields.io/badge/Version-v2.0-blue)

---

## Overview

This repository was developed during my extracurricular internship at the Nanomaterials and Nanotechnology Research Center (CINN, Spain).

It combines a complete introductory study manual on quantum sensing with an interactive Python simulator that reproduces Optically Detected Magnetic Resonance (ODMR) spectra of nitrogen-vacancy (NV) centers in diamond.

The project aims to provide an accessible introduction to NV-based quantum sensing by combining theoretical background, educational figures and numerical simulations.

---

# 📖 Quantum Sensing Study Manual

The repository includes a complete **75-page study manual** introducing the physical principles of quantum sensing based on nitrogen-vacancy (NV) centers in diamond.

The manual was written during my internship at CINN as educational material for students beginning in this research field.

Topics covered include:

- Fundamentals of quantum sensing
- Diamond as a quantum material
- Nitrogen-vacancy (NV) centers
- Spin Hamiltonian
- Zero-field splitting
- Zeeman interaction
- Electric-field and strain effects
- Spin relaxation and coherence
- Optically Detected Magnetic Resonance (ODMR)
- Magnetic sensitivity
- Noise and signal-to-noise ratio
- Microwave spin control
- RF signal generators
- Spectrum analyzers

The manual also contains:

- Original educational figures
- Technical diagrams
- Commented equations
- Figures adapted from scientific literature
- Complete bibliography
- Figure attribution table

📄 **Study Manual**

`docs/Quantum_Sensing_Study_Manual.pdf`

---

# 🐍 ODMR Python Simulator

The repository includes an interactive Python simulator that reproduces ODMR spectra from a single NV center in diamond.

The simulator numerically solves the effective spin Hamiltonian, calculates the resonance frequencies under different experimental conditions and estimates the longitudinal magnetic field from the Zeeman splitting.

All simulation results are automatically exported as publication-quality figures together with a structured simulation report.

---

## Current Features

- ✅ Interactive user input
- ✅ Complete spin Hamiltonian diagonalization
- ✅ Zero-field splitting simulation
- ✅ Zeeman interaction
- ✅ Transverse perturbation (strain/electric field parameter)
- ✅ Arbitrary magnetic-field orientation
- ✅ ODMR spectrum without magnetic field
- ✅ ODMR spectrum under an external magnetic field
- ✅ Resonance-frequency calculation
- ✅ Magnetic-field reconstruction from resonance frequencies
- ✅ Automatic generation of PNG figures
- ✅ Automatic generation of TXT reports

---

# Generated Figures

Each simulation automatically generates:

- ODMR spectrum without magnetic field
- ODMR spectrum with magnetic field
- Resonance frequencies versus magnetic field
- Resonance frequencies versus magnetic-field angle
- ODMR evolution for different magnetic fields
- Complete simulation report (.txt)

---

# Example Outputs

## ODMR Spectrum without Magnetic Field

![ODMR without magnetic field](images/odmr_without_magnetic_field.png)

---

## ODMR Spectrum with an External Magnetic Field

![ODMR with magnetic field](images/odmr_with_magnetic_field.png)

---

## Resonance Frequencies vs Magnetic Field

![Magnetic field](images/frequencies_vs_field.png)

---

## Resonance Frequencies vs Angle

![Angle](images/frequencies_vs_angle.png)

---

## ODMR Evolution for Different Magnetic Fields

![Evolution](images/evolution_odmr_with_magnetic_field.png)

---

# Repository Contents

```
Quantum-Sensing
│
├── docs/
│ └── Quantum_Sensing_Study_Manual.pdf
│
├── images/
│ ├── odmr_without_magnetic_field.png
│ ├── odmr_with_magnetic_field.png
│ ├── frequencies_vs_field.png
│ ├── frequencies_vs_angle.png
│ └── evolution_odmr_with_magnetic_field.png
│
├── results_odmr/
│ └── results_odmr.txt
│
├── odmr_simulation.py
│
└── README.md
```

---

# Physics Implemented

The current simulator includes the following physical concepts:

- Zero-field splitting (D = 2.87 GHz)
- Effective spin Hamiltonian (S = 1)
- Zeeman interaction
- Arbitrary magnetic-field orientation
- Transverse perturbation parameter (E)
- Hamiltonian diagonalization
- Lorentzian ODMR spectra
- Resonance-frequency calculation
- Longitudinal magnetic-field estimation

Future versions will progressively incorporate additional concepts presented in the study manual.

---

# Research Progress

Current project status:

- ✅ Literature review completed
- ✅ 75-page study manual completed
- ✅ Original educational figures completed
- ✅ ODMR simulator Version 2.0
- ✅ Scientific figure generation
- ✅ Automatic simulation report
- 🔄 Experimental noise simulation
- 🔄 Magnetic sensitivity estimation
- 🔄 Lorentzian fitting
- 🔄 Pulsed ODMR simulations

---

# Requirements

Python 3.11 or later

Libraries:

- NumPy
- Matplotlib

Install with:

```bash
pip install numpy matplotlib
```

---

# Running the Simulator

Clone the repository:

```bash
git clone https://github.com/sofianunezda/Quantum-sensing.git
```

Run:

```bash
python odmr_simulation.py
```

The program will ask for:

- Magnetic field
- Magnetic-field angle
- Transverse perturbation

and automatically generate all figures and the simulation report.

---

# References

The project is based on scientific literature used during my internship at CINN.

Main references include:

- Barry, J. F. et al. (2020)
- Degen, C. L., Reinhard, F. & Cappellaro, P. (2017)
- Doherty, M. W. et al. (2013)
- Dréau, A. et al. (2011)
- Rondin, L. et al. (2014)

Additional references are included in the study manual.

---

# Project Goals

The purpose of this repository is to provide an educational introduction to quantum sensing based on nitrogen-vacancy (NV) centers in diamond.

The project combines:

- theoretical background;
- educational material;
- scientific illustrations;
- numerical simulations.

Future developments will progressively extend the simulator while maintaining its educational focus and its correspondence with the accompanying study manual.

---

# License

This project is distributed under the MIT License.

If you use the study manual, figures or simulation code for academic or educational purposes, citation of this repository is appreciated.

---

# Author

**Sofía Núñez de Andrés**

Physics Student

University of Oviedo

Extracurricular Internship

Nanomaterials and Nanotechnology Research Center (CINN)

GitHub:

https://github.com/sofianunezda
