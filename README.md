# Ubiquitin Charge Masking

This project utilizes the **Wako-Saitô-Muñoz-Eaton (WSME)** statistical mechanical model combined with **Monte Carlo (MC)** simulations to investigate how masking certain charges affects the protein stability of Ubiquitin (**1UBQ**).

*Based on original code by Prof. Athi N. Naganathan - FesCalc_Block.m and cmapCalcElecBlock.m*

---

## 1. Model Calibration (`calc_H_avg.m`)

Before simulating mutants, the intrinsic energy parameter (ene) for the WSME model is calibrated to match experimental thermal stability data.

### Calibration Protocol

* **Parameter Sweep:** We generate Heat Capacity (Cp vs T) profiles across a range of 250K to 420K for ene values between -0.085 and -0.125 J/mol.
* **Entropy Weighting:** Uses STRIDE-derived secondary structure data to apply specific entropic penalties:
* **Proline:** 1.0 (no entropic cost).
* **Glycine/Disordered Loops:** Enhanced penalty (z_vec).

> **Results:**
> Optimal ene = -0.091 J/mol |
> Predicted Melting Temperature = 357K |
> Experimental Melting Temperature = 358K |
> Melting Temperature Error **0.279%**.

---

## 2. Identification of Surface Charged Residues (SCRs)

We identified residues for masking by calculating **Relative Solvent Accessibility (RSA)** by dividing **Solvent Accessible Surface Area (SASA)** values from STRIDE by **Maximum allowed solvent accessibilites** _see sources_. Residues with RSA>0.3 and non-zero charge at pH=7 were selected (Histidine was considered uncharged).

**Target Mask List (Residue Indices):**
`[6, 11, 16, 18, 24, 32, 33, 34, 39, 42, 48, 51, 52, 54, 58, 63, 64, 72, 74]`

---

## 3. Mutant Free Energy Profiling (`Generate_Mutants.m`)

The energetic impact of masking each identified SCR is calculated at physiological temperature (310K). The script utilizes the `get_Free_Energy` function to zero out the charge contribution in the contact distance matrix for the selected residue index.

_Line 5 of each Energy Profile Text file was manually commented out for MATLAB_

---

## 4. Monte Carlo Simulation (`Monte_Carlo.m`)

We simulate stochastic folding dynamics using a dual-walker approach with a **Metropolis criterion**:

1. **WT Walker (Control):** Static wild-type energy landscape.
2. **Dynamic Walker (Experiment):** Every  steps (a switch interval ), the energy landscape is swapped for a randomly selected mutant profile.

---

## Project Structure

| Directory | Description |
| --- | --- |
| `Code_Base/` | All `.m` scripts (Calibration, Generation, MC Simulation, and Plotting). |
| `Input_Files/` | Prerequisite data files (Contact maps, distances, STRIDE output). |
| `Free_Energies/` | Saved energy profiles for Wild Type and all Mutants. |
| `Output_Files/` | Heat Capacity profiles, trajectories, and a `Plots/` subfolder. |

---
