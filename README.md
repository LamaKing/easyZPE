![logo](logo.png)

# easyZPE

Zero Point Energy (and finite-temperature vibrational free energy) from frozen phonons using [Quantum ESPRESSO](https://www.quantum-espresso.org/) `pw.x` output.

The approach is deliberately simple: no symmetry reduction, no DFPT, just central finite differences on a user-selected set of atoms, a full 3N×3N dynamical matrix at Γ, and `numpy` for diagonalisation. This is appropriate when you only care about a small adsorbate on a frozen surface and want results without the overhead of a full phonon calculation.

## Method

1. **Generate displaced configurations** (`generate_conf.py`): for each selected atom and each Cartesian direction, write a QE `pw.x` input with the atom displaced by ±δ.
2. **Run QE**: standard `pw.x` SCF with `tprnfor = .true.` on all displaced configurations.
3. **Compute frequencies** (`compute_freq.py`): read forces from the QE XML outputs, build the dynamical matrix via central finite differences, diagonalise it, and compute ZPE and thermal corrections.

The dynamical matrix element is:

$$D_{ij} = -\frac{F_i(+\delta_j) - F_i(-\delta_j)}{2\,\delta_j\,\sqrt{m_i\,m_j}}$$

where $i,j$ run over (atom, Cartesian) pairs of the selected atoms and $m$ is the atomic mass. ZPE and vibrational free energy follow from the harmonic approximation:

$$\text{ZPE} = \sum_k \frac{\hbar\omega_k}{2}$$

$$E_\text{vib} = \sum_k \hbar\omega_k \left[\frac{1}{2} + \frac{1}{e^{\hbar\omega_k/k_BT}-1}\right]$$

$$A_\text{vib}(T) = \sum_k \frac{\hbar\omega_k}{2} + k_BT\ln\!\left(1 - e^{-\hbar\omega_k/k_BT}\right)$$

Modes with $\omega < 0$ (soft/translational modes due to DFT noise) are skipped in all sums.

## Dependencies

| Package | Purpose |
|---|---|
| [ASE](https://wiki.fysik.dtu.dk/ase/) | Read QE input files, geometry manipulation |
| [qeschema](https://github.com/QEF/qeschema) | Parse QE XML output |
| [numpy](https://numpy.org/) | Linear algebra |
| [load_calculator](https://github.com/andreasilva759/QE-utilities) | Parse QE input namelists (copy included) |

Install Python dependencies:
```bash
pip install ase qeschema numpy
```

## Usage

### Step 1 — generate displaced inputs

```bash
python generate_conf.py INDICES [options]
```

**Positional arguments:**

| Argument | Description |
|---|---|
| `INDICES` | Space-separated 0-based indices of the atoms to displace |

**Options:**

| Flag | Default | Description |
|---|---|---|
| `-f`, `--filename` | `pw.in` | QE input file with the equilibrium geometry |
| `-d`, `--displ` | `0.01` | Displacement magnitude [Å] |
| `--outdir` | `./` | Root directory for the displacement folders |
| `--outfile` | `pw.in` | Filename for each displaced input |
| `--extra_cp` | — | Extra files to copy into each folder (e.g. job submission script) |
| `--debug` | — | Verbose output |

This creates one folder per displacement: `dpos_<atom>_<coord>/` and `dneg_<atom>_<coord>/`, each containing a ready-to-run QE input.

**Example** — displace atoms 0, 1, 2 (an adsorbed H₂O molecule) with 0.02 Å:
```bash
python generate_conf.py 0 1 2 -f pw.in -d 0.02 --outdir vib/ --extra_cp job.sh
```

### Step 2 — run QE

Run `pw.x` in each displacement folder. The XML output file (`pwscf.xml` by default) must be present in every folder before proceeding.

### Step 3 — compute frequencies and ZPE

```bash
python compute_freq.py INDICES [options]
```

**Options:**

| Flag | Default | Description |
|---|---|---|
| `-f`, `--filename` | `pw.in` | Equilibrium QE input file (for masses and geometry check) |
| `-d`, `--displ` | `0.01` | Displacement used in step 1 [Å] — must match |
| `--indir` | `./` | Directory containing the `dpos_*/dneg_*` folders |
| `--outdir` | `./` | Directory for output files |
| `--qeout` | `pwscf.xml` | Name of the QE XML output (same for all runs) |
| `-T`, `--temperature` | `298.15` | Temperature for thermal corrections [K] |
| `--debug` | — | Verbose output |

**Example** — continuing from the example above:
```bash
python compute_freq.py 0 1 2 -f pw.in -d 0.02 --indir vib/ --outdir vib/results/
```

**Output files** (written to `--outdir`):

| File | Content |
|---|---|
| `eigenvalues.dat` | Eigenvalues, frequencies (THz, cm⁻¹), ZPE, TS, A_vib per mode |
| `eigenvectors.dat` | Eigenvectors (columns) of the dynamical matrix |
| `force_posdx.dat` | Force matrix for positive displacements [Eh/a₀] |
| `force_negdx.dat` | Force matrix for negative displacements [Eh/a₀] |
| `force_const.dat` | Force constant matrix [Eh/a₀²] |
| `dyna_matr.dat` | Mass-weighted dynamical matrix [Eh/(a₀² Da)] |

## Test cases

Three reference calculations are included in `test_case-*/`. Run:

```bash
python test_easyZPE.py
```

| System | Atoms displaced | Expected ZPE |
|---|---|---|
| H₂ molecule | 0, 1 | ~0.268 eV |
| H₂O molecule | 0, 1, 2 | ~0.566 eV |
| H₂O on Ni surface | 56, 57, 58 | ~0.60–0.65 eV |

The test for `generate_conf` re-generates all displaced geometries for H₂O and checks positions against the committed reference inputs.

## Notes and caveats

- **No symmetry**: all 3N degrees of freedom are displaced and computed independently. For a single adsorbate this is fine; for many atoms it scales poorly.
- **Surface atoms are frozen**: only the atoms you list are displaced. The surrounding surface is kept fixed, which is the usual frozen-surface approximation.
- **Near-zero modes**: translations and frustrated translations of the adsorbate will appear as near-zero or slightly imaginary modes due to DFT residual forces. They are skipped in ZPE/Avib sums automatically and flagged in the output.
- **Displacement size**: 0.01 Å is a reasonable default. If you see large imaginary frequencies on real modes, tighter SCF convergence (`conv_thr`) helps more than changing δ.
- **Thermal quantities**: `TS` and `A_vib` can be unreliable for very soft modes (< ~50 cm⁻¹) that are not true translations — cell-size artefacts produce artificially large entropy contributions. ZPE is unaffected.

