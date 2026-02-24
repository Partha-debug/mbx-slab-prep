# mbx-slab-prep

A command-line tool that generates ready-to-run LAMMPS input for **MB-pol air/water slab simulations** using the [MBX](https://github.com/paesanilab/MBX) pair style.

Given a number of water molecules, a bulk density, and a vacuum multiplier, the tool packs the molecules with [Packmol](https://m3g.github.io/packmol/), assigns the correct MB-pol charges and masses, writes an `atom_style full` LAMMPS data file with the water slab centred in a periodic box, and generates a matching `in.equil.lammps` and `mbx.json` ready to run.

---

### MB-pol and MBX

**MB-pol** is a many-body potential energy surface for water that accurately reproduces bulk liquid, vapour, and ice properties from first principles.  It is implemented in the **MBX** library, which exposes a `pair_style mbx` plugin for LAMMPS.

```lammps
bond_style  none
angle_style none
```

### MB-pol parameters hard-coded in this tool

The following values are fixed by the MB-pol potential definition and are applied automatically.  They should not be changed.

| Quantity | Value | Notes |
|---|---|---|
| Atom type 1 | O | oxygen |
| Atom type 2 | H | hydrogen |
| O charge | −1.1128 e | MB-pol point charge |
| H charge | +0.5564 e | MB-pol point charge (molecule is neutral) |
| O mass | 15.9994 amu | |
| H mass | 1.00794 amu | |
| Real-space cutoff | 9.0 Å | standard for liquid water and slab geometries |
| 2-body cutoff | 9.0 Å | 6.5 Å for bulk; 9.0 Å used here for slab safety |
| 3-body cutoff | 7.0 Å | 4.5 Å for bulk; 7.0 Å used here for slab safety |
| Dipole method | CG | conjugate-gradient induced-dipole solver |
| Dipole tolerance | 1 × 10⁻⁸ | tight convergence for vSFG accuracy |
| Ewald α (elec) | 0.60 | PME splitting parameter |
| Ewald α (disp) | 0.60 | PME splitting parameter |

### Slab geometry

The generated box has the air/water interface in the **x–y plane**:

`Lx = Ly` is the natural cubic bulk side length at the requested density:

```
Lx = ( N × M_water / (ρ × N_A) )^(1/3)
```

`Lz = z_mult × Lx`.  `z_mult = 3` gives ≈ `Lx` Å of vacuum on each side, which is standard for vSFG spectroscopy simulations.

---

## This tool generates an *unequilibrated* initial structure

**You must run NVT equilibration with MB-pol before collecting any production data.**

The tool guarantees only that the **minimum O–O distance under full 3D periodic boundary conditions is ≥ 2.0 Å**.

A typical protocol after generating the files:

1. Energy minimisation (`minimize` in the generated input) — relaxes residual close contacts
2. NVT at target temperature for **at least 0.5 ns** before any production run

The generated `in.equil.lammps` handles both steps automatically.

---

## Prerequisites

| Requirement | Version | Notes |
|---|---|---|
| Python | ≥ 3.6 | |
| NumPy | ≥ 1.20 | `pip install -r requirements.txt` |
| Packmol | ≥ 20.x | `pip install -r requirements.txt` |

All dependencies including Packmol are installed via `pip install -r requirements.txt`.  The `packmol` pip package ships pre-compiled binaries for Linux, macOS, and Windows, so no manual compilation is needed.

---

## Installation

```bash
git clone https://github.com/Partha-debug/mbx-slab-prep.git
cd mbx-slab-prep
pip install -r requirements.txt
```
---

## Quick start

### Use config defaults (125 waters, 1.0 g/cm³, z_mult = 3, 298.15 K)

```bash
bash make_slab.sh
```

### Fully specified

```bash
bash make_slab.sh -n 256 -r 1.0 -z 3 -T 298.15 -d 0.2 -S 5000000 -o ./my_sim
```

### Override individual parameters; rest come from `config.conf`

```bash
bash make_slab.sh -n 512           # change system size only
bash make_slab.sh -z 4             # wider vacuum gap
bash make_slab.sh -T 350 -S 2000000  # higher temperature, shorter run
bash make_slab.sh -d 0.5           # larger timestep
bash make_slab.sh -s 99            # different Packmol seed
```

---

## Usage reference

```
bash make_slab.sh [flags]
```

| Flag | Type | Description | Default |
|------|------|-------------|---------|
| `-n` | INT | Number of water molecules | `NWATERS` in config |
| `-r` | FLOAT | Bulk density in g/cm³ (sets Lx = Ly = Lz_cubic) | `DENSITY` in config |
| `-z` | FLOAT | Z multiplier: Lz = z_mult × Lx | `Z_MULT` in config |
| `-T` | FLOAT | NVT temperature in K | `TEMP` in config |
| `-d` | FLOAT | MD timestep in fs | `TIMESTEP` in config |
| `-S` | INT | NVT equilibration run length in steps | `EQUIL_STEPS` in config |
| `-o` | DIR | Output directory | `OUTDIR` in config (`.`) |
| `-s` | INT | Packmol random seed | `SEED` in config (`42`) |
| `-t` | PATH | Path to water molecule PDB template | `water.pdb` next to script |
| `-P` | PATH | Python interpreter | active venv, then `python3` |
| `-c` | PATH | Config file | `config.conf` next to script |
| `-h` | | Print help and exit | |

---

## Configuration

`config.conf` sets defaults for every parameter.  Edit it once for your typical workflow and use flags only if you need single time overrides.


## Output files

Each run produces three files in the output directory (default - OUTDIR):

### `water_slab_<N>H2O_rho<ρ>_z<z>x.data`

LAMMPS `atom_style full` data file.

### `in.equil.lammps`

Ready-to-run LAMMPS equilibration input.

### `mbx.json`

MBX configuration file with all parameters documented inline.  Place this file in the same directory as `in.equil.lammps` before running LAMMPS.

---

## Note

If Packmol is not found after installation, make sure the virtual environment where you ran `pip install -r requirements.txt` is active before calling `make_slab.sh`, or set `PYTHON` in `config.conf` to point to that environment's interpreter.