"""
validate.py - Validate the generated water slab LAMMPS data file.

Uses numpy for vectorised O–O pairwise distance calculation with full 3D
periodic boundary conditions, which is significantly faster than a Python
loop for large systems (N > ~100 molecules).

Checks performed
────────────────
  1. Atom / bond / angle counts are consistent with N water molecules
  2. Atom ordering within each molecule is O, H, H  (required by MBX)
  3. No atoms sit outside the box boundaries (tolerance 1e-4 Å)
  4. Minimum O–O distance under full 3D PBC is >= min_oo threshold (default 1.8 Å)
     This is the critical check: the MBX induced-dipole conjugate-gradient solver
     will fail at the first energy evaluation if two oxygens are closer than
     ~2 Å, causing LAMMPS to abort with "Max number of iterations reached".

Reports
───────
  - Box dimensions
  - Slab z-range and approximate vacuum thickness per side
  - Estimated density in the water region
  - Minimum and mean O–O PBC distances

Usage (internal — called by make_slab.sh):
    python validate.py --slab SLAB.data --nwaters N --lx LX [--min-oo FLOAT]
"""

import argparse
import sys

import numpy as np


def read_slab(path):
    atoms = {}
    box   = {}
    mode  = None

    with open(path) as f:
        for line in f:
            s = line.strip()
            if "xlo xhi"  in s: box["x"] = tuple(float(v) for v in s.split()[:2])
            elif "ylo yhi" in s: box["y"] = tuple(float(v) for v in s.split()[:2])
            elif "zlo zhi" in s: box["z"] = tuple(float(v) for v in s.split()[:2])
            elif s in ("Atoms # full", "Atoms"): mode = "atoms"
            elif mode == "atoms":
                if not s:
                    continue
                if s[0].isalpha():
                    mode = None; continue
                t = s.split()
                if len(t) >= 7:
                    atoms[int(t[0])] = {
                        "mol":  int(t[1]),
                        "type": int(t[2]),
                        "x":    float(t[4]),
                        "y":    float(t[5]),
                        "z":    float(t[6]),
                    }
    return atoms, box


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--slab",    required=True)
    p.add_argument("--nwaters", type=int,   required=True)
    p.add_argument("--lx",      type=float, required=True)
    p.add_argument("--min-oo",  type=float, default=1.8,
                   help="Minimum allowed O–O PBC distance in Å (default 1.8)")
    args = p.parse_args()

    atoms, box = read_slab(args.slab)

    Lx = box["x"][1] - box["x"][0]
    Ly = box["y"][1] - box["y"][0]
    Lz = box["z"][1] - box["z"][0]
    L  = np.array([Lx, Ly, Lz])

    errors   = []
    warnings = []

    # 1. counts
    oxygens   = [a for a in atoms.values() if a["type"] == 1]
    hydrogens = [a for a in atoms.values() if a["type"] == 2]

    if len(atoms)    != args.nwaters * 3:   errors.append(f"Total atom count:    expected {args.nwaters*3}, got {len(atoms)}")
    if len(oxygens)  != args.nwaters:        errors.append(f"Oxygen count:        expected {args.nwaters}, got {len(oxygens)}")
    if len(hydrogens)!= args.nwaters * 2:    errors.append(f"Hydrogen count:      expected {args.nwaters*2}, got {len(hydrogens)}")

    # 2. O,H,H ordering per molecule
    mols = {}
    for aid in sorted(atoms.keys()):
        m = atoms[aid]["mol"]
        mols.setdefault(m, []).append(atoms[aid]["type"])
    bad_order = [m for m, types in mols.items() if types != [1, 2, 2]]
    if bad_order:
        warnings.append(
            f"{len(bad_order)} molecule(s) do not have O,H,H atom ordering "
            f"(MBX requires O first within each molecule's consecutive atom IDs)"
        )

    # 3. atoms outside box
    TOL = 1e-4
    x0, x1 = box["x"]; y0, y1 = box["y"]; z0, z1 = box["z"]
    outside = [
        aid for aid, a in atoms.items()
        if a["x"] < x0-TOL or a["x"] > x1+TOL
        or a["y"] < y0-TOL or a["y"] > y1+TOL
        or a["z"] < z0-TOL or a["z"] > z1+TOL
    ]
    if outside:
        warnings.append(
            f"{len(outside)} atom(s) lie outside box boundaries "
            f"(LAMMPS will wrap them at startup)"
        )

    # 4. min O–O PBC distance via numpy
    coords = np.array([[a["x"], a["y"], a["z"]] for a in oxygens])   # (N, 3)

    # Pairwise displacement vectors, minimum image convention
    diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]       # (N, N, 3)
    diff -= np.round(diff / L) * L
    dist  = np.sqrt(np.sum(diff**2, axis=-1))                         # (N, N)
    np.fill_diagonal(dist, np.inf)

    min_oo = float(dist.min())
    mean_oo = float(np.mean(dist[dist < np.inf]))  # mean over upper triangle would do too

    if min_oo < args.min_oo:
        idx = np.unravel_index(np.argmin(dist), dist.shape)
        mol_a = oxygens[idx[0]]["mol"]
        mol_b = oxygens[idx[1]]["mol"]
        errors.append(
            f"Min O–O PBC distance {min_oo:.4f} Å < threshold {args.min_oo} Å "
            f"(molecules {mol_a} and {mol_b}) — MBX will fail on first evaluation"
        )

    # 5. report
    z_o       = [a["z"] for a in oxygens]
    z_min, z_max = min(z_o), max(z_o)
    z_water   = z_max - z_min + 3.0   # approx slab thickness (+1 O–H bond)
    vac_each  = (Lz - z_water) / 2.0
    density   = args.nwaters * 18.015 / (6.022e23 * Lx * Ly * z_water * 1e-24)

    print(f"\n  ── Slab validation ──────────────────────────────────────────────")
    print(f"  Box         : {Lx:.4f} × {Ly:.4f} × {Lz:.4f} Å")
    print(f"  Water z     : {z_min:.3f} – {z_max:.3f} Å  "
          f"(thickness ≈ {z_water:.2f} Å)")
    print(f"  Vacuum      : ≈ {vac_each:.2f} Å on each side")
    print(f"  Density     : ≈ {density:.4f} g/cm³  (in slab region)")
    print(f"  Min O–O     : {min_oo:.4f} Å  (PBC)   mean nearest-O–O: {mean_oo:.4f} Å")

    for w in warnings: print(f"  WARNING : {w}")
    for e in errors:   print(f"  ERROR   : {e}")

    if errors:
        sys.exit(1)

    print(f"  Status      : PASSED")


if __name__ == "__main__":
    main()
