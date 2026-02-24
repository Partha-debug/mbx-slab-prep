"""
pack.py — Generate a Packmol input file and run Packmol to pack N water molecules.

The inner packing cube is inset by BUFFER angstroms on every face so that the
minimum periodic-image atom-atom distance is guaranteed to be >= the packing
tolerance (2.0 Å by default).  This prevents the MBX induced-dipole CG solver
from failing on the first energy evaluation due to unphysically close contacts
that Packmol would otherwise place across a periodic boundary.

Requires the `packmol` pip package (pip install packmol), which places the
packmol executable on PATH inside the active Python environment.

Usage (internal — called by make_slab.sh):
    python pack.py --nwaters N --lx LX --template WATER_PDB
                   --out-inp PACK_INP --out-pdb PACK_PDB
                   [--seed INT] [--tolerance FLOAT]
"""

import argparse
import os
import shutil
import subprocess
import sys

TOLERANCE = 2.0   # Å  inter-atom distance cutoff passed to Packmol
BUFFER    = 1.0   # Å  inset from each box face  (= TOLERANCE / 2)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--nwaters",   type=int,   required=True)
    p.add_argument("--lx",        type=float, required=True)
    p.add_argument("--template",  required=True)
    p.add_argument("--out-inp",   required=True)
    p.add_argument("--out-pdb",   required=True)
    p.add_argument("--seed",      type=int,   default=42)
    p.add_argument("--tolerance", type=float, default=TOLERANCE)
    args = p.parse_args()

    # Look in the same bin/ as the running interpreter first (covers venvs
    # invoked via an absolute path without activating them), then fall back
    # to whatever is on the system PATH.
    _venv_bin = os.path.join(os.path.dirname(sys.executable), "packmol")
    if os.path.isfile(_venv_bin) and os.access(_venv_bin, os.X_OK):
        packmol_exe = _venv_bin
    else:
        packmol_exe = shutil.which("packmol")
    if packmol_exe is None:
        sys.exit(
            "ERROR: packmol not found.\n"
            "       Install it with:  pip install packmol"
        )

    buf     = BUFFER
    inner_L = args.lx - 2.0 * buf

    if inner_L <= 0:
        sys.exit(
            f"ERROR: Lx = {args.lx:.4f} Å is too small for a {buf} Å boundary buffer.\n"
            f"       Need Lx > {2 * buf:.1f} Å (equivalent to N > ~8 waters at 1 g/cm³)."
        )

    inp = (
        f"tolerance {args.tolerance}\n"
        f"output {os.path.abspath(args.out_pdb)}\n"
        f"filetype pdb\n"
        f"seed {args.seed}\n\n"
        f"structure {os.path.abspath(args.template)}\n"
        f"  number {args.nwaters}\n"
        f"  inside cube {buf:.4f} {buf:.4f} {buf:.4f} {inner_L:.4f}\n"
        f"end structure\n"
    )
    with open(args.out_inp, "w") as f:
        f.write(inp)

    print(f"  Simulation box : [0, {args.lx:.4f}] Å")
    print(f"  Packing region : [{buf:.1f}, {args.lx - buf:.4f}] Å  "
          f"(min PBC distance ≥ {args.lx - inner_L:.1f} Å)")

    with open(args.out_inp) as fh:
        result = subprocess.run(
            [packmol_exe],
            stdin=fh,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

    stdout = result.stdout.decode(errors="replace")
    stderr = result.stderr.decode(errors="replace")

    if result.returncode != 0 or not os.path.exists(args.out_pdb):
        print(stdout[-3000:])
        print(stderr[-500:])
        sys.exit("ERROR: Packmol exited with an error.")

    atom_lines = [l for l in open(args.out_pdb) if l.startswith(("ATOM", "HETATM"))]
    expected   = args.nwaters * 3
    if len(atom_lines) != expected:
        sys.exit(f"ERROR: Expected {expected} atoms in PDB output, got {len(atom_lines)}.")

    print(f"  Packmol OK — {len(atom_lines)} atoms written")


if __name__ == "__main__":
    main()
