"""
build_slab.py — Convert a cubic water box LAMMPS data file into an air/water slab.

The xy dimensions are unchanged (no rescaling).  The z box is extended to
Lz = z_mult × Lx, and the water slab is centered within that z range, creating
equal vacuum regions above and below the interface.

    Z_SHIFT = Lx × (z_mult − 1) / 2

This places both air/water interfaces in the x–y plane, which is the standard
geometry for vibrational sum-frequency generation (vSFG) spectroscopy simulations
using the MB-pol potential.

Usage (internal — called by make_slab.sh):
    python build_slab.py --input CUBIC.data --lx LX --z-mult Z --out SLAB.data
"""

import argparse
import math
import sys


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input",  required=True)
    p.add_argument("--lx",     type=float, required=True)
    p.add_argument("--z-mult", type=float, required=True)
    p.add_argument("--out",    required=True)
    args = p.parse_args()

    Lx      = args.lx
    z_mult  = args.z_mult
    Lz_new  = Lx * z_mult
    Z_SHIFT = Lx * (z_mult - 1.0) / 2.0

    # parse
    atoms  = {}
    bonds  = []
    angles = []
    masses = {}
    n_atoms = n_bonds = n_angles = 0
    n_atom_types = n_bond_types = n_angle_types = 0
    xlo = xhi = ylo = yhi = zlo = zhi = 0.0

    with open(args.input) as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        s = lines[i].strip()

        if   s.endswith(" atoms"):          n_atoms       = int(s.split()[0])
        elif s.endswith(" bonds"):          n_bonds       = int(s.split()[0])
        elif s.endswith(" angles"):         n_angles      = int(s.split()[0])
        elif s.endswith(" atom types"):     n_atom_types  = int(s.split()[0])
        elif s.endswith(" bond types"):     n_bond_types  = int(s.split()[0])
        elif s.endswith(" angle types"):    n_angle_types = int(s.split()[0])
        elif "xlo xhi" in s: xlo, xhi = float(s.split()[0]), float(s.split()[1])
        elif "ylo yhi" in s: ylo, yhi = float(s.split()[0]), float(s.split()[1])
        elif "zlo zhi" in s: zlo, zhi = float(s.split()[0]), float(s.split()[1])

        elif s == "Masses":
            i += 2
            while i < len(lines) and lines[i].strip():
                t = lines[i].split()
                masses[int(t[0])] = float(t[1])
                i += 1
            continue

        elif s in ("Atoms # full", "Atoms"):
            i += 2
            while i < len(lines) and lines[i].strip():
                t = lines[i].split()
                if len(t) >= 7:
                    aid = int(t[0])
                    atoms[aid] = {
                        "mol":  int(t[1]),  "type": int(t[2]),  "q": float(t[3]),
                        "x":    float(t[4]), "y":   float(t[5]), "z":  float(t[6]),
                        "ix": int(t[7]) if len(t) > 7 else 0,
                        "iy": int(t[8]) if len(t) > 8 else 0,
                        "iz": int(t[9]) if len(t) > 9 else 0,
                    }
                i += 1
            continue

        elif s == "Bonds":
            i += 2
            while i < len(lines) and lines[i].strip():
                t = lines[i].split()
                bonds.append((int(t[0]), int(t[1]), int(t[2]), int(t[3])))
                i += 1
            continue

        elif s == "Angles":
            i += 2
            while i < len(lines) and lines[i].strip():
                t = lines[i].split()
                angles.append((int(t[0]), int(t[1]), int(t[2]), int(t[3]), int(t[4])))
                i += 1
            continue

        i += 1

    Lx_in = xhi - xlo
    Ly_in = yhi - ylo
    Lz_in = zhi - zlo

    # ── unwrap image flags ───────────────────────────────────────────────────
    for a in atoms.values():
        a["xu"] = a["x"] + a["ix"] * Lx_in
        a["yu"] = a["y"] + a["iy"] * Ly_in
        a["zu"] = a["z"] + a["iz"] * Lz_in

    # verify O–H bond lengths after unwrapping
    bad = sum(
        1 for _, _, a1, a2 in bonds
        if math.sqrt(sum((atoms[a1][k] - atoms[a2][k])**2 for k in ("xu","yu","zu"))) > 2.0
    )
    if bad:
        print(f"  WARNING: {bad} bonds > 2.0 Å after unwrapping — check image flags")
    else:
        print(f"  Bond check OK (all O–H < 2.0 Å after unwrapping)")

    # apply z shift; xy coordinates are unchanged (z-only extension)
    for a in atoms.values():
        a["xn"] = a["xu"] - xlo
        a["yn"] = a["yu"] - ylo
        a["zn"] = a["zu"] - zlo + Z_SHIFT

    z_vals   = [a["zn"] for a in atoms.values()]
    z_center = sum(z_vals) / len(z_vals)
    print(
        f"  Slab z range : {min(z_vals):.3f} – {max(z_vals):.3f} Å  "
        f"(center {z_center:.3f} Å)"
    )
    print(
        f"  Box          : [0, {Lz_new:.3f}] Å  "
        f"(vacuum ≈ {Z_SHIFT:.2f} Å each side)"
    )

    n_mol = n_atoms // 3
    with open(args.out, "w") as f:
        f.write(
            f"LAMMPS data file — {n_mol} H2O air/water slab  "
            f"{Lx:.4f} × {Lx:.4f} × {Lz_new:.4f} Å  "
            f"z_mult={z_mult}  (generated by build_slab.py)\n\n"
        )
        f.write(f"{n_atoms} atoms\n{n_atom_types} atom types\n")
        f.write(f"{n_bonds} bonds\n{n_bond_types} bond types\n")
        f.write(f"{n_angles} angles\n{n_angle_types} angle types\n\n")
        f.write(f"0 {Lx:.6f} xlo xhi\n")
        f.write(f"0 {Lx:.6f} ylo yhi\n")
        f.write(f"0 {Lz_new:.6f} zlo zhi\n\n")

        f.write("Masses\n\n")
        for t, m in sorted(masses.items()):
            f.write(f"{t} {m}\n")

        f.write("\nAtoms # full\n\n")
        for aid in sorted(atoms.keys()):
            a = atoms[aid]
            f.write(
                f"{aid:6d} {a['mol']:6d} {a['type']} {a['q']:8.4f} "
                f"{a['xn']:14.8f} {a['yn']:14.8f} {a['zn']:14.8f} 0 0 0\n"
            )

        f.write("\nBonds\n\n")
        for b in bonds:
            f.write(f"{b[0]} {b[1]} {b[2]} {b[3]}\n")

        f.write("\nAngles\n\n")
        for a in angles:
            f.write(f"{a[0]} {a[1]} {a[2]} {a[3]} {a[4]}\n")

    print(f"  Written: {args.out}")


if __name__ == "__main__":
    main()
