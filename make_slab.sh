#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# make_slab.sh
#
# Generate a LAMMPS data file for an MB-pol air/water slab simulation.
# Outputs atom_style full with MB-pol charges, ready for use with the MBX
# pair style in LAMMPS.
#
# When called with NO flags all parameters are read from config.conf in the
# same directory as this script.  Any flag overrides the corresponding config
# value for that run without modifying the file.
#
# Usage
#   bash make_slab.sh                        # use config.conf defaults
#   bash make_slab.sh -n 256 -r 1.0 -z 3    # fully specified
#   bash make_slab.sh -n 256                 # override N only
#
# Flags
#   -n INT    Number of water molecules
#   -r FLOAT  Bulk density in g/cm³  (sets Lx = Ly = Lz_cubic)
#   -z FLOAT  Z multiplier: Lz = z_mult × Lx
#   -T FLOAT  NVT temperature in K                [default: 298.15]
#   -S INT    Equilibration run length in steps   [default: 2500000]
#   -d FLOAT  MD timestep in fs                   [default: 0.2]
#   -o DIR    Output directory                    [default: current dir]
#   -s INT    Packmol random seed                 [default: 42]
#   -t PATH   Path to water molecule template PDB [default: water.pdb]
#   -P PATH   Python interpreter to use           [default: auto-detect]
#   -c PATH   Config file to load                 [default: config.conf]
#   -h        Print this help and exit
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# helpers
usage() { grep "^#" "$0" | sed 's/^# \{0,3\}//' | sed '1d'; exit 0; }
die()   { echo "ERROR: $*" >&2; exit 1; }

# locate config
CONFIG_FILE="${SCRIPT_DIR}/config.conf"
for ((idx=1; idx<=$#; idx++)); do
    if [[ "${!idx}" == "-c" ]]; then
        nxt=$((idx+1)); CONFIG_FILE="${!nxt}"; break
    fi
done

# load config
if [[ -f "$CONFIG_FILE" ]]; then
    while IFS='=' read -r key raw_val; do
        [[ "$key" =~ ^[[:space:]]*# || -z "${key// }" ]] && continue
        key="${key//[[:space:]]/}"
        val="${raw_val%%#*}"; val="${val//[[:space:]]/}"
        case "$key" in
            NWATERS|DENSITY|Z_MULT|SEED|TEMP|EQUIL_STEPS|TIMESTEP|\
            WATER_TEMPLATE|OUTDIR|PYTHON)
                declare "$key"="$val" ;;
        esac
    done < "$CONFIG_FILE"
    echo "Loaded config: $CONFIG_FILE"
fi

# parse flags (override config values)
while getopts ":n:r:z:T:S:d:o:s:t:P:c:h" opt; do
    case "$opt" in
        n) NWATERS="$OPTARG"        ;;
        r) DENSITY="$OPTARG"        ;;
        z) Z_MULT="$OPTARG"         ;;
        T) TEMP="$OPTARG"           ;;
        S) EQUIL_STEPS="$OPTARG"    ;;
        d) TIMESTEP="$OPTARG"       ;;
        o) OUTDIR="$OPTARG"         ;;
        s) SEED="$OPTARG"           ;;
        t) WATER_TEMPLATE="$OPTARG" ;;
        P) PYTHON="$OPTARG"         ;;
        c) ;;
        h) usage                    ;;
        :) die "Flag -$OPTARG requires an argument." ;;
       \?) die "Unknown flag: -$OPTARG.  Use -h for help." ;;
    esac
done

# validate required parameters
[[ -z "${NWATERS:-}" ]] && die "Number of waters not set. Use -n or set NWATERS in config.conf."
[[ -z "${DENSITY:-}" ]] && die "Density not set. Use -r or set DENSITY in config.conf."
[[ -z "${Z_MULT:-}"  ]] && die "Z multiplier not set. Use -z or set Z_MULT in config.conf."

# apply defaults for optional parameters
SEED="${SEED:-42}"
OUTDIR="${OUTDIR:-.}"
TEMP="${TEMP:-298.15}"
EQUIL_STEPS="${EQUIL_STEPS:-2500000}"
TIMESTEP="${TIMESTEP:-0.2}"

# Water template: same directory as script
WATER_TEMPLATE="${WATER_TEMPLATE:-${SCRIPT_DIR}/water.pdb}"
[[ -f "$WATER_TEMPLATE" ]] || die "Water template not found: $WATER_TEMPLATE"

# Python: active venv > VIRTUAL_ENV > python3 on PATH > python on PATH
if [[ -z "${PYTHON:-}" ]]; then
    if [[ -n "${VIRTUAL_ENV:-}" ]]; then
        PYTHON="${VIRTUAL_ENV}/bin/python"
    elif command -v python3 &>/dev/null; then
        PYTHON="$(command -v python3)"
    elif command -v python &>/dev/null; then
        PYTHON="$(command -v python)"
    else
        die "Python not found. Activate a virtual environment, install Python 3, or set PYTHON in config.conf."
    fi
fi
[[ -x "$PYTHON" ]] || PYTHON="$(command -v "$PYTHON" 2>/dev/null || echo "$PYTHON")"

# compute Lx and Lz from N and density
read -r LX LZ < <(
    "$PYTHON" - <<PYEOF
import math
N, rho, z = $NWATERS, $DENSITY, $Z_MULT
Lx = (N * 18.015 / (rho * 6.022e23) * 1e24) ** (1.0 / 3.0)
print(f"{Lx:.6f} {Lx * z:.6f}")
PYEOF
)

# output path
mkdir -p "$OUTDIR"
OUTFILE="${OUTDIR}/water_slab_${NWATERS}H2O_rho${DENSITY}_z${Z_MULT}x.data"

# PID-scoped temp files so parallel runs don't collide
TMP="${OUTDIR}/.mbxslab_$$"
PACK_INP="${TMP}.pack.inp"
PACK_PDB="${TMP}.packed.pdb"
INIT_DATA="${TMP}.initial.data"
trap 'rm -f "$PACK_INP" "$PACK_PDB" "$INIT_DATA"' EXIT

# banner
echo ""
echo "  ┌─────────────────────────────────────────────────────────┐"
echo "  │              mbx-water-slab-prep                        │"
echo "  └─────────────────────────────────────────────────────────┘"
printf "  %-18s %s\n"  "N waters:"     "$NWATERS"
printf "  %-18s %s\n"  "density:"      "${DENSITY} g/cm³"
printf "  %-18s %s\n"  "Lx = Ly:"      "${LX} Å"
printf "  %-18s %s\n"  "Lz:"           "${LZ} Å  (z_mult = ${Z_MULT})"
printf "  %-18s %s\n"  "seed:"         "$SEED"
printf "  %-18s %s\n"  "temperature:"  "${TEMP} K"
printf "  %-18s %s\n"  "timestep:"     "${TIMESTEP} fs"
printf "  %-18s %s\n"  "equil steps:"  "$EQUIL_STEPS"
printf "  %-18s %s\n"  "output:"       "$OUTFILE"
echo ""

# step 1: pack
echo "── [1/5] Packing with Packmol ─────────────────────────────────"
"$PYTHON" "${SCRIPT_DIR}/pack.py" \
    --nwaters  "$NWATERS"        \
    --lx       "$LX"             \
    --seed     "$SEED"           \
    --template "$WATER_TEMPLATE" \
    --out-inp  "$PACK_INP"       \
    --out-pdb  "$PACK_PDB"

# step 2: pdb → lammps
echo ""
echo "── [2/5] Converting PDB → LAMMPS data file ────────────────────"
"$PYTHON" "${SCRIPT_DIR}/pdb_to_lammps.py" \
    --pdb     "$PACK_PDB"  \
    --nwaters "$NWATERS"   \
    --lx      "$LX"        \
    --out     "$INIT_DATA"

# step 3: build slab
echo ""
echo "── [3/5] Building slab (z_mult = ${Z_MULT}) ───────────────────"
"$PYTHON" "${SCRIPT_DIR}/build_slab.py" \
    --input  "$INIT_DATA" \
    --lx     "$LX"        \
    --z-mult "$Z_MULT"    \
    --out    "$OUTFILE"

# step 4: validate
echo ""
echo "── [4/5] Validating ────────────────────────────────────────────"
"$PYTHON" "${SCRIPT_DIR}/validate.py" \
    --slab    "$OUTFILE" \
    --nwaters "$NWATERS" \
    --lx      "$LX"

# step 5: generate LAMMPS equilibration input
echo ""
echo "── [5/5] Generating LAMMPS input and mbx.json ─────────────────"
"$PYTHON" "${SCRIPT_DIR}/gen_input.py" \
    --data-file   "$OUTFILE"     \
    --nwaters     "$NWATERS"     \
    --lx          "$LX"          \
    --lz          "$LZ"          \
    --temp        "$TEMP"        \
    --equil-steps "$EQUIL_STEPS" \
    --timestep    "$TIMESTEP"    \
    --outdir      "$OUTDIR"

echo ""
echo "  Data file  → ${OUTFILE}"
echo "  LAMMPS in  → ${OUTDIR}/in.equil.lammps"
echo "  MBX config → ${OUTDIR}/mbx.json"
echo ""
