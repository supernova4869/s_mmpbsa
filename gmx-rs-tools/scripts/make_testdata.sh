#!/bin/sh
# Builds the test data used by the optional parity tests from the systems that
# ship in the GROMACS source tree.
#
# Usage: scripts/make_testdata.sh [output-directory]
#
# The result contains topol.tpr and traj.xtc, which is what
# `GMXRS_PARITY_TESTDATA` expects.  A GROMACS installation (gmx + data prefix)
# is required.
set -eu

out=${1:-/tmp/gmx-rs-tools/work}
src=$(cd "$(dirname "$0")/../../.." && pwd)
gmx=${GMX_BIN:-gmx}

mkdir -p "$out"
cp "$src/src/testutils/simulationdatabase/alanine_vsite_solvated."* "$out/"
cp "$src/src/testutils/simulationdatabase/alanine_vsite.itp" "$out/"

cat > "$out/test.mdp" <<'MDPMDP'
integrator  = md
dt          = 0.002
nsteps      = 500
nstxout     = 100
nstvout     = 100
nstxout-compressed = 100
nstenergy   = 100
nstlog      = 100
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 0.9
rvdw        = 0.9
pbc         = xyz
tcoupl      = no
gen_vel     = yes
gen_temp    = 300
gen_seed    = 12345
MDPMDP

( cd "$out" && "$gmx" grompp -f test.mdp -c alanine_vsite_solvated.gro \
    -p alanine_vsite_solvated.top -o topol.tpr -maxwarn 5 >/dev/null )
( cd "$out" && printf '0\n' | "$gmx" trjconv -f alanine_vsite_solvated.xtc \
    -s topol.tpr -o traj.xtc >/dev/null 2>&1 )

echo "test data ready in $out"
