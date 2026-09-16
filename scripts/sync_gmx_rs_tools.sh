#!/usr/bin/env bash
# Refresh the vendored gmx-rs-tools sources used by s_mmpbsa.
#
# The crate lives in <repo>/gmx-rs-tools/ and is linked into the s_mmpbsa
# binary; src/gmx.rs is the code that drives it. This script copies the latest
# source files from a local gmx-rs-tools checkout and, when the copy carries
# local changes, re-applies gmx-rs-tools/local.patch (see gmx-rs-tools/VENDORED.md).
#
# Usage:
#   scripts/sync_gmx_rs_tools.sh [path-to-gmx-rs-tools]
# Example:
#   scripts/sync_gmx_rs_tools.sh /run/media/Data1/Projects/gmx-rs-tools

set -euo pipefail

gmx_rs_dir="${1:-/run/media/Data1/Projects/gmx-rs-tools}"

if [[ ! -f "${gmx_rs_dir}/Cargo.toml" || ! -d "${gmx_rs_dir}/src" ]]; then
    echo "Error: cannot find a gmx-rs-tools checkout at '${gmx_rs_dir}'." >&2
    exit 1
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"
dest="${repo_root}/gmx-rs-tools"

rm -rf "${dest}/src" "${dest}/tests" "${dest}/scripts"
cp -r "${gmx_rs_dir}/src" "${dest}/src"
cp -r "${gmx_rs_dir}/tests" "${dest}/tests"
cp -r "${gmx_rs_dir}/scripts" "${dest}/scripts"
cp -f "${gmx_rs_dir}/Cargo.toml" "${gmx_rs_dir}/README.md" "${gmx_rs_dir}/NOTES.md" "${dest}/"
echo "Updated gmx-rs-tools/{src,tests,scripts,Cargo.toml,README.md,NOTES.md}"

if [[ -s "${dest}/local.patch" ]]; then
    patch -p1 -d "${dest}" --forward < "${dest}/local.patch"
    echo "Re-applied gmx-rs-tools/local.patch"
else
    echo "No local patch to apply (the copy tracks the checkout verbatim)"
fi

echo
echo "Done. Rebuild with: cargo build --release"
