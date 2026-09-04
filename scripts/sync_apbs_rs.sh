#!/usr/bin/env bash
# Refresh the vendored Rust APBS (apbs-rs) solver sources used by s_mmpbsa.
#
# The solver crates live under <repo>/apbs/ and the in-process driver is
# src/apbs_runner.rs. This script copies the latest source files from a local
# apbs-rs checkout into the repo.
#
# Note: the vendored Cargo.toml files are kept (they are adjusted for the
# s_mmpbsa workspace), so only src/ trees and the driver file are refreshed.
#
# Usage:
#   scripts/sync_apbs_rs.sh [path-to-apbs-rs]
# Example:
#   scripts/sync_apbs_rs.sh /run/media/Data/apbs-rs

set -euo pipefail

apbs_rs_dir="${1:-/run/media/Data/apbs-rs}"

if [[ ! -f "${apbs_rs_dir}/Cargo.toml" ]]; then
    echo "Error: cannot find apbs-rs workspace at '${apbs_rs_dir}'." >&2
    exit 1
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"

for crate in apbs-generic apbs-pmgc apbs-mg; do
    if [[ ! -d "${apbs_rs_dir}/${crate}/src" ]]; then
        echo "Error: '${apbs_rs_dir}/${crate}/src' not found." >&2
        exit 1
    fi
    rm -rf "${repo_root}/apbs/${crate}/src"
    cp -r "${apbs_rs_dir}/${crate}/src" "${repo_root}/apbs/${crate}/src"
    echo "Updated apbs/${crate}/src"
done

cp -f "${apbs_rs_dir}/apbs-core/src/routines.rs" "${repo_root}/src/apbs_runner.rs"
echo "Updated src/apbs_runner.rs"

echo
echo "Done. Rebuild with: cargo build --release"
