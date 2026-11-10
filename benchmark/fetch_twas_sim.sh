#!/bin/bash
# Fetch twas_sim at a fixed commit and apply the two fixes it needs to run.
# Upstream sim.py at this commit (a) ends in 40 null bytes, so Python refuses to
# load it, and (b) calls jnp.diag in estimate_her without importing jax.numpy.
# Neither fix changes the statistical method.
set -euo pipefail
dest=${1:?usage: fetch_twas_sim.sh DEST_DIR}
commit=a395a9a9eb0d9ad7f2eba6c578fe71c99a51010f
mkdir -p "$dest"
curl -fsSL "https://raw.githubusercontent.com/mancusolab/twas_sim/${commit}/sim.py" -o "$dest/sim.upstream.py"
tr -d '\000' < "$dest/sim.upstream.py" | sed 's/jnp\.diag/np.diag/' > "$dest/sim.py"
rm "$dest/sim.upstream.py"
echo "$commit" > "$dest/COMMIT"
