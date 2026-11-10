#!/bin/bash
# Build an offline install bundle of linux-64 packages on a machine with web access,
# for clusters whose nodes cannot reach the internet.
#   Usage: make_offline_bundle.sh MICROMAMBA_BIN OUT_DIR
# Then copy OUT_DIR to the cluster and run install_offline.sh from inside it.
set -euo pipefail
mm=${1:?usage: make_offline_bundle.sh MICROMAMBA_BIN OUT_DIR}
outdir=${2:?usage: make_offline_bundle.sh MICROMAMBA_BIN OUT_DIR}
root=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$outdir/pkgs"

# Report a typical Linux host, so a solve from macOS doesn't reject glibc packages.
export CONDA_OVERRIDE_GLIBC=2.28 CONDA_OVERRIDE_LINUX=4.18 CONDA_OVERRIDE_ARCHSPEC=x86_64
export MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX:-$outdir/.mamba}
"$mm" create -n bundle --platform linux-64 --dry-run --json -y -f "$root/environment.yml" > "$outdir/solve.json"

python3 - "$outdir" "$mm" <<'EOF'
import hashlib, json, subprocess, sys, tempfile, urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

out, mm = Path(sys.argv[1]), sys.argv[2]
links = json.load(open(out / "solve.json"))["actions"]["LINK"]

def fetch(pkg):
    name = pkg["url"].rsplit("/", 1)[1]
    dest = out / "pkgs" / name
    if not dest.exists():
        tmp = dest.with_suffix(dest.suffix + ".part")
        urllib.request.urlretrieve(pkg["url"], tmp)
        tmp.rename(dest)
    digest = hashlib.md5(dest.read_bytes()).hexdigest()
    if pkg.get("md5") and digest != pkg["md5"]:
        dest.unlink()
        raise SystemExit(f"Checksum mismatch: {name}")
    return pkg, name, digest

with ThreadPoolExecutor(8) as pool:
    done = list(pool.map(fetch, links))

# A Bioconductor data package holds no data: conda downloads its content while
# installing it, which fails without web access. Put that content in r_data/, where
# install_offline.sh installs it with R, and leave the conda package off the list.
catalog = {}
for pkg, name, _ in done:
    if pkg["name"] == "bioconductor-data-packages":
        with tempfile.TemporaryDirectory() as tmp:
            subprocess.run([mm, "package", "extract", str(out / "pkgs" / name), tmp],
                           check=True, stdout=subprocess.DEVNULL)
            catalog = json.load(open(Path(tmp) / "share/bioconductor-data-packages/dataURLs.json"))

def fetch_data(entry):
    dest = out / "r_data" / entry["fn"]
    for url in entry["urls"]:
        # curl, because bioconductor.org refuses Python's default client.
        got = subprocess.run(["curl", "-fsSL", "-o", str(dest), url]).returncode == 0
        if got and hashlib.md5(dest.read_bytes()).hexdigest() == entry["md5"]:
            return
    raise SystemExit(f"Could not download {entry['fn']}")

(out / "r_data").mkdir(exist_ok=True)
listed = []
for pkg, name, digest in done:
    entry = catalog.get(f"{pkg['name'][len('bioconductor-'):]}-{pkg['version']}")
    if pkg["name"].startswith("bioconductor-") and entry:
        fetch_data(entry)
        print(f"{pkg['name']}: content saved as r_data/{entry['fn']}")
    else:
        listed.append(f"PKGDIR/{name}#{digest}\n")
(out / "explicit.txt").write_text("@EXPLICIT\n" + "".join(listed))
print(f"{len(done)} packages downloaded")
EOF

# R packages with no conda package: source tarballs, and the script that installs them.
bash "$root/tools/install_r_source_packages.sh" fetch "$outdir/r_src"
cp "$root/tools/install_r_source_packages.sh" "$outdir/"

cat > "$outdir/install_offline.sh" <<'EOF'
#!/bin/bash
# Usage: bash install_offline.sh CONDA_BIN ENV_PREFIX
set -euo pipefail
conda_bin=${1:?usage: install_offline.sh CONDA_BIN ENV_PREFIX}
prefix=${2:?usage: install_offline.sh CONDA_BIN ENV_PREFIX}
here=$(cd "$(dirname "$0")" && pwd)
# An incomplete copy of the bundle is the usual cause of a failed install, so check
# every file before starting.
bad=""
while IFS='#' read -r name sum; do
  file="$here/pkgs/${name#PKGDIR/}"
  if [ ! -f "$file" ]; then
    bad="$bad"$'\n'"pkgs/${name#PKGDIR/}: missing"
  elif [ "$(md5sum < "$file" | cut -d' ' -f1)" != "$sum" ]; then
    bad="$bad"$'\n'"pkgs/${name#PKGDIR/}: incomplete"
  fi
done < <(grep -v '^@' "$here/explicit.txt")
shopt -s nullglob
for f in "$here"/r_data/*.tar.gz "$here"/r_src/*.tar.gz; do
  gzip -t "$f" 2>/dev/null || bad="$bad"$'\n'"${f#"$here"/}: incomplete"
done
if [ -n "$bad" ]; then
  echo "These bundle files are missing or incomplete. Copy them again, then rerun:" >&2
  echo "$bad" >&2
  exit 1
fi
echo "Bundle files verified"
sed "s#^PKGDIR/#file://$here/pkgs/#" "$here/explicit.txt" > "$here/explicit.local.txt"
"$conda_bin" create -y --offline -p "$prefix" --file "$here/explicit.local.txt"
# Bioconductor data packages, which conda could only install with web access.
for f in "$here"/r_data/*.tar.gz; do
  "$conda_bin" run --no-capture-output -p "$prefix" R CMD INSTALL "$f"
done
# Install into a copy of the tarballs, so the logs land beside the environment.
mkdir -p "$prefix/r_src"
cp "$here"/r_src/*.tar.gz "$prefix/r_src/"
"$conda_bin" run --no-capture-output -p "$prefix" \
  bash "$here/install_r_source_packages.sh" install "$prefix/r_src"
EOF
rm -f "$outdir/solve.json"
du -sh "$outdir/pkgs"
