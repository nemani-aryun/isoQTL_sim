#!/bin/bash
# isotwas and five of its dependencies have no conda package, so they are installed
# from source at fixed versions, after the conda environment exists.
#   install_r_source_packages.sh fetch DIR     download the source tarballs into DIR
#   install_r_source_packages.sh install DIR   install them with the R that is on PATH
# `fetch` needs web access and `install` does not, so a cluster without web access
# can install from a DIR that was fetched elsewhere.
set -euo pipefail
action=${1:?usage: install_r_source_packages.sh fetch|install DIR}
dir=${2:?usage: install_r_source_packages.sh fetch|install DIR}

# In install order: each package comes after the ones it needs.
packages=(
  "condMVNorm https://cran.r-project.org/package=condMVNorm&version=2025.1"
  "palasso https://cran.r-project.org/package=palasso&version=1.0.0"
  "cornet https://cran.r-project.org/package=cornet&version=1.0.0"
  "joinet https://cran.r-project.org/package=joinet&version=1.0.0"
  "RMTL https://cran.r-project.org/package=RMTL&version=1.0.0"
  "ACAT https://github.com/yaowuliu/ACAT/archive/e737d7d188ba6731800c5e6971539ca6a3888405.tar.gz"
  "isotwas https://github.com/bhattacharya-a-bt/isotwas/archive/5398e5dc27f7ccfdc3eee60fb47e20906d601294.tar.gz"
)

mkdir -p "$dir"
for entry in "${packages[@]}"; do
  name=${entry%% *}
  url=${entry#* }
  tarball="$dir/$name.tar.gz"
  case "$action" in
    fetch)
      curl -fsSL -o "$tarball.part" "$url"
      mv "$tarball.part" "$tarball"
      echo "Fetched $name"
      ;;
    install)
      test -s "$tarball" || { echo "Missing $tarball; run the fetch step first" >&2; exit 1; }
      R CMD INSTALL --no-docs "$tarball" > "$dir/$name.install.log" 2>&1 \
        || { echo "Installing $name failed. End of $dir/$name.install.log:" >&2
             tail -20 "$dir/$name.install.log" >&2; exit 1; }
      echo "Installed $name"
      ;;
    *)
      echo "Unknown action: $action" >&2; exit 1
      ;;
  esac
done
if [ "$action" = install ]; then
  Rscript -e 'suppressPackageStartupMessages(library(isotwas)); cat("isotwas", as.character(packageVersion("isotwas")), "loads\n")'
fi
