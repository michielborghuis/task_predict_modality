#!/usr/bin/env bash
set -euo pipefail
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
dep="$repo_root/.dependencies/datasets-106edbde"
commit=106edbdeae9775d71bdfb4e2dcca67e187dc843b
mkdir -p "$repo_root/.dependencies"
if [ ! -d "$dep/.git" ]; then
  git clone --filter=blob:none --no-checkout https://github.com/openproblems-bio/datasets.git "$dep"
  git -C "$dep" checkout --detach "$commit"
fi
actual="$(git -C "$dep" rev-parse HEAD)"
if [ "$actual" != "$commit" ]; then
  echo "Expected datasets commit $commit, found $actual" >&2
  exit 1
fi
git -C "$dep" diff --quiet HEAD -- src/normalization/log_cp src/api _viash.yaml || {
  echo "Pinned normalization source has local changes; refusing to build" >&2
  exit 1
}
cd "$repo_root"
target="${1:-$dep/target_viash097}"
viash ns build --src "$dep/src/normalization/log_cp" --target "$target" \
  --query '^normalization/log_cp$' --runner 'executable|nextflow' --engine native \
  --config_mod '.engines += [{"type":"native"}]'
echo "Verified and built datasets dependency at $actual"
