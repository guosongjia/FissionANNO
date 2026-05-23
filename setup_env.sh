#!/usr/bin/env bash
# FissionANNO env installer — reproducible from a fresh machine.
#
# Handles every pitfall encountered during initial setup:
#   1. The user's ~/.condarc may carry third-party channels (cruizperez,
#      r/main) and `channel_priority: strict` that explode the solver.
#      We write a minimal sidecar condarc and pass it via $CONDARC for
#      this install only — your real ~/.condarc is left untouched.
#   2. `cigar` (lifton transitive dep) ships a broken `ez_setup.py`
#      bootstrap. We download the tarball, strip the bootstrap, and
#      build a wheel manually with --no-build-isolation.
#   3. lifton 1.0.2 imports `pytest` at runtime (via liftoff/tests)
#      but does not declare it; pip-install pytest explicitly.
#
# Usage:
#   bash setup_env.sh                  # default prefix /data/c/jiaguosong/conda_envs/fissionanno
#   ENV_PREFIX=/path bash setup_env.sh # override install location
#
# Requires:
#   - micromamba (or mamba) on PATH, OR set MICROMAMBA to the binary path
#   - curl, tar, sed (standard on Linux)

set -euo pipefail

ENV_PREFIX="${ENV_PREFIX:-/data/c/jiaguosong/conda_envs/fissionanno}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
YAML="$REPO_ROOT/workflow/envs/fissionanno.yaml"

# ----- locate micromamba / mamba -----
if [ -n "${MICROMAMBA:-}" ]; then
    MM="$MICROMAMBA"
elif command -v micromamba >/dev/null 2>&1; then
    MM="$(command -v micromamba)"
elif command -v mamba >/dev/null 2>&1; then
    MM="$(command -v mamba)"
else
    echo "ERROR: neither micromamba nor mamba found on PATH." >&2
    echo "Install micromamba into a writable location, e.g.:" >&2
    echo "  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \\" >&2
    echo "    | tar -xj -C /path/to/install bin/micromamba" >&2
    exit 1
fi
echo "[1/5] Using solver: $MM"

# ----- sidecar condarc to neutralize user-level channels / priority -----
TMP_CONDARC_DIR="$(mktemp -d)"
TMP_CONDARC="$TMP_CONDARC_DIR/condarc"
cat > "$TMP_CONDARC" <<EOF
channels:
  - conda-forge
  - bioconda
channel_priority: flexible
EOF
trap 'rm -rf "$TMP_CONDARC_DIR"' EXIT
echo "[2/5] Wrote sidecar condarc at $TMP_CONDARC"

# ----- create env -----
if [ -e "$ENV_PREFIX" ]; then
    echo "ERROR: $ENV_PREFIX already exists. Remove it or pass ENV_PREFIX=." >&2
    exit 1
fi
echo "[3/5] Creating env at $ENV_PREFIX (this takes ~10-15 minutes)..."
CONDARC="$TMP_CONDARC" "$MM" env create \
    -y \
    --override-channels \
    -c conda-forge -c bioconda \
    -p "$ENV_PREFIX" \
    -f "$YAML"

PIP="$ENV_PREFIX/bin/pip"
PY="$ENV_PREFIX/bin/python"
[ -x "$PIP" ] || { echo "ERROR: pip not found after conda step at $PIP" >&2; exit 1; }

# ----- step (1): patched cigar wheel -----
echo "[4/5] Building patched cigar wheel..."
CIGAR_BUILD="$(mktemp -d)"
(
    cd "$CIGAR_BUILD"
    curl -sL https://files.pythonhosted.org/packages/source/c/cigar/cigar-0.1.3.tar.gz -o cigar.tar.gz
    tar xzf cigar.tar.gz
    cd cigar-0.1.3
    # Strip the broken ez_setup bootstrap; setuptools is already in env
    sed -i '/^import ez_setup/d; /^ez_setup.use_setuptools/d' setup.py
    "$PIP" install --no-build-isolation --no-deps . >/dev/null
)
rm -rf "$CIGAR_BUILD"
"$PY" -c "import cigar; print(f'  cigar OK: {cigar.__file__}')" || { echo "ERROR: cigar import failed" >&2; exit 1; }

# ----- step (2): lifton itself (cigar dep now satisfied) -----
echo "[5/5] Installing lifton 1.0.2 + pytest (lifton runtime dep)..."
"$PIP" install --no-build-isolation lifton==1.0.2 pytest >/dev/null

# ----- final verification -----
echo ""
echo "=== Final verification ==="
"$PY" --version
"$ENV_PREFIX/bin/snakemake" --version | head -1
"$ENV_PREFIX/bin/lifton" --version 2>&1 | grep -E "^v[0-9]" || echo "  lifton --version printed banner OK"
"$ENV_PREFIX/bin/miniprot" --version
"$ENV_PREFIX/bin/diamond" --version | head -1
"$ENV_PREFIX/bin/samtools" --version | head -1
"$ENV_PREFIX/bin/liftoff" --version
"$ENV_PREFIX/bin/minimap2" --version
echo ""
echo "Activate with:"
echo "  conda activate $ENV_PREFIX"
echo "  # or: micromamba activate $ENV_PREFIX"
