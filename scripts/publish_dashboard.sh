#!/usr/bin/env bash
# Build the read-only dashboard SPA for this workspace (see the companion
# .github/workflows/publish-dashboard.yml). Scaffolded by
# `vivarium-workbench add-dashboard`; run it locally to preview the bundle.
set -euo pipefail
WS_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="${1:-$WS_ROOT/reports/published/workbench}"
BASE_PATH="/viva-simbio/workbench"
INTERACTIVE_URL="https://github.com/vivarium-collective/viva-simbio"
rm -rf "$OUT"
PYTHONPATH="$WS_ROOT${PYTHONPATH:+:$PYTHONPATH}" \
  vivarium-workbench-publish \
    --workspace "$WS_ROOT" \
    --out "$OUT" \
    --base-path "$BASE_PATH" \
    --interactive-url "$INTERACTIVE_URL"
find "$OUT" -name '*.map' -delete
touch "$OUT/.nojekyll"
echo "built read-only workbench bundle at $OUT ($(du -sh "$OUT" | cut -f1))"
