#!/bin/sh
# Created by: Kamila Riedlová (parts of the implementation were taken over from DODO), some edits were assisted by ChatGPT5.
# File: run_all.sh

set -eu # -e: fail on error, -u: no unset vars

PARAM_DIR="/data/docking_parameters"

# Optional runtime filters (can be set via CLI flags below)
MAX_FILES=""         # e.g. --max 3
PATTERN=""           # e.g. --pattern '1a4u|1a54|1g24'

# --- args ---
while [ $# -gt 0 ]; do
  case "$1" in
    --max)       MAX_FILES="${2:-}"; shift 2 ;;
    --pattern)   PATTERN="${2:-}";   shift 2 ;;
    *) echo "Unknown arg: $1" >&2; exit 2 ;;
  esac
done

echo "Looking for JSONs in: $PARAM_DIR"
if [ ! -d "$PARAM_DIR" ]; then
  echo "Parameter directory not found: $PARAM_DIR" >&2
  exit 1
fi

if [ -n "$PATTERN" ]; then
  set -- $(ls "$PARAM_DIR"/*.json 2>/dev/null | sort | grep -E "$PATTERN" || true)
else
  set -- $(ls "$PARAM_DIR"/*.json 2>/dev/null | sort || true)
fi

if [ $# -eq 0 ]; then
  echo "No JSON files found."
  exit 0
fi

if [ -n "${MAX_FILES}" ]; then
  set -- $(printf '%s\n' "$@" | head -n "$MAX_FILES")
fi

count=$#
echo "Will run $count job(s)."

kill_child() {
  if [ -n "${CHILD_PID:-}" ] && kill -0 "$CHILD_PID" 2>/dev/null; then
    echo "Stopping current job (PID $CHILD_PID) ..."
    kill "$CHILD_PID" 2>/dev/null || true
    sleep 1
    kill -9 "$CHILD_PID" 2>/dev/null || true
  fi
  echo "Aborting run_all.sh"
  exit 130
}
trap kill_child INT TERM

i=0
for j in "$@"; do
  i=$((i+1))
  echo "==> [$i/$count] Docking: $j"
  python3 /app/run_docking.py "$j" &
  CHILD_PID=$!
  wait "$CHILD_PID" || echo "[WARN] Job failed for $j — continuing."
  unset CHILD_PID
done

echo "All done."
