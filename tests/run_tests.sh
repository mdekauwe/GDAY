#!/bin/bash
#
# GDAY regression tests.
#
# Builds the model with CHECK_WATER_BALANCE (the hydraulics water balance
# is checked every timestep and the run aborts if it doesn't close), then
# runs the daily Duke example and a set of sub-daily cases (bucket,
# hydraulics, drought, cascading drainage). A case fails if the model exits
# with an error or the output contains NaN/inf.
#
# Usage: tests/run_tests.sh [nyears]   (default 3, needs python3)
#
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(dirname "$HERE")"
NYEARS="${1:-3}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/gday_tests.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

# build out of tree so the normal build isn't touched
mkdir -p "$WORK/build"
cp -r "$ROOT/src/"*.c "$ROOT/src/include" "$ROOT/src/Makefile" "$WORK/build/"
( cd "$WORK/build" && rm -f version.c && \
  make -s CFLAGS="-O2 -Wall -Wno-unused-parameter -DCHECK_WATER_BALANCE -DCHECK_NUTRIENT_BALANCE" \
  > build.log 2>&1 ) || { echo "BUILD FAILED"; cat "$WORK/build/build.log"; exit 1; }
GDAY="$WORK/build/gday"

CASES=$(python3 "$HERE/make_test_cases.py" "$ROOT/example" "$WORK/cases" \
        "$NYEARS") || { echo "could not make test cases"; exit 1; }

nfail=0
for c in $CASES; do
    if ! "$GDAY" -p "$WORK/cases/$c.cfg" > "$WORK/$c.log" 2>&1; then
        echo "FAIL  $c (model error)"; sed 's/^/      /' "$WORK/$c.log" | tail -5
        nfail=$((nfail + 1)); continue
    fi
    bad=$(grep -v '^#' "$WORK"/cases/out_"$c"*.csv | grep -ciE 'nan|inf')
    if [ "$bad" -gt 0 ]; then
        echo "FAIL  $c ($bad lines with NaN/inf)"
        nfail=$((nfail + 1)); continue
    fi
    echo "ok    $c"
done

[ "$nfail" -eq 0 ] && echo "all tests passed" || echo "$nfail test(s) failed"
exit "$nfail"
