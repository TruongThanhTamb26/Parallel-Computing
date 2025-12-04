#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

# Build both programs with a small grid (fast)
GRID=64
CSV=radioactive_matrix.csv
echo "Building with GRID_N=${GRID} (fast smoke test)"
make clean >/dev/null 2>&1 || true
make seq GRID_N=${GRID}
make parallel GRID_N=${GRID}

echo "Generating ${CSV} (${GRID}x${GRID})"
python3 - <<PY
import csv
N = ${GRID}
mid = N // 2
with open("${CSV}", "w", newline="") as f:
  writer = csv.writer(f)
  for i in range(N):
    row = [0.0] * N
    if i == mid:
      row[mid] = 1e6
    writer.writerow(row)
PY

# Run sequential
echo "Running sequential..."
./Lab2_sequence ${CSV} | tee seq.out

# Run parallel with 4 processes (or fewer if not available)
NP=4
echo "Running parallel (np=${NP})..."
mpirun -np ${NP} ./Lab2_parallel ${CSV} | tee par.out

# Show outputs
echo "--- sequential output ---"
cat seq.out
echo "--- parallel output ---"
cat par.out

# Extract checksums
seqsum=$(grep -m1 "checksum(sample)" seq.out | awk -F= '{print $2}' | tr -d ' ' || echo "")
parsum=$(grep -m1 "checksum(sample)" par.out | awk -F= '{print $2}' | tr -d ' ' || echo "")

if [[ -z "$seqsum" || -z "$parsum" ]]; then
  echo "Could not find checksum in outputs; inspect seq.out and par.out"
  exit 1
fi

if [[ "$seqsum" == "$parsum" ]]; then
  echo "Checksums match: $seqsum"
  exit 0
else
  echo "Checksum mismatch: seq=$seqsum par=$parsum"
  exit 2
fi
