Lab2 - Improved sequential and MPI parallel implementations

## Files

- `Lab2_sequence.cpp` – optimized sequential version (no MPI dependency). Uses flat arrays, CSV initialization, and a single swap per timestep.
- `Lab2_parallel.cpp` – MPI version using a 2D Cartesian communicator, halo exchange (point-to-point), MPI_Scatter/Gather, MPI_Reduce, and MPI_Bcast.
- `Makefile` – builds both targets. Accepts `GRID_N=<size>` to compile with a smaller grid (passes `-DN=<size>`).
- `run_smoke.sh` – optional helper script that builds with a tiny grid, generates a sample CSV, runs both executables, and compares checksums.

## Input format (CSV)

Both executables now require a CSV file describing the initial concentration map:

- Exactly `N` rows and `N` comma-separated values per row (where `N` is the compile-time grid size, default 4000).
- Values are interpreted as doubles. Whitespace around commas is allowed.
- Example row: `0,0,0,1000000,0,...`

If the CSV shape is invalid, the program exits with an error.

## Build

Sequential (no MPI needed):

```bash
make seq              # builds with default N=4000
make seq GRID_N=64    # builds with a smaller compile-time grid
```

Parallel (requires OpenMPI/MPICH):

```bash
make parallel
make parallel GRID_N=64
```

## Run

Provide the CSV path as the first argument in both modes:

```bash
./Lab2_sequence radioactive_matrix.csv
mpirun -np 4 ./Lab2_parallel radioactive_matrix.csv
```

The parallel program assumes `N` is divisible by the inferred process grid dimensions (`MPI_Dims_create`). Adjust `GRID_N` or the process count if needed.

## Smoke test (quick check)

The script builds everything with `GRID_N=64`, generates a centered source CSV, and compares the sequential/parallel checksums:

```bash
chmod +x run_smoke.sh
./run_smoke.sh
```

The script generates a `radioactive_matrix.csv` test map (with a central source) and feeds it to both executables.

## Notes

- Sequential improvements: flat arrays, better cache locality, CSV-driven initialization, and checksum reporting to prevent dead-code elimination.
- Parallel improvements: safe scatter/gather buffer handling, derived datatype for halo columns, ghost-layer swapping, per-iteration barriers, and clean communicator/datatype teardown.
- `GRID_N` controls both the compiled domain size and the expected CSV dimensions; keep them in sync.
