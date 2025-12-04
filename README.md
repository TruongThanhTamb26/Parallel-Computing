Lab2 - Improved sequential and MPI parallel implementations

Files:

- `Lab2_sequence.cpp` - optimized sequential version (no MPI dependency). Uses flat arrays and a single swap per timestep.
- `Lab2_parallel.cpp` - improved MPI-parallel implementation using a 2D Cartesian communicator, halo exchange, scatter/gather and reductions.
- `Makefile` - targets `seq` and `parallel` to build the two executables.

Build:

1. Sequential (no MPI required):

   make seq

   This produces `Lab2_sequence`.

2. Parallel (requires MPI, e.g., OpenMPI or MPICH):

   make parallel

   This produces `Lab2_parallel`.

Run:

1. Run sequential program:

   ./Lab2_sequence

2. Run parallel program (example with 4 processes):

   mpirun -np 4 ./Lab2_parallel

Notes / improvements made:

- Removed unnecessary MPI dependency from the sequential implementation and improved memory/cache usage.
- In the parallel code: validated Cartesian communicator creation, made scatter/gather safe by only providing root buffers on the root process, used unique tags for halo exchanges, swapped buffers to avoid expensive element-wise copies, and freed MPI types/communicators cleanly.
- The parallel code assumes `N` is divisible by the process grid dimensions. If not, adjust `N` or the process grid.

If you want, I can attempt to compile and run both here to validate; tell me which one to try first.
