#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

// Constants (N can be overridden at compile time: -DN=...)
#ifndef N
#define N 4000
#endif
const int D = 1000;
const double k = 0.0001;
const double lambda = 0.00003;
const double ux = 3.3;
const double uy = 1.4;
const double dx = 10.0;
const double dy = 10.0;
const double dt = 1.0;
const int ITERATIONS = 100;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // -------------------------------------------------------------------------
    // 1. Setup Cartesian Topology (Grid of Processes)
    // -------------------------------------------------------------------------
    
    // Let MPI decide best dimensions (rows x cols), attempt near-square grid
    int dims[2] = {0, 0};
    MPI_Dims_create(world_size, 2, dims);

    int periods[2] = {0, 0}; // non-periodic
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm);

    if (cart_comm == MPI_COMM_NULL) {
        if (world_rank == 0) std::cerr << "Error creating Cartesian communicator" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int cart_rank;
    int coords[2]; // coords[0] = row, coords[1] = col
    MPI_Comm_rank(cart_comm, &cart_rank);
    MPI_Cart_coords(cart_comm, cart_rank, 2, coords);

    // Find neighbors (North, South, West, East)
    // MPI_Cart_shift handles boundary conditions automatically (returns MPI_PROC_NULL)
    int top, bottom, left, right;
    MPI_Cart_shift(cart_comm, 0, 1, &top, &bottom); // Shift along 0 (Row) -> Top/Bottom
    MPI_Cart_shift(cart_comm, 1, 1, &left, &right); // Shift along 1 (Col) -> Left/Right

    // Calculate Local Block Size
    if (dims[0] <= 0 || dims[1] <= 0 || dims[0] * dims[1] != world_size) {
        if (cart_rank == 0) std::cerr << "Invalid process grid dims." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (N % dims[0] != 0 || N % dims[1] != 0) {
        if (cart_rank == 0) std::cerr << "Error: N must be divisible by grid dims (" << dims[0] << "x" << dims[1] << ")" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int block_rows = N / dims[0];
    int block_cols = N / dims[1];

    // Local grid includes Halo (Ghost) cells on ALL 4 sides
    // Size: (block_rows + 2) x (block_cols + 2)
    int local_rows_halo = block_rows + 2;
    int local_cols_halo = block_cols + 2;
    std::vector<double> local_grid(local_rows_halo * local_cols_halo, 0.0);
    std::vector<double> local_new_grid(local_rows_halo * local_cols_halo, 0.0);

    // -------------------------------------------------------------------------
    // 2. Prepare and Scatter Data (Data Rearrangement)
    // -------------------------------------------------------------------------
    
    // Buffer on root to hold the re-ordered data for scattering
    std::vector<double> send_buffer;
    if (cart_rank == 0) {
        std::vector<double> global_grid(static_cast<size_t>(N) * N, 0.0);
        // place source in the center
        int mid = N / 2;
        global_grid[mid * N + mid] = 1e6;

        send_buffer.resize(static_cast<size_t>(N) * N);

        // Pack blocks (row-major over process rows/cols)
        size_t p = 0;
        for (int pr = 0; pr < dims[0]; ++pr) {
            for (int pc = 0; pc < dims[1]; ++pc) {
                for (int r = 0; r < block_rows; ++r) {
                    int global_r = pr * block_rows + r;
                    size_t base = static_cast<size_t>(global_r) * N + pc * block_cols;
                    for (int c = 0; c < block_cols; ++c) {
                        send_buffer[p++] = global_grid[base + c];
                    }
                }
            }
        }
    }

    // Buffer to receive one scattered block (real data only, no ghosts)
    std::vector<double> recv_block(block_rows * block_cols);

    // REQUIREMENT: MPI_Scatter
    // Provide nullptr on non-root for clarity
    double* send_ptr = (cart_rank == 0 && !send_buffer.empty()) ? send_buffer.data() : nullptr;
    MPI_Scatter(send_ptr, block_rows * block_cols, MPI_DOUBLE,
                recv_block.data(), block_rows * block_cols, MPI_DOUBLE,
                0, cart_comm);

    // Copy scattered data into the center of local_grid (surrounded by halo)
    // We map 1D recv_block to the 2D local_grid structure
    for (int r = 0; r < block_rows; r++) {
        for (int c = 0; c < block_cols; c++) {
            // local_grid indices start at (1, 1) to skip top and left halo
            local_grid[(r + 1) * local_cols_halo + (c + 1)] = recv_block[r * block_cols + c];
        }
    }

    // -------------------------------------------------------------------------
    // 3. Define MPI Datatype for Column Exchange (Non-contiguous memory)
    // -------------------------------------------------------------------------
    MPI_Datatype col_type;
    MPI_Type_vector(block_rows, 1, local_cols_halo, MPI_DOUBLE, &col_type);
    MPI_Type_commit(&col_type);

    // Synchronize before timing
    MPI_Barrier(cart_comm);
    double start_time = MPI_Wtime();

    // ------------------ Simulation Loop -----------------------
    for (int t = 0; t < ITERATIONS; t++) {

        // --- REQUIREMENT: Exchange Boundary Values (Halo Exchange) ---
        // We need to exchange 4 directions: Top, Bottom, Left, Right
    MPI_Request reqs[8];
    int num_reqs = 0;

    // Use unique tags per direction to avoid any matching surprises
    const int TAG_TOP = 101, TAG_BOTTOM = 102, TAG_LEFT = 103, TAG_RIGHT = 104;

    // Up/Down (rows contiguous) - post only to real neighbors
    if (top != MPI_PROC_NULL) {
        MPI_Isend(&local_grid[1 * local_cols_halo + 1], block_cols, MPI_DOUBLE, top, TAG_TOP, cart_comm, &reqs[num_reqs++]);
        MPI_Irecv(&local_grid[0 * local_cols_halo + 1], block_cols, MPI_DOUBLE, top, TAG_BOTTOM, cart_comm, &reqs[num_reqs++]);
    }

    if (bottom != MPI_PROC_NULL) {
        MPI_Isend(&local_grid[block_rows * local_cols_halo + 1], block_cols, MPI_DOUBLE, bottom, TAG_BOTTOM, cart_comm, &reqs[num_reqs++]);
        MPI_Irecv(&local_grid[(block_rows + 1) * local_cols_halo + 1], block_cols, MPI_DOUBLE, bottom, TAG_TOP, cart_comm, &reqs[num_reqs++]);
    }

    // Left/Right (columns, non-contiguous)
    if (left != MPI_PROC_NULL) {
        MPI_Isend(&local_grid[1 * local_cols_halo + 1], 1, col_type, left, TAG_LEFT, cart_comm, &reqs[num_reqs++]);
        MPI_Irecv(&local_grid[1 * local_cols_halo + 0], 1, col_type, left, TAG_RIGHT, cart_comm, &reqs[num_reqs++]);
    }

    if (right != MPI_PROC_NULL) {
        MPI_Isend(&local_grid[1 * local_cols_halo + block_cols], 1, col_type, right, TAG_RIGHT, cart_comm, &reqs[num_reqs++]);
        MPI_Irecv(&local_grid[1 * local_cols_halo + block_cols + 1], 1, col_type, right, TAG_LEFT, cart_comm, &reqs[num_reqs++]);
    }

    if (num_reqs > 0) MPI_Waitall(num_reqs, reqs, MPI_STATUS_IGNORE);

        // --- Computation ---
        long long local_uncontaminated = 0;

        for (int i = 1; i <= block_rows; i++) {
            for (int j = 1; j <= block_cols; j++) {
                int idx = i * local_cols_halo + j;
                double c_curr = local_grid[idx];

                // Access neighbors (Halo cells are populated now)
                double c_up    = local_grid[(i - 1) * local_cols_halo + j];
                double c_down  = local_grid[(i + 1) * local_cols_halo + j];
                double c_left  = local_grid[i * local_cols_halo + (j - 1)];
                double c_right = local_grid[i * local_cols_halo + (j + 1)];

                // Diffusion
                double d2C_dx2 = (c_down - 2 * c_curr + c_up) / (dx * dx); // Be careful with orientation, usually i is Y (rows), j is X (cols)
                double d2C_dy2 = (c_right - 2 * c_curr + c_left) / (dy * dy);
                
                // Advection (Upwind scheme)
                double dC_dx = (c_curr - c_up) / dx; 
                double dC_dy = (c_curr - c_left) / dy;

                double delta_C = D * (d2C_dx2 + d2C_dy2) - (ux * dC_dx + uy * dC_dy) - (lambda + k) * c_curr;
                
                double val = c_curr + dt * delta_C;
                local_new_grid[idx] = val;

                if (val <= 1e-20) local_uncontaminated++;
            }
        }

        // --- REQUIREMENT: MPI Reduce ---
        long long global_uncontaminated = 0;
        MPI_Reduce(&local_uncontaminated, &global_uncontaminated, 1, MPI_LONG_LONG, MPI_SUM, 0, cart_comm);

        // --- Update ---
        // Fast swap entire buffers (ghosts will be overwritten on next exchange)
        local_grid.swap(local_new_grid);

        // --- REQUIREMENT: MPI Barrier ---
        MPI_Barrier(cart_comm);
    }

    // Broadcast stop (informational)
    int stop = 1;
    MPI_Bcast(&stop, 1, MPI_INT, 0, cart_comm);

    // -------------------------------------------------------------------------
    // 4. Gather Results
    // -------------------------------------------------------------------------
    
    // Pack real data back into contiguous block for Gathering (remove halo)
    for (int r = 0; r < block_rows; r++) {
        for (int c = 0; c < block_cols; c++) {
            recv_block[r * block_cols + c] = local_grid[(r + 1) * local_cols_halo + (c + 1)];
        }
    }

    // Gather blocks back to send_buffer on Rank 0
    // Prepare pointer for gather: recvbuf only valid on root
    double* gather_recv_ptr = (cart_rank == 0 && !send_buffer.empty()) ? send_buffer.data() : nullptr;
    MPI_Gather(recv_block.data(), block_rows * block_cols, MPI_DOUBLE,
               gather_recv_ptr, block_rows * block_cols, MPI_DOUBLE,
               0, cart_comm);

    // Reconstruct global map from blocks on Rank 0
    if (cart_rank == 0) {
        std::vector<double> final_grid(N * N);
        int p = 0;
        // The data in send_buffer is ordered by Process Rank (Block 0, Block 1...)
        // We need to map it back to Global Grid (Row 0...N)
        for (int pr = 0; pr < dims[0]; pr++) {
            for (int pc = 0; pc < dims[1]; pc++) {
                for (int r = 0; r < block_rows; r++) {
                    for (int c = 0; c < block_cols; c++) {
                        int global_r = pr * block_rows + r;
                        int global_c = pc * block_cols + c;
                        final_grid[global_r * N + global_c] = send_buffer[p++];
                    }
                }
            }
        }
        
        double end_time = MPI_Wtime();
        std::cout << "Parallel Time (2D Blocks): " << end_time - start_time << " s" << std::endl;
    }

    MPI_Type_free(&col_type);
    MPI_Comm_free(&cart_comm);
    MPI_Finalize();
    return 0;
}