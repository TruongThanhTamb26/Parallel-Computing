// Sequential version (improved): no MPI dependency, flat arrays, fewer copies
#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>

// Grid / physics constants
// Grid size can be overridden at compile-time with -DN=<value>
#ifndef N
#define N 4000
#endif
static const int D = 1000;
static const double k = 0.0001;
static const double lambda_ = 0.00003;
static const double ux = 3.3;
static const double uy = 1.4;
static const double dx = 10.0;
static const double dy = 10.0;
static const double dt = 1.0;
static const int ITERATIONS = 100;

inline size_t idx(int i, int j) { return static_cast<size_t>(i) * N + j; }

bool load_grid_from_csv(const std::string& path, std::vector<double>& grid) {
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "Failed to open CSV file: " << path << std::endl;
        return false;
    }

    grid.assign(static_cast<size_t>(N) * N, 0.0);
    std::string line;
    int row = 0;

    while (row < N && std::getline(file, line)) {
        if (line.empty()) {
            continue;
        }

        std::stringstream line_stream(line);
        std::string cell;
        int col = 0;

        while (col < N && std::getline(line_stream, cell, ',')) {
            std::stringstream cell_stream(cell);
            double value = 0.0;
            if (!(cell_stream >> value)) {
                std::cerr << "Invalid numeric value at row " << row << ", column " << col
                          << " in " << path << std::endl;
                return false;
            }
            grid[idx(row, col)] = value;
            ++col;
        }

        if (col != N) {
            std::cerr << "Row " << row << " in " << path << " has " << col
                      << " columns, expected " << N << std::endl;
            return false;
        }
        ++row;
    }

    if (row != N) {
        std::cerr << "CSV file " << path << " has " << row << " rows, expected " << N << std::endl;
        return false;
    }
    return true;
}

int main(int argc, char** argv) {
    // Use flat 1D arrays for better cache behavior and a single swap per iteration
    std::vector<double> C_old(static_cast<size_t>(N) * N, 0.0);
    std::vector<double> C_new(static_cast<size_t>(N) * N, 0.0);

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input.csv>" << std::endl;
        return 1;
    }

    if (!load_grid_from_csv(argv[1], C_old)) {
        return 1;
    }

    auto t0 = std::chrono::steady_clock::now();

    for (int t = 0; t < ITERATIONS; ++t) {
        // i = row (0..N-1), j = col (0..N-1)
        for (int i = 0; i < N; ++i) {
            size_t base = static_cast<size_t>(i) * N;
            for (int j = 0; j < N; ++j) {
                double c_curr = C_old[base + j];

                // neighbors (use 0 outside boundaries)
                double c_up    = (i > 0)     ? C_old[base - N + j] : 0.0;
                double c_down  = (i + 1 < N) ? C_old[base + N + j] : 0.0;
                double c_left  = (j > 0)     ? C_old[base + j - 1] : 0.0;
                double c_right = (j + 1 < N) ? C_old[base + j + 1] : 0.0;

                double d2C_dx2 = (c_down - 2.0 * c_curr + c_up) / (dx * dx);
                double d2C_dy2 = (c_right - 2.0 * c_curr + c_left) / (dy * dy);
                double dC_dx   = (c_curr - c_up) / dx;   // upwind simplified
                double dC_dy   = (c_curr - c_left) / dy; // upwind simplified

                double delta_C = D * (d2C_dx2 + d2C_dy2) - (ux * dC_dx + uy * dC_dy) - (lambda_ + k) * c_curr;
                C_new[base + j] = c_curr + dt * delta_C;
            }
        }

        // single swap (fast) instead of copying
        C_old.swap(C_new);
    }

    auto t1 = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = t1 - t0;
    std::cout << "Sequential elapsed time: " << elapsed.count() << " s\n";

    // print a small checksum to ensure results are used (prevents some optimizers from discarding work)
    double checksum = 0.0;
    for (size_t i = 0; i < 1000 && i < C_old.size(); ++i) checksum += C_old[i];
    std::cout << "checksum(sample) = " << checksum << "\n";

    return 0;
}