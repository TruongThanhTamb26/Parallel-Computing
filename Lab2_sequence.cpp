// Sequential version (improved): no MPI dependency, flat arrays, fewer copies
#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>

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

int main() {
    // Use flat 1D arrays for better cache behavior and a single swap per iteration
    std::vector<double> C_old(static_cast<size_t>(N) * N, 0.0);
    std::vector<double> C_new(static_cast<size_t>(N) * N, 0.0);

    // Initialize source in the middle
    const int mid = N / 2;
    C_old[idx(mid, mid)] = 1e6;

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