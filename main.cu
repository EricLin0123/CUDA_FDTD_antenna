#include <cuda_runtime.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define IE 60
#define JE 100
#define KE 16

const float T = 15e-12;     // 15 picoseconds
const float dt = 0.441e-12; // 0.441 picoseconds
const float t0 = 3 * T;     // 45 picoseconds
const float epsilon = 2.2;  // dielectric constant of the Duroid substrates
const int steps = 1000;

/**
 * @brief Computes the value of a Gaussian function.
 *
 * This function calculates the value of a Gaussian function at a given time `t`,
 * with a specified width `T` and center `t0`.
 *
 * @param t The time at which to evaluate the Gaussian function.
 * @return The value of the Gaussian function at time `t`.
 *
 * The Gaussian function is defined as exp(-((t - t0) / T)^2).
 */
float Gaussian(float t)
{
    return exp(-pow((t - t0) / T, 2));
}

int main()
{
    int i, j, k;
    // D = epsilon E
    // permeability of the dielectric medium is around 1, so skip it.
    float dx[IE][JE][KE], dy[IE][JE][KE], dz[IE][JE][KE];
    float ex[IE][JE][KE], ey[IE][JE][KE], ez[IE][JE][KE];
    float hx[IE][JE][KE], hy[IE][JE][KE], hz[IE][JE][KE];
    dx[i][j][k] = dx[i][j][k] + .5 * (hz[i][j][k] - hz[i][j - 1][k] - hy[i][j][k] + hy[i][j][k - 1]);
    dy[i][j][k] = dz[i][j][k] + .5 * (hx[i][j][k] - hx[i][j][k - 1] - hz[i][j][k] + hz[i - 1][j][k]);
    dz[i][j][k] = dz[i][j][k] + .5 * (hy[i][j][k] - hy[i - 1][j][k] - hx[i][j][k] + hx[i][j - 1][k]);
    hx[i][j][k] = hx[i][j][k] + .5 * (ey[i][j][k + 1] - ey[i][j][k] - ez[i][j + 1][k] + ez[i][j][k]);
    hy[i][j][k] = hy[i][j][k] + .5 * (ez[i + 1][j][k] - ez[i][j][k] - ex[i][j][k + 1] + ex[i][j][k]);
    hz[i][j][k] = hz[i][j][k] + .5 * (ex[i][j + 1][k] - ex[i][j][k] - ey[i + 1][j][k] + ey[i][j][k]);
}