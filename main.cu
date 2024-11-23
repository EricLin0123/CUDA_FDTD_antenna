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
    int t = 0;
    int nsteps = 1;
    int n = 0;
    float pulse = 0.0f;
    // D = epsilon E
    // permeability of the dielectric medium is around 1, so skip it.
    float ex[IE][JE][KE], ey[IE][JE][KE], ez[IE][JE][KE];
    float hx[IE][JE][KE], hy[IE][JE][KE], hz[IE][JE][KE];

    // Initialize the fields to zero
    for (i = 0; i < IE; i++)
    {
        for (j = 0; j < JE; j++)
        {
            for (k = 0; k < KE; k++)
            {
                ex[i][j][k] = 0.0f;
                ey[i][j][k] = 0.0f;
                ez[i][j][k] = 0.0f;
                hx[i][j][k] = 0.0f;
                hy[i][j][k] = 0.0f;
                hz[i][j][k] = 0.0f;
            }
        }
    }

    // source excitation under the microstrip line
    for (int i = 21; i <= 27; i++)
    {
        for (int k = 0; k <= 3; k++)
        {
            ez[i][0][k] = Gaussian(t);
        }
    }
    while (nsteps > 0)
    {
        printf("nsteps --> ");
        scanf("%d", &nsteps);
        printf("%d \n", nsteps);
        for (n = 1; n <= nsteps; n++)
        {
            t += dt;
            // Update electric fields
            for (i = 1; i < IE; i++)
            {
                for (j = 1; j < JE; j++)
                {
                    for (k = 1; k < KE; k++)
                    {
                        ex[i][j][k] = ex[i][j][k] + .5 * (hz[i][j - 1][k] - hz[i][j][k] - hy[i][j][k] + hy[i][j][k - 1]);
                        ey[i][j][k] = ey[i][j][k] + .5 * (hx[i][j][k] - hx[i][j][k - 1] - hz[i][j][k] + hz[i - 1][j][k]);
                        ez[i][j][k] = ez[i][j][k] + .5 * (hy[i][j][k] - hy[i - 1][j][k] - hx[i][j][k] + hx[i][j - 1][k]);
                    }
                }
            }

            pulse = Gaussian(t);
            // Source excitation
            for (i = 21; i <= 27; i++)
            {
                for (k = 0; k <= 3; k++)
                {
                    ez[i][0][k] = pulse;
                }
            }

            // Update magnetic fields
            for (i = 0; i < IE - 1; i++)
            {
                for (j = 0; j < JE - 1; j++)
                {
                    for (k = 0; k < KE - 1; k++)
                    {
                        hx[i][j][k] = hx[i][j][k] + .5 * (ey[i][j][k + 1] - ey[i][j][k] - ez[i][j + 1][k] + ez[i][j][k]);
                        hy[i][j][k] = hy[i][j][k] + .5 * (ez[i + 1][j][k] - ez[i][j][k] - ex[i][j][k + 1] + ex[i][j][k]);
                        hz[i][j][k] = hz[i][j][k] + .5 * (ex[i][j + 1][k] - ex[i][j][k] - ey[i + 1][j][k] + ey[i][j][k]);
                    }
                }
            }
        }
    }
}