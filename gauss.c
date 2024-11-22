#include <math.h>

// Gaussian function: f(x) = (1 / (σ * sqrt(2π))) * exp(-0.5 * ((x - μ) / σ)^2)
float gaussian(float t, float t0, float T)
{
    float exponent = -pow((t - t0) / T, 2);
    return exp(exponent);
}
