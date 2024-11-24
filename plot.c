#include <stdio.h>
#include <math.h>

#define PI 3.141592653589793
#define SIGMA 1.0
#define MU 0.0

double gaussian(double x)
{
    return (1.0 / (SIGMA * sqrt(2 * PI))) * exp(-0.5 * pow((x - MU) / SIGMA, 2));
}

int main()
{
    FILE *fp = fopen("junk/gaussian_data.txt", "w");
    if (fp == NULL)
    {
        perror("Error opening file");
        return 1;
    }

    double x;
    for (x = -5.0; x <= 5.0; x += 0.1)
    {
        fprintf(fp, "%f %f\n", x, gaussian(x));
    }
    fclose(fp);

    // Plot using Gnuplot
    FILE *gp = popen("gnuplot -persistent", "w"); // open a pipe to gnuplot
    if (gp == NULL)
    {
        perror("Error opening Gnuplot");
        return 1;
    }

    fprintf(gp, "plot 'junk/gaussian_data.txt' with lines title 'Gaussian Curve'\n");
    pclose(gp);

    return 0;
}
