/* www.antenna-theory.com */
/* Pete Bevelacqua - EE 517 */

/* This program is a 3D FDTD simulation that will
    model the fields coming off of a microstrip, onto a patch antenna.
    ABC's are 1st order.
*/

// #define STRIP

#include <stdio.h>
#include <math.h>

/* dimensions in the X, Y, and Z directions */
#define LIMX 60
#define LIMY 100
#define LIMZ 16

/* time at which source is switched off and ABC turned on */
#define SWITCH1 225 // 405
#define DELAY 0

/* Total number of time steps */
#define totalT 8000

#define PI 3.14159265358979
#define MU0 1.25663706e-6
#define EPS0 8.854e-12
#define EPSR 2.2

/* globally declare fields */
__device__ double Ex[LIMX][LIMY][LIMZ], Ey[LIMX][LIMY][LIMZ], Ez[LIMX][LIMY][LIMZ];
double h_Ez[LIMX][LIMY][LIMZ];
__device__ double Hx[LIMX][LIMY][LIMZ], Hy[LIMX][LIMY][LIMZ], Hz[LIMX][LIMY][LIMZ];

/* globally declare stored field arrays for ABCs */
__device__ double HxABC1[LIMX][LIMZ], HzABC1[LIMX][LIMZ], HyABC2[LIMY][LIMZ], HzABC2[LIMY][LIMZ];
__device__ double HyABC3[LIMY][LIMZ], HzABC3[LIMY][LIMZ], HxABC4[LIMX][LIMZ], HzABC4[LIMX][LIMZ];
__device__ double HxABC5[LIMX][LIMY], HyABC5[LIMX][LIMY], ExABC6[LIMX][LIMZ], EzABC6[LIMX][LIMZ];
__device__ double ExABC5[LIMX][LIMZ], EzABC5[LIMX][LIMZ];

/* Storing the output to calculate S-parameters */
double *d_EzOut;
double *h_EzOut;

/*  I want all variables declared globally */
__device__ int i, j, k, ntime, frame = 0;

/*  Variables defining lattice and time steps, from Sheen, 1990 */
__device__ double delX, delY, delZ, delT;
__device__ double T, T0, temp;
double h_T, h_T0, h_delT;

/*  ABC Coefficients....and the FDTD coefficients */
__device__ double abcFSx, abcFSy, abcFSz, abcDIx, abcDIy, abcDIz, abcBx, abcBy, abcBz, cF, cB, cD;
__device__ double tMUX, tMUY, tMUZ, tEPX, tEPY, tEPZ, tERX, tERY, tERZ, tEBX, tEBY, tEBZ;
int hi, hj, hk, h_ntime, h_frame = 0;
FILE *out;
FILE *outGauss;
/* declaration of functions */
void Initialize();
void UpdateEfields();
void Conductors();
void Source();
void FirstABC();
void UpdateHfields();
void SecondABC();

__global__ void InitializeData(void)
{
    /* Define the Space */
    delX = 0.389e-3;
    delY = 0.400e-3;
    delZ = 0.265e-3;
    delT = 0.441e-12;

    /*  The source parameters */
    T = 15.e-12;

    T0 = 3. * T;

    /* Define Free Space ABC coefficients */
    cF = 1 / sqrt(MU0 * EPS0);
    abcFSx = (delT * cF - delX) / (delT * cF + delX);
    abcFSy = (delT * cF - delY) / (delT * cF + delY);
    abcFSz = (delT * cF - delZ) / (delT * cF + delZ);

    /* Define Dielectric ABC coefficients */
    cD = 1 / sqrt(MU0 * EPS0 * EPSR);
    abcDIx = (delT * cD - delX) / (delT * cD + delX);
    abcDIy = (delT * cD - delY) / (delT * cD + delY);
    abcDIz = (delT * cD - delZ) / (delT * cD + delZ);

    /* Define Boundary ABC coefficients */
    cB = 1 / sqrt(MU0 * EPS0 * (EPSR + 1.) / 2.);
    abcBx = (delT * cB - delX) / (delT * cB + delX);
    abcBy = (delT * cB - delY) / (delT * cB + delY);
    abcBz = (delT * cB - delZ) / (delT * cB + delZ);

    printf("abcBx = %lf, abcBy = %lf, abcBz = %lf\n", abcBx, abcBy, abcBz);

    /* Define H coefficients */
    tMUX = delT / MU0 / delX;
    tMUY = delT / MU0 / delY;
    tMUZ = delT / MU0 / delZ;

    /* E coefficients (Free Space)*/
    tEPX = delT / EPS0 / delX;
    tEPY = delT / EPS0 / delY;
    tEPZ = delT / EPS0 / delZ;

    /* E Coefficients (Dielectric) */
    tERX = delT / EPS0 / EPSR / delX;
    tERY = delT / EPS0 / EPSR / delY;
    tERZ = delT / EPS0 / EPSR / delZ;

    /* E Coefficients (Boundary) */
    tEBX = delT / EPS0 * 2. / (EPSR + 1) / delX;
    tEBY = delT / EPS0 * 2. / (EPSR + 1) / delY;
    tEBZ = delT / EPS0 * 2. / (EPSR + 1) / delZ;
    printf("Pete Rules %lf %lf %lf\n", tEBX, tEBY, tEBZ);
}

__global__ void ntimeplus1()
{
    ntime++;
}
int main()
{
    h_delT = 0.441e-12;
    h_T = 15.e-12;
    h_T0 = 3. * h_T;
    FILE *dump;
    cudaError_t err;
    char basename[80] = "junk", filename[100];
#ifdef STRIP
    char outputF[20] = "Incident_strip.txt";
#else
    char outputF[20] = "Incident.txt";
#endif
    char outputGauss[20] = "Gauss.txt";
    outGauss = fopen(outputGauss, "w+");
    out = fopen(outputF, "w+");

    InitializeData<<<1, 1>>>();
    Initialize();
    err = cudaMalloc((void **)&d_EzOut, totalT * sizeof(double));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device vector d_EzOut (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    h_EzOut = (double *)malloc(totalT * sizeof(double));

    /*Do time stepping */
    cudaEvent_t start, stop;
    float milliseconds = 0;

    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start);

    for (h_ntime = 0; h_ntime < totalT; h_ntime++)
    {
        printf("Doing time step %d\r", h_ntime);

        UpdateEfields();
        FirstABC();
        Conductors();
        Source();
        UpdateHfields();
        SecondABC();

        /* Write out E-field */
        hk = 2;
        if (h_ntime % 5 == 0)
        {
            err = cudaMemcpyFromSymbol(h_Ez, Ez, LIMX * LIMY * LIMZ * sizeof(double), 0, cudaMemcpyDeviceToHost);
            if (err != cudaSuccess)
            {
                fprintf(stderr, "Failed to copy vector Ez from device to host (error code %s)!\n", cudaGetErrorString(err));
                exit(EXIT_FAILURE);
            }
            sprintf(filename, "%s.%d", basename, h_frame++);
            dump = fopen(filename, "w");
            for (hi = 0; hi < LIMX; hi++)
                for (hj = 0; hj < LIMY; hj++)
                {
                    fprintf(dump, "%lf\n", h_Ez[hi][hj][hk]);
                }

            fclose(dump);
        }
        ntimeplus1<<<1, 1>>>();
    } /*End of time stepping*/

    cudaEventRecord(stop);
    cudaEventSynchronize(stop);

    cudaEventElapsedTime(&milliseconds, start, stop);
    printf("Time for the loop: %f ms\n", milliseconds);

    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    err = cudaMemcpy(h_EzOut, d_EzOut, totalT * sizeof(double), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy vector d_EzOut from device to host (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    for (h_ntime = 0; h_ntime < totalT; h_ntime++)
    {
        fprintf(out, "%lf\n", h_EzOut[h_ntime]);
    }

    fclose(out);
    fclose(outGauss);
}

/*  Function:  Initialize Fields   */
/**********************************
 *  Zeros all fields and ABC storage arrays *
 *****************************************/
__global__ void InitializeFields(void)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (i < LIMX && j < LIMY && k < LIMZ)
    {
        Ex[i][j][k] = 0.;
        Ey[i][j][k] = 0.;
        Ez[i][j][k] = 0.;
        Hx[i][j][k] = 0.;
        Hy[i][j][k] = 0.;
        Hz[i][j][k] = 0.;
    }
}

__global__ void InitializeABCs(void)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (i < LIMX && k < LIMZ)
    {
        if (j == 0)
        {
            HxABC1[i][k] = 0.;
            HzABC1[i][k] = 0.;
            HxABC4[i][k] = 0.;
            HzABC4[i][k] = 0.;
            ExABC6[i][k] = 0.;
            EzABC6[i][k] = 0.;
            ExABC5[i][k] = 0.;
            EzABC5[i][k] = 0.;
        }
    }

    if (j < LIMY && k < LIMZ)
    {
        if (i == 0)
        {
            HyABC2[j][k] = 0.;
            HzABC2[j][k] = 0.;
            HyABC3[j][k] = 0.;
            HzABC3[j][k] = 0.;
        }
    }

    if (i < LIMX && j < LIMY)
    {
        if (k == 0)
        {
            HxABC5[i][j] = 0.;
            HyABC5[i][j] = 0.;
        }
    }
}

void Initialize()
{
    // 2 x 100 x 4 = 800 <= 1024
    dim3 threadsPerBlock3D(2, LIMY, 4);
    // 30 x 1 x 4
    dim3 numBlocks3D((LIMX + threadsPerBlock3D.x - 1) / threadsPerBlock3D.x, (LIMY + threadsPerBlock3D.y - 1) / threadsPerBlock3D.y, (LIMZ + threadsPerBlock3D.z - 1) / threadsPerBlock3D.z);

    InitializeFields<<<numBlocks3D, threadsPerBlock3D>>>();
    InitializeABCs<<<numBlocks3D, threadsPerBlock3D>>>();
}
/*  End Initialize Function **********/

/*  Function:  UpdateEfields()    */
/**********************************
/*  Updates Ex, Ey, and Ez.
 *
 ***********************************/
__global__ void UpdateEx(void)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (i < LIMX && 1 <= j && j < LIMY - 1 && 1 <= k && k < LIMZ)
    {
        if (k > 3)
        {
            Ex[i][j][k] += tEPY * (Hz[i][j][k] - Hz[i][j - 1][k]) - tEPZ * (Hy[i][j][k] - Hy[i][j][k - 1]);
        }
        else if (k == 3)
        {
            Ex[i][j][k] += tEBY * (Hz[i][j][k] - Hz[i][j - 1][k]) - tEBZ * (Hy[i][j][k] - Hy[i][j][k - 1]);
        }
        else
        {
            Ex[i][j][k] += tERY * (Hz[i][j][k] - Hz[i][j - 1][k]) - tERZ * (Hy[i][j][k] - Hy[i][j][k - 1]);
        }
    }
}

__global__ void UpdateExSource(void)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int k = blockIdx.y * blockDim.y + threadIdx.y;
    int j = 0;
    if (i < LIMX && 1 <= k && k < LIMZ)
    {
        if (k > 3)
        {
            Ex[i][j][k] += tEPY * 2. * Hz[i][j][k] - tEPZ * (Hy[i][j][k] - Hy[i][j][k - 1]);
        }
        else if (k == 3)
        {
            Ex[i][j][k] += tEBY * 2. * Hz[i][j][k] - tEBZ * (Hy[i][j][k] - Hy[i][j][k - 1]);
        }
        else
        {
            Ex[i][j][k] += tERY * 2. * Hz[i][j][k] - tERZ * (Hy[i][j][k] - Hy[i][j][k - 1]);
        }
    }
}
__global__ void UpdateEy(void)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (1 <= i && i < LIMX && j < LIMY - 1 && 1 <= k && k < LIMZ)
    {
        if (k > 3)
        {
            Ey[i][j][k] += tEPZ * (Hx[i][j][k] - Hx[i][j][k - 1]) - tEPX * (Hz[i][j][k] - Hz[i - 1][j][k]);
        }
        else if (k == 3)
        {
            Ey[i][j][k] += tEBZ * (Hx[i][j][k] - Hx[i][j][k - 1]) - tEBX * (Hz[i][j][k] - Hz[i - 1][j][k]);
        }
        else
        {
            Ey[i][j][k] += tERZ * (Hx[i][j][k] - Hx[i][j][k - 1]) - tERX * (Hz[i][j][k] - Hz[i - 1][j][k]);
        }
    }
}

__global__ void UpdateEz(double *d_EzOut)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (1 <= i && i < LIMX && 1 <= j && j < LIMY - 1 && k < LIMZ)
    {
        if (k > 2)
        {
            Ez[i][j][k] += tEPX * (Hy[i][j][k] - Hy[i - 1][j][k]) - tEPY * (Hx[i][j][k] - Hx[i][j - 1][k]);
        }
        else if (k < 3)
        {
            Ez[i][j][k] += tERX * (Hy[i][j][k] - Hy[i - 1][j][k]) - tERY * (Hx[i][j][k] - Hx[i][j - 1][k]);
        }
    }
    __syncthreads();
    d_EzOut[ntime] = Ez[22][40][2] + Ez[22][40][1] + Ez[22][40][0];
}

__global__ void UpdateEzSource(double *d_EzOut)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = 0;
    int k = blockIdx.y * blockDim.y + threadIdx.y;
    if (1 <= i && i < LIMX && k < LIMZ)
    {
        if (k >= 3)
        {
            Ez[i][j][k] += tEPX * (Hy[i][j][k] - Hy[i - 1][j][k]) - tEPY * 2. * Hx[i][j][k];
        }
        else if (k < 3)
        {
            Ez[i][j][k] += tERX * (Hy[i][j][k] - Hy[i - 1][j][k]) - tERY * 2. * Hx[i][j][k];
        }
    }
    __syncthreads();
    d_EzOut[ntime] = Ez[22][40][2] + Ez[22][40][1] + Ez[22][40][0];
}
void UpdateEfields()
{
    // 2 x 100 x 4 = 800 <= 1024
    dim3 threadsPerBlock3D(2, LIMY, 4);
    // 30 x 1 x 4
    dim3 numBlocks3D((LIMX + threadsPerBlock3D.x - 1) / threadsPerBlock3D.x, (LIMY + threadsPerBlock3D.y - 1) / threadsPerBlock3D.y, (LIMZ + threadsPerBlock3D.z - 1) / threadsPerBlock3D.z);
    // 60 x 16 = 960 <= 1024
    dim3 threadPerBlock2D(LIMX, LIMZ);
    // 1 x 1
    dim3 numBlocks2D((LIMX + threadPerBlock2D.x - 1) / threadPerBlock2D.x, (LIMZ + threadPerBlock2D.y - 1) / threadPerBlock2D.y);
    /* Update Electric Fields */
    UpdateEx<<<numBlocks3D, threadsPerBlock3D>>>();
    if (h_ntime < SWITCH1)
    {
        UpdateExSource<<<numBlocks2D, threadPerBlock2D>>>();
    }
    UpdateEy<<<numBlocks3D, threadsPerBlock3D>>>();
    UpdateEz<<<numBlocks3D, threadsPerBlock3D>>>(d_EzOut);
    if (h_ntime < SWITCH1)
    {
        UpdateEzSource<<<numBlocks2D, threadPerBlock2D>>>(d_EzOut);
    }
}
/* End UpdateEfields function ********************************/

/*  Function:  Conductors ()   */
/******************************
/*  Zeros the tangential (Ex, Ey) fields on the conductor
    surfaces (ground plane, microstrip, antenna)
/**************************************/

__global__ void GroundPlane(void)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i < LIMX && j < LIMY)
    {
        Ex[i][j][0] = 0.;
        Ey[i][j][0] = 0.;
    }
}

__global__ void uStrip(void)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = 3;
#ifdef STRIP
    if (19 <= i && i <= 25 && j < LIMY)
    {
        if (i < 25)
            Ex[i][j][k] = 0.;
        Ey[i][j][k] = 0.; // i = 25
    }
#else
    if (19 <= i && i <= 25 && j < 50)
    {
        if (i < 25)
            Ex[i][j][k] = 0.;
        Ey[i][j][k] = 0.; // i = 25
    }
#endif
}

__global__ void PatchAntenna(void)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = 3;
    if (14 <= i && i < 46 && 50 <= j && j < 89)
    {
        Ex[i][j][k] = 0.;
        Ey[i][j][k] = 0.;
    }
    if (i == 46 && 50 <= j && j < 89)
    {
        Ey[i][j][k] = 0.;
    }
    if (j == 89 && 14 <= i && i <= 46)
    {
        Ex[i][j][k] = 0.;
    }
}
void Conductors()
{
    // 10 x 100 = 1000 <= 1024
    dim3 threadPerBlock2D(10, LIMY);
    // 6 x 1
    dim3 numBlocks2D((LIMX + threadPerBlock2D.x - 1) / threadPerBlock2D.x, (LIMY + threadPerBlock2D.y - 1) / threadPerBlock2D.y);
    GroundPlane<<<numBlocks2D, threadPerBlock2D>>>();
    uStrip<<<numBlocks2D, threadPerBlock2D>>>();
#ifndef STRIP
    PatchAntenna<<<numBlocks2D, threadPerBlock2D>>>();
#endif
}

/* End function:  Conductors **************************/

/* Function:  Source ********************************/
/*
/*  Adds in the source *******************************/
__global__ void GaussianSource(void)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x; // from 0 to 6 offset by 19
    int k = blockIdx.y * blockDim.y + threadIdx.y; // from 0 to 2 no offset
    int j = 0;

    temp = (ntime * delT - T0) / T;
    Ez[i + 19][j][k] = exp(-temp * temp);
}

__global__ void NoSource(void)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x; // from 0 to 5 offset by 19
    int k = blockIdx.y * blockDim.y + threadIdx.y; // from 0 to 2 no offset
    int j = 0;

    Ex[i + 19][j][k] = 0.;
}

void Source()
{
    dim3 threadPerBlock(7, 3);
    dim3 threadPerBlock2(6, 3);
    if (h_ntime < SWITCH1)
    {
        GaussianSource<<<1, threadPerBlock>>>(); // 7x3
    }
    else
    {
        NoSource<<<1, threadPerBlock2>>>(); // 6x3
    }

    fprintf(outGauss, "%lf\n", exp(-((h_ntime * h_delT - h_T0) / h_T) * ((h_ntime * h_delT - h_T0) / h_T)));
}

/* End Function:   Source **********************************/

/* Function:  FirstABC() **********************************/
/* ************************************************       */
/* This first ABC is the only one applied to the E-fields.*/
/* Implementation details are in Scheen, 1990.  Performed */
/* after the source is turned off.  Also stores fields    */
/* needed for next round.                                 */

__global__ void ABCY0(void)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int k = blockIdx.y * blockDim.y + threadIdx.y;
    int j = 0;

    if (i < LIMX && k < LIMZ)
    {
        if (k > 3)
        {
            Ex[i][j][k] = ExABC6[i][k] + abcFSy * (Ex[i][j + 1][k] - Ex[i][j][k]);
            Ez[i][j][k] = EzABC6[i][k] + abcFSy * (Ez[i][j + 1][k] - Ez[i][j][k]);
        }
        else if (k == 3)
        {
            Ex[i][j][k] = ExABC6[i][k] + abcBy * (Ex[i][j + 1][k] - Ex[i][j][k]);
            Ez[i][j][k] = EzABC6[i][k] + abcFSy * (Ez[i][j + 1][k] - Ez[i][j][k]);
        }
        else
        {
            Ex[i][j][k] = ExABC6[i][k] + abcDIy * (Ex[i][j + 1][k] - Ex[i][j][k]);
            Ez[i][j][k] = EzABC6[i][k] + abcDIy * (Ez[i][j + 1][k] - Ez[i][j][k]);
        }
    }
}

__global__ void StoreFieldsY0(void)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int k = blockIdx.y * blockDim.y + threadIdx.y;
    int j = 0;

    if (i < LIMX && k < LIMZ)
    {
        ExABC6[i][k] = Ex[i][j + 1][k];
        EzABC6[i][k] = Ez[i][j + 1][k];
    }
}

__global__ void ABCYLIMY(void)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int k = blockIdx.y * blockDim.y + threadIdx.y;
    int j = LIMY - 1;

    if (i < LIMX && k < LIMZ)
    {
        if (k > 3)
        {
            Ex[i][j][k] = ExABC5[i][k] + abcFSy * (Ex[i][j - 1][k] - Ex[i][j][k]);
            Ez[i][j][k] = EzABC5[i][k] + abcFSy * (Ez[i][j - 1][k] - Ez[i][j][k]);
        }
        else if (k == 3)
        {
            Ex[i][j][k] = ExABC5[i][k] + abcBy * (Ex[i][j - 1][k] - Ex[i][j][k]);
            Ez[i][j][k] = EzABC5[i][k] + abcFSy * (Ez[i][j - 1][k] - Ez[i][j][k]);
        }
        else
        {
            Ex[i][j][k] = ExABC5[i][k] + abcDIy * (Ex[i][j - 1][k] - Ex[i][j][k]);
            Ez[i][j][k] = EzABC5[i][k] + abcDIy * (Ez[i][j - 1][k] - Ez[i][j][k]);
        }
    }
}

__global__ void StoreFieldsYLIMY(void)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int k = blockIdx.y * blockDim.y + threadIdx.y;
    int j = LIMY - 1;

    if (i < LIMX && k < LIMZ)
    {
        ExABC5[i][k] = Ex[i][j - 1][k];
        EzABC5[i][k] = Ez[i][j - 1][k];
    }
}

void FirstABC()
{
    // 60 x 16 = 960 <= 1024
    dim3 threadPerBlock2D(LIMX, LIMZ);
    // 1 x 1
    dim3 numBlocks2D((LIMX + threadPerBlock2D.x - 1) / threadPerBlock2D.x, (LIMZ + threadPerBlock2D.y - 1) / threadPerBlock2D.y);
    /* ABC on the wall y=0 */
    if (h_ntime >= SWITCH1 + DELAY)
    {
        ABCY0<<<numBlocks2D, threadPerBlock2D>>>();
    }

    StoreFieldsY0<<<numBlocks2D, threadPerBlock2D>>>();
    ABCYLIMY<<<numBlocks2D, threadPerBlock2D>>>();
    StoreFieldsYLIMY<<<numBlocks2D, threadPerBlock2D>>>();
}
/* End Function:   FirstABC *******************************/

/* Function:  UpdateHfields() *****************************/
/* Updates H-fields.   Nothing special here.  *************/
__global__ void UpdateHx(void)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (i < LIMX && j < LIMY - 1 && k < LIMZ - 1)
    {
        Hx[i][j][k] += tMUZ * (Ey[i][j][k + 1] - Ey[i][j][k]) - tMUY * (Ez[i][j + 1][k] - Ez[i][j][k]);
    }
}

__global__ void UpdateHy(void)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (i < LIMX - 1 && j < LIMY && k < LIMZ - 1)
    {
        Hy[i][j][k] += tMUX * (Ez[i + 1][j][k] - Ez[i][j][k]) - tMUZ * (Ex[i][j][k + 1] - Ex[i][j][k]);
    }
}

__global__ void UpdateHz(void)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (i < LIMX - 1 && j < LIMY - 1 && k < LIMZ)
    {
        Hz[i][j][k] += tMUY * (Ex[i][j + 1][k] - Ex[i][j][k]) - tMUX * (Ey[i + 1][j][k] - Ey[i][j][k]);
    }
}

void UpdateHfields()
{
    // 2 x 100 x 4 = 800 <= 1024
    dim3 threadsPerBlock3D(2, LIMY, 4);
    // 30 x 1 x 4
    dim3 numBlocks3D((LIMX + threadsPerBlock3D.x - 1) / threadsPerBlock3D.x, (LIMY + threadsPerBlock3D.y - 1) / threadsPerBlock3D.y, (LIMZ + threadsPerBlock3D.z - 1) / threadsPerBlock3D.z);

    UpdateHx<<<numBlocks3D, threadsPerBlock3D>>>();
    UpdateHy<<<numBlocks3D, threadsPerBlock3D>>>();
    UpdateHz<<<numBlocks3D, threadsPerBlock3D>>>();
}
/* End Function:   UpdateHfields() ***********************/

/* Function:  SecondABC() *********************************/
/* Implements the remaining ABCs on the walls X = 0, LIMX */
/* and Y = LIMY, Z = LIMZ.   Also, the required fields are*/
/* then stored.                                           */
__global__ void SecondABC_1()
{
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    int k = blockIdx.y * blockDim.y + threadIdx.y;
    int i = 0;

    if (j < LIMY && k < LIMZ)
    {
        if (k > 3)
        {
            Hy[i][j][k] = HyABC2[j][k] + abcFSx * (Hy[i + 1][j][k] - Hy[i][j][k]);
            Hz[i][j][k] = HzABC2[j][k] + abcFSx * (Hz[i + 1][j][k] - Hz[i][j][k]);
        }
        else if (k == 3)
        {
            Hy[i][j][k] = HyABC2[j][k] + abcFSx * (Hy[i + 1][j][k] - Hy[i][j][k]);
            Hz[i][j][k] = HzABC2[j][k] + abcBx * (Hz[i + 1][j][k] - Hz[i][j][k]);
        }
        else
        {
            Hy[i][j][k] = HyABC2[j][k] + abcDIx * (Hy[i + 1][j][k] - Hy[i][j][k]);
            Hz[i][j][k] = HzABC2[j][k] + abcDIx * (Hz[i + 1][j][k] - Hz[i][j][k]);
        }
    }
}

__global__ void SecondABC_2()
{
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    int k = blockIdx.y * blockDim.y + threadIdx.y;
    int i = LIMX - 1;

    if (j < LIMY && k < LIMZ)
    {
        if (k > 3)
        {
            Hy[i][j][k] = HyABC3[j][k] + abcFSx * (Hy[i - 1][j][k] - Hy[i][j][k]);
            Hz[i][j][k] = HzABC3[j][k] + abcFSx * (Hz[i - 1][j][k] - Hz[i][j][k]);
        }
        else if (k == 3)
        {
            Hy[i][j][k] = HyABC3[j][k] + abcFSx * (Hy[i - 1][j][k] - Hy[i][j][k]);
            Hz[i][j][k] = HzABC3[j][k] + abcBx * (Hz[i - 1][j][k] - Hz[i][j][k]);
        }
        else
        {
            Hy[i][j][k] = HyABC3[j][k] + abcDIx * (Hy[i - 1][j][k] - Hy[i][j][k]);
            Hz[i][j][k] = HzABC3[j][k] + abcDIx * (Hz[i - 1][j][k] - Hz[i][j][k]);
        }
    }
}

__global__ void SecondABC_3()
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = LIMZ - 1;

    if (i < LIMX && j < LIMY)
    {
        if (k > 3)
        {
            Hx[i][j][k] = HxABC5[i][j] + abcFSz * (Hx[i][j][k - 1] - Hx[i][j][k]);
            Hy[i][j][k] = HyABC5[i][j] + abcFSz * (Hy[i][j][k - 1] - Hy[i][j][k]);
        }
        else if (k == 3)
        {
            Hx[i][j][k] = HxABC5[i][j] + abcFSz * (Hx[i][j][k - 1] - Hx[i][j][k]);
            Hy[i][j][k] = HyABC5[i][j] + abcFSz * (Hy[i][j][k - 1] - Hy[i][j][k]);
        }
        else
        {
            Hx[i][j][k] = HxABC5[i][j] + abcDIz * (Hx[i][j][k - 1] - Hx[i][j][k]);
            Hy[i][j][k] = HyABC5[i][j] + abcDIz * (Hy[i][j][k - 1] - Hy[i][j][k]);
        }
    }
}

__global__ void SaveFields_1()
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int k = blockIdx.y * blockDim.y + threadIdx.y;
    int j = 0;

    if (i < LIMX && k < LIMZ)
    {
        HxABC4[i][k] = Hx[i][j + 1][k];
        HzABC4[i][k] = Hz[i][j + 1][k];
    }
}

__global__ void SaveFields_2()
{
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    int k = blockIdx.y * blockDim.y + threadIdx.y;
    int i = 0;

    if (j < LIMY && k < LIMZ)
    {
        HyABC2[j][k] = Hy[i + 1][j][k];
        HzABC2[j][k] = Hz[i + 1][j][k];
    }
}

__global__ void SaveFields_3()
{
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    int k = blockIdx.y * blockDim.y + threadIdx.y;
    int i = LIMX - 1;

    if (j < LIMY && k < LIMZ)
    {
        HyABC3[j][k] = Hy[i - 1][j][k];
        HzABC3[j][k] = Hz[i - 1][j][k];
    }
}

__global__ void SaveFields_4()
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int k = blockIdx.y * blockDim.y + threadIdx.y;
    int j = LIMY - 1;

    if (i < LIMX && k < LIMZ)
    {
        HxABC1[i][k] = Hx[i][j - 1][k];
        HzABC1[i][k] = Hz[i][j - 1][k];
    }
}

__global__ void SaveFields_5()
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = LIMZ - 1;

    if (i < LIMX && j < LIMY)
    {
        HxABC5[i][j] = Hx[i][j][k - 1];
        HyABC5[i][j] = Hy[i][j][k - 1];
    }
}

void SecondABC()
{
    dim3 threadsPerBlock(4, 4);
    dim3 numBlocks1((LIMY + threadsPerBlock.x - 1) / threadsPerBlock.x, (LIMZ + threadsPerBlock.y - 1) / threadsPerBlock.y);
    dim3 numBlocks2((LIMX + threadsPerBlock.x - 1) / threadsPerBlock.x, (LIMY + threadsPerBlock.y - 1) / threadsPerBlock.y);
    dim3 numBlocks3((LIMX + threadsPerBlock.x - 1) / threadsPerBlock.x, (LIMZ + threadsPerBlock.y - 1) / threadsPerBlock.y);

    SecondABC_1<<<numBlocks1, threadsPerBlock>>>();
    SecondABC_2<<<numBlocks1, threadsPerBlock>>>();
    SecondABC_3<<<numBlocks2, threadsPerBlock>>>();

    SaveFields_1<<<numBlocks3, threadsPerBlock>>>(); // y = 0
    SaveFields_2<<<numBlocks1, threadsPerBlock>>>(); // x = 0
    SaveFields_3<<<numBlocks1, threadsPerBlock>>>(); // x = LIMX - 1
    SaveFields_4<<<numBlocks3, threadsPerBlock>>>(); // y = LIMY - 1
    SaveFields_5<<<numBlocks2, threadsPerBlock>>>(); // z = LIMZ - 1 (top)
}
/* End Function:   SecondABC() *****************************/