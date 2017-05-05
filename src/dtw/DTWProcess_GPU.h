#pragma once

#ifndef DTW_Process_GPU_h
	#define DTW_Process_GPU_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* CUDA-C includes */
#include <cuda_runtime.h>

/* CuFFT includes */
#include <cufft.h>

/* CuBLAS includes */
#include <cublas_v2.h>

#include "DTWCommon_GPU.h"

#define CUDAERR(x) do { if((x)!=cudaSuccess) { \
   printf("CUDA error: %s : %s, line %d\n", cudaGetErrorString(x), __FILE__, __LINE__);\
   return EXIT_FAILURE;}} while(0)

#define CUBLASERR(x) do { if((x)!=CUBLAS_STATUS_SUCCESS) { \
   printf("CUBLAS error: %s, line %d\n", __FILE__, __LINE__);\
   return EXIT_FAILURE;}} while(0)

#define CUFFTERR(x) do { if((x)!=CUFFT_SUCCESS) { \
   printf("CUFFT error: %s, line %d\n", __FILE__, __LINE__);\
   return EXIT_FAILURE;}} while(0)

                                                       
/* ******************************** Preproceso Functions Prototypes  **************************** */
void         BlocksAndThreads(int*, int*, int*, const int, const int);
int          FFTGPU(MyType*, MyType*, MyFFTType*);
bool         HaveCompatibleGPU(int &);
inline bool  IsPow2(unsigned int);
int          MyImin(MyType *, int*, MyType *, const int, const int);
void         MySumPows(MyType *, MyType *, const int, const int);
unsigned int NextPow2(unsigned int);

bool AllocUnifiedAuxi(MyType **,  MyType **, MyType **, MyType **,    const int,  const int, const int);
bool AllocUnifiedData(MyType **,  int **,    int **,    int **,       const int,  const int, DTWfiles);
bool AllocUnifiedDTW (MyType **,  MyType **, int **,    MyType **,    const int,  const int);
bool AllocUnifiedFFT(MyFFTType *, MyType **, MyType **, MyType **,    int*, int*, const int, DTWfiles);
bool AllocUnifiedShar(MyType **,  MyType **, int **,    const int,    const int,  const int);
bool AllocUnifiedS_fk(MyType **,  MyType **, MyType **, const MyType, const int,  const int, DTWfiles);


/* ********************************* Preproceso kernels Prototypes  ***************************** */
__global__ void kernel_InitDTW(MyType* __restrict__, MyType* __restrict__, int* __restrict__, const int, const int);

__global__ void kernel_CompNorB0(MyType* __restrict__, const MyType,               const int);
__global__ void kernel_CompNorB1(MyType* __restrict__, const MyType* __restrict__, const int,    const int);
__global__ void kernel_CompNorBG(MyType* __restrict__, const MyType* __restrict__, const int,    const MyType, const int);
__global__ void kernel_PowToReal(MyType* __restrict__, const MyType* __restrict__, const MyType, const int);

__global__ void kernel_ApplyWin  (MyType* __restrict__, const MyType* __restrict__, const int);
__global__ void kernel_Cfreq     (MyType* __restrict__, const MyType* __restrict__);
__global__ void kernel_CopyAndSet(MyType* __restrict__, const MyType* __restrict__, const int, const int);
__global__ void kernel_Modul     (MyType* __restrict__, const MyType* __restrict__, const int);

__global__ void kernel_Reduction (MyType* __restrict__, const int);

__global__ void kernel_SumPows(MyType* __restrict__, const MyType* __restrict__, const int, const bool, const int);
__global__ void kernel_Sum    (MyType* __restrict__, const MyType* __restrict__, const int, const bool, const int);
__global__ void kernel_Vnorm  (MyType* __restrict__);

__global__ void kernel_MyImin(MyType* __restrict__, int* __restrict__, const MyType* __restrict__,
                              const int, const bool, const int);
                              
__global__ void kernel_MyIminLast(MyType* __restrict__, int* __restrict__, const MyType* __restrict__,
                                  const int* __restrict__, const int, const bool, const int);

__global__ void kernel_CalSxD1(MyType* __restrict__, const MyType*, const int* __restrict__,
                               const int* __restrict__, const int* __restrict__ , const int);                               
__global__ void kernel_CalSxD2(MyType* __restrict__, const MyType,  const MyType* __restrict__, const int);

__global__ void kernel_CalSxCuBlas1(MyType* __restrict__, const MyType* __restrict__, const int* __restrict__,
                                    const int* __restrict__, const int* __restrict__, const int);                                    
__global__ void kernel_CalSxCuBlas2(MyType* __restrict__, const MyType,  const MyType, const int);


__global__ void kernel_DTW(const MyType* __restrict__, MyType* __restrict__, MyType* __restrict__, int* __restrict__,
                           const int, const int, const int);

__global__ void kernel_CompDisB0(MyType* __restrict__, const MyType* __restrict__, const MyType* __restrict__,
                                 const MyType* __restrict__, const int, const int);

__global__ void kernel_CompDisB1(MyType* __restrict__, const MyType* __restrict__, const MyType* __restrict__,
                                 const MyType* __restrict__, const int, const int);

__global__ void kernel_CompDisBG(MyType* __restrict__, const MyType* __restrict__, const MyType* __restrict__,
                                 const MyType* __restrict__, const MyType* __restrict__, const MyType* __restrict__,
                                 const MyType, const int, const int);

#endif
