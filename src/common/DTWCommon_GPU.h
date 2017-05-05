#pragma once

#ifndef DTW_Common_GPU_h
	#define DTW_Common_GPU_h

#define N_MIDI     114 /* For now, N_MIDI     = 114             */
#define N_MIDI_PAD 128 /* Para que esté alineado                */
#define N_FFT    16384 /* For now, N_FFT      = 16384           */
#define TAMTRAMA  5700 /* For now, TAMTRAMA   = 10*TAMMUESTRA   */
#define TAMMUESTRA 570 /* For now, TAMMUESTRA = 570             */

#define N_COSTS            4  /* For now 2 dimensins and 4 costs per dimension    */
#define T_COSTS (2*N_COSTS-1) /* Thereby total number of costs used is this value */
#define TBLOCK    (N_COSTS+1) /* another used derivate value                      */

/* GPU constants */
#define sizeWarp    32 
#define maxThreads 512
#define maxBlocks  256

/* For execution */
#define UnifiedCuBLAS   1
#define UnifiedNoCuBLAS 2
#define ClassicCuBLAS   3
#define ClassicNoCuBLAS 4

#ifdef SIMPLE
   #define MyType float
#else
   #define MyType double
#endif

#define MyFFTType cufftHandle

#ifndef min
   #define min(x,y) ((x < y) ? x : y)
#endif
   

/* Members of DTWparam and DTWfiles are the Preprocess and DTW */
/* parameters. Each composition needs a file with values for   */
/* these parameters. The structure of this text file should be */
/* the same as this DTWparam and DTWfiles structs, one line    */
/* per parameter. No empty lines and one empty line to the end */
/* Example of this file:         */
/* 57                            */
/* 77                            */
/* 1                             */
/* 11.5129                       */
/* Datos/30sec/Input/hanning.bin */
/* .....                         */
/* Datos/30sec/Verif/v_SxD.bin   */
/*                               */
/* BETA is not in these files. BETA is a program argument      */
struct DTWconst
{
   int N_BASES;    /* no default value                         */
   int N_STATES;   /* no default value                         */
   int N_TRAMAS;   /* no default value                         */
   MyType ALPHA;   /* default value is 11.5129                 */
};


struct DTWfiles
{
   /* Files with input data for the problem */
   char *file_hanning;
   char *file_trama;
   char *file_partitura;
   char *file_kmax;
   char *file_kmin;
   char *fileStates_Time_e;
   char *fileStates_Time_i;
   char *fileStates_seq;

   /* Files with data for verification */
   char *FileVerifyHanning;
   char *FileVerifyFFT;
   char *FileVerifyVecDist;
   char *FileVerifyDistor;
};

typedef struct DTWconst  DTWconst;
typedef struct DTWfiles  DTWfiles;

#endif
