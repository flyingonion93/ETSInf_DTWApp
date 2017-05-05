__constant__ MyType CCosts[T_COSTS]; /* Placed within GPU Const Memory. Using by DTW */
__constant__ int  Ckmax_fft[N_MIDI]; /* Placed within GPU Const Memory. Using by FFT */
__constant__ int  Ckmin_fft[N_MIDI]; /* Placed within GPU Const Memory. Using by FFT */

#include "DTWExec_GPU.h"
#include "../common/DTWIO_GPU.h"
#include "../dtw/DTWProcess_GPU.h"
#include "../common/kernels.cuh"

int DTWExec (int argc , char *argv[])
{  
   /* for the GPU check, control and ... */
   int
     maxGrid,    TPBlock,
     GridTTRAMA, GridNBASES,
     GridNFFTd2, GridNFFT,
     GridNMID32, BPGrid5, BPGrid6,
     BPGrid9;

   cudaEvent_t
     start, stop;

   cublasHandle_t
     handle;
            
   /* Using only by DTW */
   int
     *pS=NULL, pos_min,
     DTWWhere, DTWSize,
     DTWSizeMinusOne,
     DTWSizePlusPad;

   MyType
     *pD = NULL,
     *pV = NULL,
     Costs[T_COSTS];


   /* Using only by FFT */
   MyFFTType
     plan;

   MyType
     *X_fft  =NULL,
     *Out_fft=NULL,
     *Mod_fft=NULL;
     
   int
     kmax_fft[N_MIDI],
     kmin_fft[N_MIDI],
     N_FFTdiv2;


   /* Using by other preprocessing steps */
   MyType
     *norms     = NULL,
     *s_fk      = NULL,
     *v_SxD     = NULL,
     *v_cfreq   = NULL,
     *v_hanning = NULL,
     *sdata     = NULL,
     *mdata     = NULL,
     *trama     = NULL,
     *tauxi     = NULL,
     *ts_fk     = NULL,
     *v_dxState = NULL,
     BETA=(MyType)1.5, vnorm;

   int
     *states_time_i = NULL,
     *states_time_e = NULL,
     *states_seq    = NULL,
     *mpos          = NULL;

   DTWconst  Param;
   DTWfiles  NameFiles;

   /* General/others varibles */
   int    i, mode;
   FILE   *fp;
   float  time;

   /* Reading global paramentes */
   if (argc != 5) {
      printf("Usage: %s <BETA> <configuration file> <threadsPerBlock> <mode>. ", argv[0]);
      printf("Example: main 1.5 parametes.dat 64 1, where the meaning of the last parameter is:\n");
      printf("   1: Unified Memory with    CuBLAS\n");
      printf("   2: Unified Memory without CuBLAS for reduction operations inside GPU\n");
      printf("   3: Classic Memory with    CuBLAS\n");
      printf("   4: Classic Memory without CuBLAS for reduction operations inside GPU\n");
      return -1;
   }
   #ifdef SIMPLE
     BETA=strtof(argv[1], NULL);
   #else
     BETA=strtod(argv[1], NULL);
   #endif
   TPBlock=atoi(argv[3]);
   mode   =atoi(argv[4]);

   /* Have we a compatible GPU? We assume that we only use one GPU (with ID = 0) */
   if (!HaveCompatibleGPU(maxGrid)) return -1;
   
   #ifndef SIMPLE
     CUDAERR(cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte));
   #endif

   CUBLASERR(cublasCreate(&handle));
   cudaEventCreate(&start);
   cudaEventCreate(&stop);

   /* Reading general information and file names */
   if (!ReadParameters(&Param, &NameFiles, argv[2])) return -1;

   /* Allocating memory and reading some structures */
   if (!AllocUnifiedData(&v_hanning, &states_time_i, &states_time_e, &states_seq, TAMTRAMA, Param.N_STATES, NameFiles))
      return -1;

   /* Allocating memory for s_fk or auxiliar structures when Beta !=0.0 and !=1.0 */
   if (!AllocUnifiedS_fk(&s_fk, &tauxi, &ts_fk, BETA, N_MIDI_PAD, Param.N_BASES, NameFiles)) return -1;

   /* Allocating memory for FFT memory and reading data */
   if (!AllocUnifiedFFT(&plan, &X_fft, &Out_fft, &Mod_fft, kmin_fft, kmax_fft, N_FFT, NameFiles)) return -1;

   /* Initializing vector of costs and const memory GPU with this information */
   Costs[0]=(MyType)1.00; Costs[1]=(MyType)1.05; Costs[2]=(MyType)1.10; Costs[3]=(MyType)1.20;
   Costs[4]=(MyType)1.05; Costs[5]=(MyType)1.10; Costs[6]=(MyType)1.20;

   /* Moving info to GPU const memory */
   CUDAERR(cudaMemcpyToSymbol(Ckmax_fft, kmax_fft, sizeof(kmax_fft)));
   CUDAERR(cudaMemcpyToSymbol(Ckmin_fft, kmin_fft, sizeof(kmin_fft)));
   CUDAERR(cudaMemcpyToSymbol(CCosts,    Costs,    sizeof(Costs)));

   /* Allocating memory for and setting some DTW constants */
   DTWSize        = states_time_e[Param.N_STATES - 1] + 1;
   DTWSizeMinusOne= DTWSize - 1;
   DTWSizePlusPad = (DTWSizeMinusOne +  N_COSTS) * (N_COSTS + 1);
   if (!AllocUnifiedDTW(&pD, &pV, &pS, &v_SxD, DTWSize, DTWSizePlusPad)) return -1;

   /* Allocating memory for the rest of general structures */
   if (!AllocUnifiedAuxi(&norms, &trama, &v_cfreq, &v_dxState, Param.N_BASES, TAMTRAMA, N_MIDI+1)) return -1;
 
   /* Allocating memory for GPU strutuctures related with shared memory */
   if (!AllocUnifiedShar(&sdata, &mdata, &mpos, maxGrid, DTWSize, DTWSizeMinusOne)) return -1;

   /* Some helper constants */
   N_FFTdiv2 = N_FFT/2;
   GridNBASES = (Param.N_BASES + TPBlock - 1) / TPBlock;
   GridTTRAMA = (TAMTRAMA      + TPBlock - 1) / TPBlock;
   GridNFFT   = (N_FFT         + TPBlock - 1) / TPBlock;   
   GridNFFTd2 = (N_FFTdiv2     + TPBlock)     / TPBlock;
   GridNMID32 = (N_MIDI        + 31)          / 32;
   #ifdef TESLA
     dim3 GridCompDi((Param.N_BASES + 31) / 32, 1);
     dim3 TPBCompDi(32,32);
   #else
     dim3 GridCompDi((Param.N_BASES + 15) / 16, 1);
     dim3 TPBCompDi(32,16);
   #endif

   BPGrid5 = (Param.N_STATES + TPBlock - 1) / TPBlock;
   BPGrid6 = (DTWSize        + TPBlock - 1) / TPBlock;
   BPGrid9 = (DTWSizeMinusOne+ TPBlock - 1) / TPBlock;

   // Init DTW, auxiliar structures and compute norms */
   kernel_InitDTW<<<(DTWSizePlusPad+TPBlock-1)/TPBlock, TPBlock>>>(pD, pV, pS, DTWSizePlusPad-DTWSize+1, DTWSizePlusPad);
   if (BETA>=0.0 && BETA<=0.0)
      kernel_CompNorB0<<<GridNBASES, TPBlock>>>(norms, (MyType)N_MIDI, Param.N_BASES);
   else if (BETA>=1.0 && BETA<=1.0)
      kernel_CompNorB1<<<GridNBASES, TPBlock>>>(norms, s_fk, N_MIDI, Param.N_BASES);
   else {
      kernel_CompNorBG<<<GridNBASES, TPBlock>>>(norms, s_fk, N_MIDI_PAD, BETA, Param.N_BASES);
      kernel_PowToReal<<<(N_MIDI_PAD*Param.N_BASES+TPBlock-1)/TPBlock, TPBlock>>>(ts_fk, s_fk, BETA-1.0, N_MIDI_PAD*Param.N_BASES);
   }   
   
   /* Opening file with frames */
   fp=fopen(NameFiles.file_trama, "rb");
   if(fp==NULL)
   {
      printf("Error opening %s input file\n", NameFiles.file_trama);
      return -1;
   }

   switch(mode) {
   case UnifiedCuBLAS:
     cudaDeviceSynchronize();
     cudaEventRecord(start, 0);
     for(i=0; i<Param.N_TRAMAS; i++)
     {
       // Read new frame (trama) 
       ReadFrame(trama, TAMTRAMA, fp);

       // Apply Hanning
       kernel_ApplyWin<<<GridTTRAMA, TPBlock>>>(trama, v_hanning, TAMTRAMA);

       // Compute FFT
       kernel_CopyAndSet<<<GridNFFT, TPBlock>>>(X_fft, trama, TAMTRAMA, N_FFT);
       FFTGPU(X_fft, Out_fft, &plan);
       kernel_Modul<<<GridNFFTd2, TPBlock>>>(Mod_fft, Out_fft, N_FFTdiv2);
       kernel_Cfreq<<<N_MIDI, 32>>>(v_cfreq, Mod_fft);

       // Compute distortion
       if (BETA>=0.0 && BETA<=0.0) {
         kernel_CompDisB0<<<GridNBASES, TPBlock>>>(v_dxState, v_cfreq, norms, s_fk, N_MIDI, Param.N_BASES);
         //kernel_CompDisB0<<<GridCompDi, TPBCompDi>>>(v_dxState, v_cfreq, norms, s_fk, N_MIDI, Param.N_BASES);
       }
       else if (BETA>=1.0 && BETA<=1.0) {
         kernel_Reduction<<<1, sizeWarp>>>(v_cfreq, N_MIDI);
         kernel_CompDisB1<<<GridNBASES, TPBlock>>>(v_dxState, v_cfreq, norms, s_fk, N_MIDI, Param.N_BASES);
       }
       else {
         kernel_PowToReal<<<GridNMID32, 32>>>(tauxi, v_cfreq, BETA, N_MIDI);
         kernel_CompDisBG<<<GridCompDi, TPBCompDi>>>(v_dxState, v_cfreq, norms, s_fk, ts_fk, tauxi, BETA, N_MIDI, Param.N_BASES);
       }

       // Apply distortion
       if (mode == UnifiedCuBLAS) {
         kernel_CalSxCuBlas1<<<BPGrid5, TPBlock>>>(v_SxD, v_dxState, states_time_i, states_time_e,states_seq, Param.N_STATES);
         #ifdef SIMPLE
           CUBLASERR(cublasSasum(handle, DTWSize, v_SxD, 1, &vnorm));
           vnorm = 1.0f / (sqrtf(vnorm) + FLT_EPSILON);
         #else
           CUBLASERR(cublasDasum(handle, DTWSize, v_SxD, 1, &vnorm));
           vnorm = 1.0  / (sqrt(vnorm)  + DBL_EPSILON);
         #endif                                                                                  
         kernel_CalSxCuBlas2<<<BPGrid6, TPBlock>>>(v_SxD, Param.ALPHA, vnorm, DTWSize);
       }
       else {
         kernel_CalSxD1<<<BPGrid5, TPBlock>>>(v_SxD, v_dxState, states_time_i, states_time_e, states_seq, Param.N_STATES);
         MySumPows(sdata, v_SxD, maxGrid, DTWSize);
         kernel_CalSxD2<<<BPGrid6, TPBlock>>>(v_SxD, Param.ALPHA, sdata, DTWSize);
       }

       // DTW      
       DTWWhere=(i % TBLOCK) * (N_COSTS + DTWSizeMinusOne) + N_COSTS;
       if (mode == UnifiedCuBLAS) {
         #ifdef SIMPLE
           CUBLASERR(cublasIsamin(handle, DTWSizeMinusOne, &v_SxD[1], 1, &pos_min));
         #else
           CUBLASERR(cublasIdamin(handle, DTWSizeMinusOne, &v_SxD[1], 1, &pos_min));
         #endif
         
         if (pos_min != 1) {
           kernel_DTW<<<BPGrid9, TPBlock>>>(&v_SxD[1], pD, pV, pS, i, DTWWhere, DTWSizeMinusOne);

           #ifdef SIMPLE
             CUBLASERR(cublasIsamin(handle, DTWSizeMinusOne, &pV[DTWWhere], 1, &pos_min));
           #else
             CUBLASERR(cublasIdamin(handle, DTWSizeMinusOne, &pV[DTWWhere], 1, &pos_min));
           #endif
         }
         pos_min--;
       }
       else {
         pos_min = MyImin(mdata, mpos, &v_SxD[1], maxGrid, DTWSizeMinusOne);
         if (pos_min !=0) {
           kernel_DTW<<<BPGrid9, TPBlock>>>(&v_SxD[1], pD, pV, pS, i, DTWWhere, DTWSizeMinusOne);
           pos_min = MyImin(mdata, mpos, &pV[DTWWhere], maxGrid, DTWSizeMinusOne);
         }
       }
       // END. The result is:
       printf("%d,%d\n", i, pos_min);
     }
     cudaEventRecord(stop, 0);
     cudaEventSynchronize(stop);
     cudaEventElapsedTime(&time, start, stop);
     break;

   case ClassicCuBLAS:
     break;
   default:
     printf("Default case does not defined\n");
   }
   printf("%f msec.\n", time);   



   /* CLOSE and Free ALL */
   fclose(fp);
   FreeFiles(&NameFiles);

   cudaEventDestroy(start);
   cudaEventDestroy(stop);

   CUBLASERR(cublasDestroy(handle));

   CUFFTERR(cufftDestroy(plan));

   if (!(BETA>=(MyType)0.0 && BETA<=(MyType)0.0) && !(BETA>=(MyType)1.0 && BETA<=(MyType)1.0)) {
      CUDAERR(cudaFree(tauxi));
      CUDAERR(cudaFree(ts_fk));
   }      

   CUDAERR(cudaFree(mdata));
   CUDAERR(cudaFree(mpos));
   CUDAERR(cudaFree(Mod_fft));
   CUDAERR(cudaFree(norms));
   CUDAERR(cudaFree(Out_fft));
   CUDAERR(cudaFree(pD));
   CUDAERR(cudaFree(pS));
   CUDAERR(cudaFree(pV));
   CUDAERR(cudaFree(sdata));
   CUDAERR(cudaFree(states_seq));
   CUDAERR(cudaFree(states_time_i));
   CUDAERR(cudaFree(states_time_e));
   CUDAERR(cudaFree(v_cfreq));
   CUDAERR(cudaFree(v_dxState));
   CUDAERR(cudaFree(v_hanning));
   CUDAERR(cudaFree(s_fk));
   CUDAERR(cudaFree(trama));
   CUDAERR(cudaFree(v_SxD));
   CUDAERR(cudaFree(X_fft));

   return 0;
}
