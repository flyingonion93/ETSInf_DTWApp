#include "funcionesPreprocesado.h"
#include "funcionesFicheros.h"

unsigned int NextPow2(unsigned int x)
{
   --x;
   x |= x >> 1;
   x |= x >> 2;
   x |= x >> 4;
   x |= x >> 8;
   x |= x >> 16;
   return ++x;
}

inline bool IsPow2(unsigned int x) { return ((x&(x-1))==0); }
                            
bool HaveCompatibleGPU(int &maxGrid)
{
   int deviceCount, driverVersion;
   
   cudaDeviceProp deviceProp;

   CUDAERR(cudaGetDeviceCount(&deviceCount));
   if (deviceCount == 0) {
      printf("Your GPU does not support CUDA\n");
      return false;
   }

   CUDAERR(cudaGetDeviceProperties(&deviceProp, 0));
   if (deviceProp.major < 3) {
      printf("Sorry, we need CUDA Capability >=3\n");
      return false;
   }
   else
     maxGrid=deviceProp.maxGridSize[0];
   
   CUDAERR(cudaDriverGetVersion(&driverVersion));
   if ((driverVersion/1000) < 6) {
      printf("Sorry, we need CUDA Version >=6\n");
      return false;
   }

   if (!deviceProp.unifiedAddressing) {
      printf("Your system does not support Unified Memory (32 bit or MW Vista, 7, etc.?)\n");
      return false;
   }

   return true;
}

bool AllocUnifiedS_fk(MyType **s_fk, MyType **tauxi, MyType **ts_fk, const MyType BETA,
                      const int nmidi, const int nbases, DTWfiles NameFiles)
{
   CUDAERR(cudaMallocManaged((void **)s_fk, sizeof(MyType)*nmidi*nbases, cudaMemAttachGlobal));
   if (s_fk==NULL)
   {
      printf("Error allocating memory for vector s_fk.\n");
      return false;
   }
   if (!ReadS_fk((*s_fk), nbases, NameFiles.file_partitura))
   {
      printf("Error reading vector s_fk\n");
      return false;
   }
   
   if (!(BETA>=(MyType)0.0 && BETA<=(MyType)0.0) && !(BETA>=(MyType)1.0 && BETA<=(MyType)1.0))
   {
      CUDAERR(cudaMallocManaged((void **)tauxi, sizeof(MyType)*nmidi, cudaMemAttachGlobal));
      if (tauxi==NULL)
      {
         printf("Error allocating memory for auxiliar vector tauxi.\n");
         return false;
      }
   
      CUDAERR(cudaMallocManaged((void **)ts_fk, sizeof(MyType)*nmidi*nbases, cudaMemAttachGlobal));
      if (ts_fk==NULL)
      {
         printf("Error allocating memory for auxiliar vector ts_fk.\n");
         return false;
      }
   }

   return true;
}

bool AllocUnifiedData(MyType **v_hanning, int **states_time_i,  int **states_time_e, int **states_seq,
                      const int tamtrama, const int nstates, DTWfiles NameFiles)
{
   CUDAERR(cudaMallocManaged((void **)v_hanning,  sizeof(MyType)*tamtrama, cudaMemAttachGlobal));
   if ((*v_hanning) == NULL)
   {
      printf("Error allocating memory for vector v_hanning.\n");
      return false;
   }
   CUDAERR(cudaMallocManaged((void **)states_time_i, sizeof(int)*nstates, cudaMemAttachGlobal));
   if ((*states_time_i) == NULL)
   {
      printf("Error allocating memory for vector states_time_i.\n");
      return false;
   }
   CUDAERR(cudaMallocManaged((void **)states_time_e, sizeof(int)*nstates, cudaMemAttachGlobal));
   if ((*states_time_e) == NULL)
   {
      printf("Error allocating memory for vector states_time_e.\n");
      return false;
   }
   CUDAERR(cudaMallocManaged((void **)states_seq,    sizeof(int)*nstates, cudaMemAttachGlobal));
   if ((*states_seq) == NULL)
   {
      printf("Error allocating memory for vector states_seq.\n");
      return false;
   }

   if (!ReadVector((*v_hanning), tamtrama, NameFiles.file_hanning)){
      printf("Error reading vector v_hanning %s\n", NameFiles.file_hanning);
      return -1;
   }
   if (!ReadVectorInt64((*states_seq), nstates, NameFiles.fileStates_seq)){
      printf("Error reading vector states_seq\n");
      return -1;
   }
   if (!ReadVectorInt64((*states_time_i), nstates, NameFiles.fileStates_Time_i)){
      printf("Error reading vector states_time_i\n");
      return -1;
   }
   if (!ReadVectorInt64((*states_time_e), nstates, NameFiles.fileStates_Time_e)){
      printf("Error reading vector states_time_e\n");
      return -1;
   }

   return true;
}

bool AllocUnifiedFFT(MyFFTType *plan, MyType **X_fft, MyType **Out_fft, MyType **Mod_fft,
                     int *kmin_fft, int *kmax_fft, const int nfft, DTWfiles NameFiles)
{
   CUDAERR(cudaMallocManaged((void **)X_fft, sizeof(MyType)*2*nfft+1, cudaMemAttachGlobal));
   if ((*X_fft )==NULL)
   {
      printf("Error reading vector X_fft\n");
      return false;
   }

   /* Puede ser suficiente con nfft/2+1? */
   CUDAERR(cudaMallocManaged((void **)Mod_fft, sizeof(MyType)*nfft, cudaMemAttachGlobal));
   if ((*Mod_fft )==NULL)
   {
      printf("Error reading vector Mod_fft\n");
      return false;
   }

   #ifdef SIMPLE
      CUDAERR(cudaMallocManaged((void **)Out_fft, sizeof(cufftComplex)*nfft, cudaMemAttachGlobal));
      if ((*Out_fft )==NULL)
      {
         printf("Error reading vector Out_fft\n");
         return false;
      }
      CUFFTERR(cufftPlan1d(plan, nfft, CUFFT_R2C, 1));
   #else
      CUDAERR(cudaMallocManaged((void **)Out_fft, sizeof(cufftDoubleComplex)*nfft, cudaMemAttachGlobal));
      if ((*Out_fft )==NULL)
      {
         printf("Error reading vector Out_fft\n");
         return false;
      }
      CUFFTERR(cufftPlan1d(plan, nfft, CUFFT_D2Z, 1));
   #endif
   if (plan==NULL)
   {
      printf("Error creating FFT plan Out_fft\n");
      return false;
   }

   if (!ReadVectorInt64(kmax_fft, N_MIDI, NameFiles.file_kmax))
   {
      printf("Error reading vector kmax_fft\n");
      return false;
   }
   if (!ReadVectorInt64(kmin_fft, N_MIDI, NameFiles.file_kmin))
   {
      printf("Error reading vector kmin_fft\n");
      return false;
   }

   return true;
}

bool AllocUnifiedDTW(MyType **pD, MyType **pV, int **pS, MyType **v_SxD, const int DTWSize, const int DTWSizePlusPad)
{
   CUDAERR(cudaMallocManaged((void **)pD, sizeof(MyType)*DTWSizePlusPad, cudaMemAttachGlobal));
   if ((*pD)==NULL)
   {
      printf("Error allocating memory for vector pD.\n");
      return false;
   }
   CUDAERR(cudaMallocManaged((void **)pV, sizeof(MyType)*DTWSizePlusPad, cudaMemAttachGlobal));
   if ((*pV)==NULL)
   {
      printf("Error allocating memory for vector pV.\n");
      return false;
   }
   CUDAERR(cudaMallocManaged((void **)pS, sizeof(int)*DTWSizePlusPad, cudaMemAttachGlobal));
   if ((*pS)==NULL)
   {
      printf("Error allocating memory for vector pS.\n");
      return false;
   }
   CUDAERR(cudaMallocManaged((void **)v_SxD, sizeof(MyType)*DTWSize, cudaMemAttachGlobal));
   if ((*pS)==NULL)
   {
      printf("Error allocating memory for vector v_SxD.\n");
      return false;
   }

   return true;
}

bool AllocUnifiedAuxi(MyType **norms, MyType **trama, MyType **v_cfreq, MyType **v_dxState,
                      const int nbases, const int tamtrama, const int nmidiplusone)
{
   CUDAERR(cudaMallocManaged((void **)norms, sizeof(MyType)*nbases, cudaMemAttachGlobal));
   if ((*norms)==NULL)
   {
      printf("Error allocating memory for vector norms.\n");
      return false;
   }
   CUDAERR(cudaMallocManaged((void **)v_dxState, sizeof(MyType)*nbases, cudaMemAttachGlobal));
   if ((*v_dxState)==NULL)
   {
      printf("Error allocating memory for vector v_dxState.\n");
      return false;
   }
   CUDAERR(cudaMallocManaged((void **)trama, sizeof(MyType)*tamtrama, cudaMemAttachGlobal));
   if ((*trama)==NULL)
   {
      printf("Error allocating memory for vector trama.\n");
      return false;
   }
   CUDAERR(cudaMallocManaged((void **)v_cfreq, sizeof(int)*nmidiplusone, cudaMemAttachGlobal));
   if ((*v_cfreq)==NULL)
   {
      printf("Error allocating memory for vector v_cfreq.\n");
      return false;
   }

   return true;
}

bool AllocUnifiedShar(MyType **sdata, MyType **mdata, int **mpos, const int maxGrid, const int DTWSize,
                      const int DTWSizeMinusOne)
{
   int numThreads, numBlocks, sharedSize;

   BlocksAndThreads(&numBlocks, &numThreads, &sharedSize, maxGrid, DTWSize);
   CUDAERR(cudaMallocManaged((void **)sdata, sizeof(MyType)*numBlocks, cudaMemAttachGlobal));
   if ((*sdata)==NULL)
   {
      printf("Error allocating memory for vector sdata.\n");
      return false;
   }

   BlocksAndThreads(&numBlocks, &numThreads, &sharedSize, maxGrid, DTWSizeMinusOne);
   CUDAERR(cudaMallocManaged((void **)mdata, sizeof(MyType)*numBlocks, cudaMemAttachGlobal));
   if ((*mdata)==NULL)
   {
      printf("Error allocating memory for vector mdata.\n");
      return false;
   }
   CUDAERR(cudaMallocManaged((void **)mpos, sizeof(int)*numBlocks, cudaMemAttachGlobal));
   if ((*mpos)==NULL)
   {
      printf("Error allocating memory for vector mpos.\n");
      return false;
   }

   return true;
}

void BlocksAndThreads(int *blocks, int *threads, int *sharedsize, const int maxGrid, const int size)
{
   (*threads) = (size < maxThreads*2) ? NextPow2((size + 1)/ 2) : maxThreads;
   (*blocks)  = (size + ((*threads) * 2 - 1)) / ((*threads) * 2);

   if ((*blocks) > maxGrid)
   {
      (*blocks)  /= 2;
      (*threads) *= 2;
   }

   (*blocks)     = min(maxBlocks, (*blocks));
   (*sharedsize) = ((*threads) <= 32) ? 2*(*threads)*sizeof(MyType) : (*threads)*sizeof(MyType);
}

int FFTGPU(MyType *X_fft, MyType *Out_fft, MyFFTType *plan)
{   
   #ifdef SIMPLE
      CUFFTERR(cufftExecR2C(*plan, (cufftReal *)X_fft,       (cufftComplex *)Out_fft));
   #else
      CUFFTERR(cufftExecD2Z(*plan, (cufftDoubleReal *)X_fft, (cufftDoubleComplex *)Out_fft));
   #endif

   return 0;
}

int MyImin(MyType *odata, int *opos, MyType *idata, const int maxGrid, const int size)
{
   int numBlocks=0, numThreads=0, sharedSize=0, s;

   BlocksAndThreads(&numBlocks, &numThreads, &sharedSize, maxGrid, size);

   kernel_MyImin<<<numBlocks, numThreads, sharedSize*2>>>(odata, opos, idata, numThreads, IsPow2(size), size);

   s = numBlocks;
   while (s > 1)
   {
      BlocksAndThreads(&numBlocks, &numThreads, &sharedSize, maxGrid, s);

      kernel_MyIminLast<<<numBlocks, numThreads, sharedSize*2>>>(odata, opos, odata, opos, numThreads, IsPow2(s), s);
      s = (s + (numThreads*2-1)) / (numThreads*2);
   }
   cudaDeviceSynchronize(); 
   
   return opos[0];
}


void MySumPowsA(MyType *odata, MyType *idata, const int maxGrid, const int size)
{
   int numBlocks=0, numThreads=0, sharedSize=0, s;

   BlocksAndThreads(&numBlocks, &numThreads, &sharedSize, maxGrid, size);

   kernel_SumPows<<<numBlocks, numThreads, sharedSize>>>(odata, idata, numThreads, IsPow2(size), size);

   s = numBlocks;
   while (s > 1)
   {
      BlocksAndThreads(&numBlocks, &numThreads, &sharedSize, maxGrid, s);

      kernel_Sum<<<numBlocks, numThreads, sharedSize>>>(odata, odata, numThreads, IsPow2(s), s);
      s = (s + (numThreads*2-1)) / (numThreads*2);
 
   }
   kernel_Vnorm<<<1, 1>>>(odata);
}
