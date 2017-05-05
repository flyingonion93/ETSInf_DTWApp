__device__ inline
double __shfl_downD(double var, unsigned int srcLane, int width=32)
{
  int2 a = *reinterpret_cast<int2*>(&var);

  a.x = __shfl_down(a.x, srcLane, width);
  a.y = __shfl_down(a.y, srcLane, width);

  return *reinterpret_cast<double*>(&a);
}

__inline__ __device__
double warpReduceSumD(double val)
{
  for (int offset = warpSize/2; offset > 0; offset /= 2)
    val += __shfl_downD(val, offset);
  return val;
}

__inline__ __device__
float warpReduceSumS(float val)
{
  for (int offset = warpSize/2; offset > 0; offset /= 2)
    val += __shfl_down(val, offset, 32);
  return val;
}


__global__  void kernel_InitDTW(MyType* __restrict__ pD, MyType* __restrict__ pV, int* __restrict__ pS,
                                const int pos, const int size)
{
   unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
     
   if (tid < size)
   {
      if (tid==pos)
      {
         pD[tid] = 0.0;
         pV[tid] = 0.0;
      }
      else
      {
        #ifdef SIMPLE
          pD[tid] = FLT_MAX;
          pV[tid] = FLT_MAX;      
        #else
          pD[tid] = DBL_MAX;
          pV[tid] = DBL_MAX;
        #endif
      }
      pS[tid] = 0.0;
   }
}


__global__ void kernel_DTW(const MyType* __restrict__ Sequence, MyType* __restrict__ pD, MyType* __restrict__ pV,
                           int* __restrict__ pS, const int NSeq, const int Where, const int NST) 
{
   int j=threadIdx.x + blockIdx.x * blockDim.x;
   int NSTplusNC, s, s2, k, Pos;

   MyType d, d2, v2;
   #ifdef SIMPLE
      MyType v=FLT_MAX;
   #else
      MyType v=DBL_MAX;
   #endif

   if( j<NST )
   {
      NSTplusNC = N_COSTS + NST;
      Pos       =((NSeq + N_COSTS) % TBLOCK) * NSTplusNC + N_COSTS + j - 1;
      for(k=0; k<N_COSTS; k++)
      {
         d2 = Sequence[j]*CCosts[k]+pD[Pos-k];
         s2 = 1  + pS[Pos-k];
         v2 = d2 / s2;

         if (v2 < v)  { d=d2; v=v2; s=s2; }
      }

      for (k=N_COSTS; k<T_COSTS; k++)
      {
         Pos=((NSeq + (T_COSTS-k)) % TBLOCK) * NSTplusNC + N_COSTS + j - 1;
         d2 = Sequence[j]*CCosts[k]+pD[Pos];
         s2 = 1  + pS[Pos];
         v2 = d2 / s2;

         if (v2 < v)  { d=d2; v=v2; s=s2; }
      }

      pD[Where+j] = d;
      pV[Where+j] = v;
      pS[Where+j] = s;
   }
}


__global__ void kernel_MyImin(MyType* __restrict__ odata, int* __restrict__ opos, const MyType* __restrict__ idata,
                              const int blockSize, const bool SizeIsPow2, const int size)
{
   extern __shared__ MyType ss[];
   MyType *sdata=ss;
   int    *pdata=(int *)&sdata[blockSize];

   int   tid = threadIdx.x;
   int     i = blockIdx.x*blockSize*2 + threadIdx.x;
   int gSize = blockSize*2*gridDim.x;
   int myPos;

   #ifdef SIMPLE
      MyType myMin=FLT_MAX;
   #else
      MyType myMin=DBL_MAX;
   #endif

   while (i < size)
   {
      myMin=idata[i];
      myPos=i;

      if (SizeIsPow2 || i + blockSize < size)
         if (idata[i+blockSize] < myMin) { myMin=idata[i+blockSize]; myPos=i+blockSize; }
         
      i += gSize;
   }
   sdata[tid]=myMin;
   pdata[tid]=myPos;
   __syncthreads();

   for (unsigned int s=256; s>=1; s/=2)
   {
      if ((blockSize >= s*2) && (tid < s))
         if (sdata[tid + s] < myMin) { sdata[tid]=myMin=sdata[tid+s]; pdata[tid]=myPos=pdata[tid+s]; }
      __syncthreads();
   }

   if (tid == 0) { odata[blockIdx.x]=myMin; opos[blockIdx.x]=myPos; }
}


__global__ void kernel_MyIminLast(MyType* __restrict__ odata, int* __restrict__ opos, const MyType* __restrict__ idata,
                                  const int* __restrict__ ipos, const int blockSize, const bool SizeIsPow2, const int size)
{
   extern __shared__ MyType ss[];
   MyType *sdata=ss;
   int    *pdata=(int *)&sdata[blockSize];

   int   tid = threadIdx.x;
   int     i = blockIdx.x*blockSize*2 + threadIdx.x;
   int gSize = blockSize*2*gridDim.x;
   int myPos;

   #ifdef SIMPLE
      MyType myMin=FLT_MAX;
   #else
      MyType myMin=DBL_MAX;
   #endif

   while (i < size)
   {
      myMin=idata[i];
      myPos=ipos[i];

      if (SizeIsPow2 || i + blockSize < size)
         if (idata[i+blockSize] < myMin) { myMin=idata[i+blockSize]; myPos=ipos[i+blockSize]; }
         
      i += gSize;
   }
   sdata[tid]=myMin;
   pdata[tid]=myPos;
   __syncthreads();

   for (unsigned int s=256; s>=1; s/=2)
   {
      if ((blockSize >= s*2) && (tid < s))
         if (sdata[tid + s] < myMin) { sdata[tid]=myMin=sdata[tid+s]; pdata[tid]=myPos=pdata[tid+s]; }
      __syncthreads();
   }

   if (tid == 0) { odata[blockIdx.x]=myMin; opos[blockIdx.x]=myPos; }
}


__global__ void kernel_SumPows(MyType* __restrict__ odata, const MyType* __restrict__ idata,
                               const int blockSize, const bool SizeIsPow2, const int size)
{
   extern __shared__ MyType sdata[];

   unsigned int      tid = threadIdx.x;
   unsigned int        i = blockIdx.x*blockSize*2 + threadIdx.x;
   unsigned int gridSize = blockSize*2*gridDim.x;

   MyType mySum=0.0;

   while (i < size)
   {
      mySum += idata[i]*idata[i];

      if (SizeIsPow2 || i + blockSize < size)
         mySum += idata[i+blockSize]*idata[i+blockSize];

      i += gridSize;
   }
   sdata[tid] = mySum;
   __syncthreads();

   if ((blockSize >= 512) && (tid < 256))
      sdata[tid] = mySum = mySum + sdata[tid + 256];
   __syncthreads();

   if ((blockSize >= 256) &&(tid < 128))
      sdata[tid] = mySum = mySum + sdata[tid + 128];
   __syncthreads();

   if ((blockSize >= 128) && (tid <  64))
      sdata[tid] = mySum = mySum + sdata[tid +  64];
   __syncthreads();

   if (tid < 32)
   {
      if (blockSize >=  64)
         mySum += sdata[tid + 32];

      for (int offset = sizeWarp/2; offset > 0; offset /= 2) 
         mySum += __shfl_down(mySum, offset);
   }

   if (tid == 0) odata[blockIdx.x] = mySum;
}


__global__ void kernel_Sum(MyType* __restrict__ odata, const MyType* __restrict__ idata,
                           const int blockSize, const bool SizeIsPow2, const int size)
{
   extern __shared__ MyType sdata[];

   unsigned int      tid = threadIdx.x;
   unsigned int        i = blockIdx.x*blockSize*2 + threadIdx.x;
   unsigned int gridSize = blockSize*2*gridDim.x;

   MyType mySum=0.0;

   while (i < size)
   {
      mySum += idata[i];

      if (SizeIsPow2 || i + blockSize < size)
         mySum += idata[i+blockSize];

      i += gridSize;
   }
   sdata[tid] = mySum;
   __syncthreads();

   if ((blockSize >= 512) && (tid < 256))
      sdata[tid] = mySum = mySum + sdata[tid + 256];
   __syncthreads();

   if ((blockSize >= 256) &&(tid < 128))
      sdata[tid] = mySum = mySum + sdata[tid + 128];
   __syncthreads();

   if ((blockSize >= 128) && (tid <  64))
      sdata[tid] = mySum = mySum + sdata[tid +  64];
   __syncthreads();

   if (tid < 32)
   {
      if (blockSize >=  64)
         mySum += sdata[tid + 32];

      for (int offset = sizeWarp/2; offset > 0; offset /= 2) 
         mySum += __shfl_down(mySum, offset);
   }

   if (tid == 0) odata[blockIdx.x] = mySum;
}


__global__ void kernel_Vnorm(MyType* __restrict__ odata)
{
   #ifdef SIMPLE
      odata[0] = 1.0f / (sqrtf(odata[0]) + FLT_EPSILON);
   #else
      odata[0] = 1.0  / ( sqrt(odata[0]) + DBL_EPSILON);
   #endif
}


__global__ void kernel_ApplyWin(MyType* __restrict__ dest, const MyType* __restrict__ src, const int size)
{
   unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;

   if (tid < size)
      dest[tid]=dest[tid] * src[tid];
}


__global__ void kernel_CalSxD1(MyType* __restrict__ dest, const MyType* __restrict__ src, const int* __restrict__ timei,
                               const int* __restrict__ timee, const int* __restrict__ seq, const int size)
{
   unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;

   MyType tmp;
   if (tid < size)
   {
      tmp = src[seq[tid]];
      for (unsigned int i=timei[tid]; i<=timee[tid]; i++)
         dest[i] = tmp;
   }
}


__global__ void kernel_CalSxD2(MyType* __restrict__ dest, const MyType ALPHA, const MyType* __restrict__ norm,
                               const int size)
{
   unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;

   if (tid < size)
      #ifdef SIMPLE
         dest[tid] = 1.0f - expf(ALPHA*fabsf(dest[tid]*norm[0]));
      #else
         dest[tid] = 1.0  -  exp(ALPHA* fabs(dest[tid]*norm[0]));
     #endif
}


__global__ void kernel_CalSxCuBlas1(MyType* __restrict__ dest, const MyType* __restrict__ src, const int* __restrict__ timei,
                                    const int* __restrict__ timee, const int* __restrict__ seq, const int size)
{
   unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;

   MyType tmp;
   if (tid < size)
   {
      tmp = src[seq[tid]];
      for (unsigned int i=timei[tid]; i<=timee[tid]; i++)
         dest[i] = tmp*tmp;
   }
}


__global__ void kernel_CalSxCuBlas2(MyType* __restrict__ dest, const MyType ALPHA, const MyType norm, const int size)
{
   unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;

   if (tid < size)
      #ifdef SIMPLE
         dest[tid] = 1.0f - expf(ALPHA*fabsf(sqrtf(dest[tid])*norm));
      #else
         dest[tid] = 1.0  -  exp(ALPHA* fabs( sqrt(dest[tid])*norm));
     #endif
}


__global__ void kernel_CopyAndSet(MyType* __restrict__ dest, const MyType* __restrict__ src,
                                  const int pos, const int size)
{
   unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;

   if (tid < size)
      dest[tid] = (tid < pos) ? src[tid] : 0.0;
}



__global__ void kernel_CompNorB0(MyType* __restrict__ norms, const MyType value, const int size)
{
   unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
   
   if (tid < size)
      norms[tid]=value;
}


__global__ void kernel_CompNorB1(MyType* __restrict__ norms, const MyType* __restrict__ s_fk,
                                 const int NMIDI, const int size)
{
   unsigned int    tid =  blockIdx.x * blockDim.x + threadIdx.x;
   unsigned int stride = (blockIdx.x * blockDim.x + threadIdx.x)*NMIDI;

   if (tid < size)
   {
      norms[tid]=0.0; 
      for(unsigned int j=0; j<NMIDI; j++)
         norms[tid] += s_fk[stride+j];
   }
}


__global__ void kernel_CompNorBG(MyType* __restrict__ norms, const MyType* __restrict__ s_fk,
                                 const int NMIDI, const MyType BETA, const int size)
{
   unsigned int    tid =  blockIdx.x * blockDim.x + threadIdx.x;
   unsigned int stride = (blockIdx.x * blockDim.x + threadIdx.x)*NMIDI;

   if (tid < size)
   {
      norms[tid]=0.0;
      for(unsigned int j=0; j<NMIDI; j++) 
         #ifdef SIMPLE
           norms[tid] += powf(s_fk[stride+j], BETA);
         #else
           norms[tid] +=  pow(s_fk[stride+j], BETA);
         #endif
   }   
}


__global__ void kernel_Modul(MyType* __restrict__ dest, const MyType* __restrict__ src, const int size)
{
  unsigned int    tid =  threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int stride = (threadIdx.x + blockIdx.x * blockDim.x)*2;
  
  MyType tmp1, tmp2;

  if (tid <= size)
  {
     tmp1 = src[stride];
     tmp2 = src[stride + 1];
     
     dest[tid]=tmp1*tmp1 + tmp2*tmp2;
  }
}

__global__ void kernel_Cfreq(MyType* __restrict__ dest, const MyType* __restrict__ src)
{
  unsigned int i = blockIdx.x;
  unsigned int j = threadIdx.x;

  MyType tmp = 0.0;
  for( unsigned int k=Ckmin_fft[i]+j; k<=Ckmax_fft[i]; k+=32 ) {
    tmp += src[k];
  }

  #ifdef SIMPLE
     tmp = warpReduceSumS(tmp);
  #else
     tmp = warpReduceSumD(tmp);
  #endif

  if( j==0 ) {
    #ifdef SIMPLE
       dest[i] = sqrtf(tmp);
    #else
       dest[i] =  sqrt(tmp);
    #endif
  }
}

__global__ void  kernel_Reduction(MyType* __restrict__ dest, const int size)
{
   unsigned int tid = threadIdx.x;

   MyType sum = dest[tid];

   for (unsigned int step=1; step<=size/blockDim.x; step++)
      if ((tid + step*blockDim.x) < size)
         sum = sum + dest[tid + step*blockDim.x];
   __syncthreads();

  #ifdef SIMPLE
     sum = warpReduceSumS(sum);
  #else
     sum = warpReduceSumD(sum);
  #endif
   
   if (tid == 0) dest[size] = sum;
}


__global__ void kernel_PowToReal(MyType* __restrict__ des, const MyType* __restrict__ src, const MyType ex, const int size)
{
   unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
   if (tid < size)
   {
      #ifdef SIMPLE
        des[tid]=powf(src[tid], ex);
      #else
        des[tid]= pow(src[tid], ex);
      #endif
   }
}

/* Paralelizacion bloques de 32x32 / 16*32  y sufle registros */
__global__ void kernel_CompDisB0new(MyType* __restrict__ dest, const MyType* __restrict__ v_cfreq,
                                 const MyType* __restrict__ norms, const MyType* __restrict__ s_fk,
                                 const int NMIDI, const int size)
{
   unsigned int      i =  blockIdx.x * blockDim.y + threadIdx.y;
   unsigned int      j;
   unsigned int stride = i * N_MIDI_PAD;
   unsigned int th_row = threadIdx.y;
   unsigned int th_col = threadIdx.x;
   unsigned int    row = i + threadIdx.x; /* This is useful only for the first row */
   bool          guard = th_row == 0 && row < size && th_col < blockDim.y;
   MyType a, b, tmp1;

   __shared__ MyType sh[32];

   if ( i < size )
   {
      a = 0.0;
      for(j=th_col;j<NMIDI;j+=32) {
        a += v_cfreq[j] * s_fk[stride+j];
      }
      #ifdef SIMPLE
         a = warpReduceSumS(a);
      #else
         a = warpReduceSumD(a);
      #endif
      if( guard ) {
        sh[th_col] = norms[row];
      }
      __syncthreads();
      if (th_col == 0) { b = a / sh[th_row]; }
      b = __shfl(b, 0);

      a = 0.0;
      for(j=th_col;j<NMIDI;j+=32) 
      {
         tmp1 = v_cfreq[j] / (s_fk[stride + j] * b);
         #ifdef SIMPLE
            a += tmp1 - logf(tmp1) - 1.0f;
         #else
            a += tmp1 -  log(tmp1) - 1.0;
         #endif
      }                                                          
      #ifdef SIMPLE
         a = warpReduceSumS(a);
      #else
         a = warpReduceSumD(a);
      #endif
      if( th_col == 0 ) {
        sh[th_row] = a;
      }
      __syncthreads();
      if( guard ) {
        dest[row] = sh[th_col];
      }
   }
}



__global__ void kernel_CompDisB0(MyType* __restrict__ dest, const MyType* __restrict__ v_cfreq,
                                 const MyType* __restrict__ norms, const MyType* __restrict__ s_fk,
                                 const int NMIDI, const int size)
{
   unsigned int      i =  blockIdx.x * blockDim.x + threadIdx.x, j;
   unsigned int stride = (blockIdx.x * blockDim.x + threadIdx.x)*N_MIDI_PAD;

   MyType A_kt=0.0, tmp1, tmp2=0.0;

   if (i < size)
   {
      for (j=0; j<NMIDI; j++)
         A_kt += v_cfreq[j] / s_fk[stride + j];
      A_kt = A_kt / norms[i];                              

      for (j=0; j<NMIDI; j++)
      {
         tmp1 = v_cfreq[j] / (s_fk[stride + j] * A_kt);
         #ifdef SIMPLE
           tmp2 += tmp1 - logf(tmp1) - 1.0f;
         #else
           tmp2 += tmp1 -  log(tmp1) - 1.0;
         #endif
      }
      dest[i]=tmp2;
   }
}


__global__ void kernel_CompDisB1(MyType* __restrict__ dest, const MyType* __restrict__ v_cfreq,
                                 const MyType* __restrict__ norms, const MyType* __restrict__ s_fk, 
                                 const int NMIDI, const int size)
{
   unsigned int      i =  blockIdx.x * blockDim.x + threadIdx.x, j;
   unsigned int stride = (blockIdx.x * blockDim.x + threadIdx.x)*N_MIDI_PAD;

   MyType tmp1, tmp2, tmp3, tmp4=0.0;

   if (i < size)
   {
      tmp1 = v_cfreq[NMIDI]/norms[i];

      for (j=0; j<NMIDI; j++)
      {
         tmp2 = s_fk[stride+j] * tmp1;
         tmp3 = v_cfreq[j];
         #ifdef SIMPLE
           tmp4 += tmp3*logf(tmp3/tmp2) + tmp2 - tmp3;
         #else
           tmp4 += tmp3* log(tmp3/tmp2) + tmp2 - tmp3;
         #endif
      }
      dest[i] = tmp4;
   }
}


__global__ void kernel_CompDisBGold(MyType* __restrict__ dest, const MyType* __restrict__ v_cfreq,
                                 const MyType* __restrict__ norms, const MyType* __restrict__ s_fk, 
                                 const MyType* __restrict__ ts_fk, const MyType* __restrict__ tauxi, 
                                 const MyType BETA, const int NMIDI, const int size)
{
   unsigned int      i =  blockIdx.x * blockDim.x + threadIdx.x, j;
   unsigned int stride = (blockIdx.x * blockDim.x + threadIdx.x)*N_MIDI_PAD;

   MyType tmp1, tmp2, tmp3;

   if (i < size)
   {
      tmp3 = 0.0;
      for (j=0; j<NMIDI; j++)
         tmp3 += v_cfreq[j] * ts_fk[stride+j];
      tmp3 = tmp3 / norms[i];

      #ifdef SIMPLE
         tmp1=powf(tmp3, BETA-1.0);
      #else
         tmp1= pow(tmp3, BETA-1.0);
      #endif

      tmp2=BETA*tmp1;
      tmp1=tmp1*tmp3*(BETA-1.0);
      tmp3=1.0/(BETA*(BETA-1.0));
      
      dest[i]= 0.0;
      for (j=0;j<NMIDI;j++)
        dest[i] += (tauxi[j] + ts_fk[stride+j]*(s_fk[stride+j]*tmp1 - v_cfreq[j]*tmp2))*tmp3;
   }
}

/* Paralelizacion bloques de 32x32 / 16*32  y sufle registros */
__global__ void kernel_CompDisBG(MyType* __restrict__ dest, const MyType* __restrict__ v_cfreq,
                                 const MyType* __restrict__ norms, const MyType* __restrict__ s_fk,
                                 const MyType* __restrict__ ts_fk, const MyType* __restrict__ tauxi,
                                 const MyType BETA, const int NMIDI, const int size)
{
   unsigned int      i =  blockIdx.x * blockDim.y + threadIdx.y;
   unsigned int      j;
   unsigned int stride = i * N_MIDI_PAD;
   unsigned int th_row = threadIdx.y;
   unsigned int th_col = threadIdx.x;
   unsigned int    row = i + threadIdx.x; /* This is useful only for the first row */
   bool          guard = th_row == 0 && row < size && th_col < blockDim.y;
   MyType a, b, tmp1, tmp2, tmp3;

   __shared__ MyType sh[32];

   if ( i < size )
   {
      a = 0.0;
      for(j=th_col;j<NMIDI;j+=32) {
        a += v_cfreq[j] * ts_fk[stride+j];
      }
      #ifdef SIMPLE
         a = warpReduceSumS(a);
      #else
         a = warpReduceSumD(a);
      #endif
      if( guard ) {
        sh[th_col] = norms[row];
      }
      __syncthreads();
      if( th_col == 0 ) {
        a = a / sh[th_row];
        #ifdef SIMPLE
           b = powf(a, BETA-1.0);
        #else
           b =  pow(a, BETA-1.0);
        #endif
        tmp1 = BETA * b;
        tmp2 = b * a * (BETA - 1.0);
      }
      tmp1 = __shfl(tmp1, 0);
      tmp2 = __shfl(tmp2, 0);
      tmp3 = 1.0 / (BETA * (BETA - 1.0));
      a = 0.0;
      for(j=th_col;j<NMIDI;j+=32) {
        a += ((tauxi[j] + ts_fk[stride+j] * (s_fk[stride+j] * tmp2 - v_cfreq[j] * tmp1)) * tmp3);
      }
      #ifdef SIMPLE
         a = warpReduceSumS(a);
      #else
         a = warpReduceSumD(a);
      #endif
      if( th_col == 0 ) {
        sh[th_row] = a;
      }
      __syncthreads();
      if( guard ) {
        dest[row] = sh[th_col];
      }
   }
}
