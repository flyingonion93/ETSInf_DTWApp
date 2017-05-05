#include "funcionesFicheros.h"

bool escribirVectorEnFicheroTxtd(const MyType *v, const int size, const char *filename)
{
   int i;

   FILE *fp=fopen(filename,"a");
   if (fp==NULL) return false;

   for(i=0; i<size; i++)
     #ifdef SIMPLE
        fprintf(fp, "%1.8E\n", v[i]);     
     #else
        fprintf(fp, "%1.15E\n",v[i]);
     #endif  
   fclose(fp);
   
   return true;
}

bool escribirVectorEnFicheroTxti(const int *v, const int size, const char *filename)
{
   int i;

   FILE *fp=fopen(filename,"a");
   if (fp==NULL) return false;

   for(i=0; i<size; i++)
     fprintf(fp, "%10d\n",v[i]);
   fclose(fp);
   
   return true;
}

bool ReadParameters(DTWconst *Param, DTWfiles *NameFiles, const char *filename)
{
   FILE *fp;

   int leidos=0;
   
   fp = fopen(filename,"r");
   if (fp==NULL) {
     printf("Error Opening file %s file\n", filename);   
     return false;
   }

   leidos += fscanf(fp, "%d\n", &Param->N_BASES);
   leidos += fscanf(fp, "%d\n", &Param->N_STATES);
   leidos += fscanf(fp, "%d\n", &Param->N_TRAMAS);
   #ifdef SIMPLE
     leidos += fscanf(fp, "%f\n",  &Param->ALPHA);
   #else
     leidos += fscanf(fp, "%lf\n", &Param->ALPHA);
   #endif
   Param->ALPHA = -(Param->ALPHA);

   NameFiles->file_hanning     =(char *)malloc(1024);
   NameFiles->file_trama       =(char *)malloc(1024);
   NameFiles->file_partitura   =(char *)malloc(1024);
   NameFiles->file_kmax        =(char *)malloc(1024);
   NameFiles->file_kmin        =(char *)malloc(1024);
   NameFiles->fileStates_Time_e=(char *)malloc(1024);
   NameFiles->fileStates_Time_i=(char *)malloc(1024);
   NameFiles->fileStates_seq   =(char *)malloc(1024);
   NameFiles->FileVerifyHanning=(char *)malloc(1024);
   NameFiles->FileVerifyFFT    =(char *)malloc(1024);
   NameFiles->FileVerifyVecDist=(char *)malloc(1024);
   NameFiles->FileVerifyDistor =(char *)malloc(1024);

   leidos += fscanf(fp, "%s\n", NameFiles->file_hanning);
   leidos += fscanf(fp, "%s\n", NameFiles->file_trama);
   leidos += fscanf(fp, "%s\n", NameFiles->file_partitura);
   leidos += fscanf(fp, "%s\n", NameFiles->file_kmax);
   leidos += fscanf(fp, "%s\n", NameFiles->file_kmin);
   leidos += fscanf(fp, "%s\n", NameFiles->fileStates_Time_e);
   leidos += fscanf(fp, "%s\n", NameFiles->fileStates_Time_i);
   leidos += fscanf(fp, "%s\n", NameFiles->fileStates_seq);
   leidos += fscanf(fp, "%s\n", NameFiles->FileVerifyHanning);
   leidos += fscanf(fp, "%s\n", NameFiles->FileVerifyFFT);
   leidos += fscanf(fp, "%s\n", NameFiles->FileVerifyVecDist);
   leidos += fscanf(fp, "%s\n", NameFiles->FileVerifyDistor);

   fclose(fp);
   if (leidos != 16) {
     printf("Error Rading from file %s.\n", filename);   
     return false;
   }
   else
     return true;
}

bool ReadVector(MyType *vector, const int size, const char *filename)
{
   FILE *fp;
   MyType valor;
   int contLineas, leidos;
   
   fp = fopen(filename,"rb");
   if (fp==NULL) return false;

   contLineas=0;

   leidos=fread(&valor, sizeof(MyType), 1, fp);
   if (leidos != 1) return false;
   
   while (!feof(fp)) {
      if (contLineas<size) vector[contLineas]=valor;
      contLineas++;
      
      leidos=fread(&valor, sizeof(MyType), 1, fp);     
      if ((leidos != 1) && (!feof(fp))) return false;
   }
   fclose(fp);

   if (contLineas == size)
      return true;
   else
      return false;
}

bool ReadVectorInt64(int *vector, const int size, const char *filename)
{
   int contLineas, leidos;
   FILE *fp;

   /* MODIFIED BY RANILLA 03-03-2016 13:54. Dirty, we need better solution */
   #ifdef ARM32
     long long int valorLong;
     int nbytes=sizeof(long long int);
   #else
     long int valorLong;
     int nbytes=sizeof(long int);
   #endif
  
   fp=fopen(filename, "rb");
   if (fp==NULL) return false;
   
   contLineas=0;
   
   leidos=fread(&valorLong, nbytes, 1, fp);
   if (leidos != 1) return false;
   
   while (!feof(fp))
   {
      if (contLineas < size)  vector[contLineas]=(int)valorLong;
      contLineas++;

      leidos=fread(&valorLong, nbytes, 1, fp);
      if ((leidos != 1) && (!feof(fp))) return false;
   }
   fclose(fp);    

   if (contLineas == size)
      return true;
   else
      return false;  
}

bool ReadS_fk(MyType *s_fk, const int BASES, const char *filename)
{
   long i, k;

   long size = N_MIDI_PAD*BASES;
   MyType data;

   FILE *fp;

   fp=fopen(filename, "rb");
   if (fp==NULL) return false;

   i=0;
   k=fread(&data, sizeof(MyType), 1, fp);
   if (k != 1) return false;

   while(!feof(fp))
   {
      if (i<size)
      {
        s_fk[i]=data;
        if ((i%N_MIDI_PAD)< (N_MIDI-1))
           i++;
        else
           i += (N_MIDI_PAD-N_MIDI+1);
      }

      k=fread(&data, sizeof(MyType), 1, fp);
      if ((k != 1) && (!feof(fp))) return false;
   }
   fclose(fp);

   if (i == size)
      return true;
   else
      return false;
}

bool ReadFrame(MyType *trama, const int size, FILE *fp)
{
   int leidos;

   leidos = fread(trama, sizeof(MyType), size, fp);
   
   if (leidos == size)
      return true;
   else
      return false;
}

void FreeFiles(DTWfiles *NameFiles)
{
   free(NameFiles->file_hanning);
   free(NameFiles->file_trama);
   free(NameFiles->file_partitura);
   free(NameFiles->file_kmax);
   free(NameFiles->file_kmin);
   free(NameFiles->fileStates_Time_e);
   free(NameFiles->fileStates_Time_i);
   free(NameFiles->fileStates_seq);
   free(NameFiles->FileVerifyHanning);
   free(NameFiles->FileVerifyFFT);
   free(NameFiles->FileVerifyVecDist);
   free(NameFiles->FileVerifyDistor);
}


bool CompareOutput(const MyType *vector, const int size, const char *filename)
{
   FILE *fp;
   int contValores, leidos;
   MyType valor, Error=0.0;

   fp=fopen(filename,"rb");
   if(fp==NULL) return false;

   Error      =0.0;
   contValores=0;

   leidos=fread(&valor, sizeof(MyType), 1, fp);
   if (leidos != 1) return false;

   while(!feof(fp) && (contValores < size))
   {
      Error+=(valor - vector[contValores])*(valor - vector[contValores]);
      contValores++;
      
      leidos=fread(&valor, sizeof(MyType), 1, fp);
      if ((leidos != 1) && (!feof(fp))) return false;
   }
   fclose(fp);

   #ifdef SIMPLE
      Error = sqrtf(Error)/(MyType)contValores;
   #else
      Error =  sqrt(Error)/(MyType)contValores;
   #endif
   printf("El error es: %1.14E\n", Error);
   return true;
}
