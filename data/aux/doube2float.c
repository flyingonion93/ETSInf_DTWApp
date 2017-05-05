#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

int crearFichero(char *fileNameInput, char *fileNameOutput)
{
 FILE *fpInput;
 FILE *fpOutput;
 
 double DatoI;
 float  DatoO;
 int    Leidos;
 long int cantidad; 
 
 fpInput=fopen(fileNameInput,"rb");
 if(fpInput==NULL){
   printf("Error en la lectura del fichero %s \n",fileNameInput);
   return -1;
 } 
 
 fpOutput=fopen(fileNameOutput,"wb");
 if(fpOutput==NULL){
   printf("Error en la apertura del fichero %s \n",fileNameOutput);
   return -1;
 }
 
 cantidad = 0; 
 Leidos=fread(&DatoI, sizeof(double), 1, fpInput);
 while(!feof(fpInput))
 {
    DatoO=(float)DatoI;
    Leidos=fwrite(&DatoO, sizeof(float), 1, fpOutput);
    Leidos=fread (&DatoI, sizeof(double), 1, fpInput);
    cantidad ++;
 }
 
  fclose(fpInput);
  fclose(fpOutput);  
  
  printf("Cantidad de elementos leidos/escritos %ld \n", cantidad);
  
  return 0;
}

int main(int argc, char *argv[]){
  if (argc != 3) {
     printf("Uso: %s <Input file> <Output file>. Ejemplo: muestra_input.bin muestras_output.bin\n", argv[0]);
     return -1;
  }
  crearFichero(argv[1], argv[2]);
  return 0;
}
