#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

int crearFicheroTramas(char *fileNameInput, char *fileNameOutput,int tamTrama, int tamMuestra){
 FILE *fpInput;
 FILE *fpOutput;
 
 
 long int tramaActual;
 long int totalTramas;
 int i;
 int datosLeidos;
 
 double *trama=(double *)malloc(tamTrama*sizeof(double));
 double *muestra=(double *)malloc (tamMuestra*sizeof(double));
 
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
 
 
 tramaActual=0;
 while(!feof(fpInput)){
   if(tramaActual==0){
     datosLeidos=fread(trama,sizeof(double),tamTrama,fpInput);
     if(datosLeidos==tamTrama){       
       tramaActual++;
     }    
   }
   else{
     datosLeidos=fread(muestra,sizeof(double),tamMuestra,fpInput);
     if(datosLeidos==tamMuestra){
       tramaActual++;
     }
   }   
 }
 
 totalTramas=tramaActual;
 fclose(fpInput);
 fpInput=fopen(fileNameInput,"rb");
 
// fwrite(&totalTramas,sizeof(long int),1,fpOutput);
// fwrite(&tamTrama,sizeof(long int),1,fpOutput);
 
 
 tramaActual=0;
 while(!feof(fpInput)){
   if(tramaActual==0){
     datosLeidos=fread(trama,sizeof(double),tamTrama,fpInput);
     if(datosLeidos==tamTrama){
       fwrite(trama,sizeof(double),tamTrama,fpOutput);       
       tramaActual++;
     }    
   }
   else{
     datosLeidos=fread(muestra,sizeof(double),tamMuestra,fpInput);
     if(datosLeidos==tamMuestra){
       for(i=0;i<tamTrama-tamMuestra;i++){
          trama[i]=trama[i+tamMuestra];
       }
       for(i=0;i<tamMuestra;i++){
         trama[(tamTrama-tamMuestra)+i]=muestra[i];
       }
       fwrite(trama,sizeof(double),tamTrama,fpOutput);           
       tramaActual++;
     }
   }   
   
 }
 
  
  fclose(fpInput);
  fclose(fpOutput);  
  free(trama);
  free(muestra);
  
  printf("Numero de tramas %ld \n",tramaActual);

  return tramaActual;
 
}

int crearFicheroTramasTxt(char *fileNameInput, char *fileNameOutput,char *fileNameOutputText,int tamTrama, int tamMuestra){
 FILE *fpInput;
 FILE *fpOutput;
 FILE *fpOutputText;
 
 long int tramaActual;
 long int totalTramas;
 int i;
 int datosLeidos;
 
 double *trama=(double *)malloc(tamTrama*sizeof(double));
 double *muestra=(double *)malloc (tamMuestra*sizeof(double));
 
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
 
 fpOutputText=fopen(fileNameOutputText,"w");
 if(fpOutputText==NULL){
   printf("Error en la apertura del fichero %s \n",fileNameOutputText);
   return -1;
 } 
 
 tramaActual=0;
 while(!feof(fpInput)){
   if(tramaActual==0){
     datosLeidos=fread(trama,sizeof(double),tamTrama,fpInput);
     if(datosLeidos==tamTrama){       
       tramaActual++;
     }    
   }
   else{
     datosLeidos=fread(muestra,sizeof(double),tamMuestra,fpInput);
     if(datosLeidos==tamMuestra){
       tramaActual++;
     }
   }   
 }
 
 totalTramas=tramaActual;
 fclose(fpInput);
 fpInput=fopen(fileNameInput,"rb");
 
 fwrite(&totalTramas,sizeof(long int),1,fpOutput);
 fprintf(fpOutputText,"%ld\n",totalTramas);
 fwrite(&tamTrama,sizeof(long int),1,fpOutput);
 fprintf(fpOutputText,"%d\n",tamTrama);
 
 tramaActual=0;
 while(!feof(fpInput)){
   if(tramaActual==0){
     datosLeidos=fread(trama,sizeof(double),tamTrama,fpInput);
     if(datosLeidos==tamTrama){
       fwrite(trama,sizeof(double),tamTrama,fpOutput);       
       for(i=0;i<tamTrama;i++) fprintf(fpOutputText,"%.8f\n",trama[i]);
       tramaActual++;
     }    
   }
   else{
     datosLeidos=fread(muestra,sizeof(double),tamMuestra,fpInput);
     if(datosLeidos==tamMuestra){
       for(i=0;i<tamTrama-tamMuestra;i++){
          trama[i]=trama[i+tamMuestra];
       }
       for(i=0;i<tamMuestra;i++){
         trama[(tamTrama-tamMuestra)+i]=muestra[i];
       }
       fwrite(trama,sizeof(double),tamTrama,fpOutput);
       for(i=0;i<tamTrama;i++) fprintf(fpOutputText,"%.8f\n",trama[i]);    
       tramaActual++;
     }
   }   
   
 }
 
  
  fclose(fpInput);
  fclose(fpOutput);
  fclose(fpOutputText);
  free(trama);
  free(muestra);
  
  printf("Numero de tramas %ld \n",tramaActual);

  return tramaActual;
 
}


void comprobarFichero(char *nameFileTramas){
  FILE *fpGijon;
  FILE *fpJaen;
  
  int tramaActual,i;
  
  int datosLeidos;
  
  int contErrores=0;
  
  int totalTramas;
  long int tamTrama;
  long int numTramas;
  
  double valorGijon,valorJaen;
  double *trama;
  
  char ficheroTrama[100];
  char numero[12];
  
  fpGijon=fopen(nameFileTramas,"rb");
  if(fpGijon==NULL){
     printf("Error en fichero %s \n",nameFileTramas);
     return;
  }
  
  datosLeidos=fread(&numTramas,sizeof(long int),1,fpGijon);
  datosLeidos=fread(&tamTrama,sizeof(long int),1,fpGijon);
  
  printf("Numero de tramas %ld, tamanio tramas %ld\n",numTramas,tamTrama);
  
  trama=(double *)malloc(tamTrama*sizeof(double));
  for(tramaActual=0;tramaActual<numTramas;tramaActual++){
    datosLeidos=fread(trama,sizeof(double),tamTrama,fpGijon);
        strcpy(ficheroTrama,"Jaen/trama");
        snprintf(numero, 12,"%d",tramaActual);
        strcat(ficheroTrama,numero);
        strcat(ficheroTrama,".bin");
        printf("%s\n",ficheroTrama);
    fpJaen=fopen(ficheroTrama,"rb");
    if(fpJaen==NULL){
       printf("Error en el fichero salida %s \n",ficheroTrama);
       return;
    }
    
    contErrores=0;
    for(i=0;i<10;i++){
      datosLeidos=fread(&valorJaen,sizeof(double),1,fpJaen);
      if(trama[i]!=valorJaen){
          if(contErrores==0) printf("\n\n\n\npos\tGijon\tJaen\n");      
          printf("[%d]\t%.8f\t%.8f\n",i,trama[i],valorJaen);
          contErrores++;
      }
    }
    fclose(fpJaen);
    printf("trama %d errores %d\n\n",tramaActual,contErrores);
    
  }
  
  
  fclose(fpGijon);
}

void pruebaTrama0(){
  FILE * fpJaen;
  FILE * fpGijon;
  
  int i,contError;
  
  double *tramaJaen=(double *)malloc(5700*sizeof(double));
  double *tramaGijon=(double *)malloc(5700*sizeof(double));
  
  fpJaen=fopen("Jaen/trama0.bin","rb");
  fpGijon=fopen("Datos/trama_0_a_9.bin","rb");
  
  if(fpJaen==NULL){
    printf("Error fichero jaen \n");
    return;
  }
  if(fpGijon==NULL){
    printf("error fichero gijon\n");
    return ;
  }
  
  fread(tramaJaen,sizeof(double),5700,fpJaen);
  fread(tramaGijon,sizeof(double),5700,fpGijon);
  printf("\n\nGijon  Jaen\n");
  contError=0;
  for(i=0;i<10;i++){
    if(tramaGijon[i]!=tramaJaen[i]){
       printf("[%d] %.8f %.8f \n",i,tramaGijon[i],tramaJaen[i]);
       contError++;
    }
    
  }
  
  printf("total Errores %d \n",contError);
  
  fclose(fpGijon);
  fclose(fpJaen);
  
  free(tramaJaen);
  free(tramaGijon);
  
}

int main(int argc, char *argv[]){
  // int tamMuestra=570, tamTrama  =5700;

  if (argc != 5) {
     printf("Uso: %s <tamTrama> <tamMuestra> <file muestras> <file tramas>. Ejemplo: 5700 570 muestra_input.bin trama_output.bin\n", argv[0]);
     return -1;
  }
  //crearFicheroTramas("FicherosMuestras/muestras15minutos.bin","FicherosTramas/tramas15minutos.bin", tamTrama, tamMuestra);
  crearFicheroTramas(argv[3], argv[4], atol(argv[1]), atol(argv[2]));
  return 0;
}

