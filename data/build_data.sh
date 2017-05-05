#!/bin/bash

RAIZ=`pwd`

#Name of the file with parameters
NOMBRE=(Parameters_Lineal.dat Parameters_Circular.dat Parameters_GPU.dat)

#Directory for the double files and float files
TIPOS=(Datos Datosf)

#Normal data cases and their specific parameters
CASOS=(10min 150sec 15min 30min 5min 30sec)
BASES=(2114 713 2882 5241 1121 57)
ESTADOS=(5818 1528 9793 19496 2926 77)
TRAMAS=20000

#Integer Input files
INPUTINT=(k_max.bin k_min.bin states_seq.bin states_time_e.bin states_time_i.bin)

#Real Input files
INPUTREAL=(s_fk.bin v_hanning.bin)

#Verify files
VERIFY=(trama_hanning.bin v_cfreq.bin v_distortionxState.bin v_SxD.bin)


echo "Compiling doube2float.c and build ${TIPOS[1]}"
gcc $RAIZ/aux/doube2float.c -o doube2float

mkdir -p ${TIPOS[1]}
for i in ${CASOS[@]};
do
  mkdir -p ${TIPOS[1]}/$i
  mkdir -p ${TIPOS[1]}/$i/Input
  mkdir -p ${TIPOS[1]}/$i/Verif
  for j in ${INPUTINT[@]};
  do
    cp $RAIZ/${TIPOS[0]}/$i/Input/$j $RAIZ/${TIPOS[1]}/$i/Input/$j
  done
  for j in ${INPUTREAL[@]};
  do
    $RAIZ/doube2float $RAIZ/${TIPOS[0]}/$i/Input/$j $RAIZ/${TIPOS[1]}/$i/Input/$j > /dev/null
  done
  for j in ${VERIFY[@]};
  do
    $RAIZ/doube2float $RAIZ/${TIPOS[0]}/$i/Verif/$j $RAIZ/${TIPOS[1]}/$i/Verif/$j > /dev/null
  done
done

echo "Compiling bufferDtw.c and build trama files. WAIT...."
gcc $RAIZ/aux/bufferDtw.c -o bufferDtw

$RAIZ/bufferDtw 5700 570 $RAIZ/${TIPOS[0]}/20000muestras.bin $RAIZ/${TIPOS[0]}/20000tramas.bin > /dev/null

$RAIZ/doube2float $RAIZ/${TIPOS[0]}/20000tramas.bin $RAIZ/${TIPOS[1]}/20000tramas.bin > /dev/null
$RAIZ/doube2float $RAIZ/${TIPOS[0]}/9tramas.bin     $RAIZ/${TIPOS[1]}/9tramas.bin     > /dev/null

rm -f $RAIZ/bufferDtw
rm -f $RAIZ/doube2float

echo "Building parameters file for working"
i=0
for loop_i in ${TIPOS[@]};
do
  j=0
  for loop_j in ${CASOS[@]};
  do
    echo "114"        >  $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/${NOMBRE[0]}
    echo "114"        >  $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/${NOMBRE[1]}    
    echo ${BASES[$j]} >> $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/${NOMBRE[0]}
    echo ${BASES[$j]} >> $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/${NOMBRE[1]}
    echo ${BASES[$j]} >  $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/${NOMBRE[2]}

    k=0
    for loop_k in ${NOMBRE[@]};
    do
      echo ${ESTADOS[$j]} >> $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/${NOMBRE[k]}
      k=$(( $k + 1 ))
    done
    
    echo "5700" >> $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/${NOMBRE[0]}
    echo "5700" >> $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/${NOMBRE[1]}
    echo "570"  >> $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/${NOMBRE[0]}
    echo "570"  >> $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/${NOMBRE[1]}

    k=0
    for loop_k in ${NOMBRE[@]};
    do
      echo $TRAMAS >> $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/${NOMBRE[k]}
      k=$(( $k + 1 ))
    done

    echo "16384" >> $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/${NOMBRE[0]}
    echo "16384" >> $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/${NOMBRE[1]}

    echo "50"    >> $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/${NOMBRE[0]}

    echo "4"     >> $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/${NOMBRE[0]}
    echo "4"     >> $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/${NOMBRE[1]}

    k=0
    for loop_k in ${NOMBRE[@]};
    do
      echo "11.5129" >> $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/${NOMBRE[k]}
      k=$(( $k + 1 ))
    done

    k=0
    for loop_k in ${NOMBRE[@]};
    do
      echo $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/Input/v_hanning.bin     >> $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/${NOMBRE[k]}
      echo $RAIZ/${TIPOS[$i]}/20000tramas.bin                      >> $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/${NOMBRE[k]}
      echo $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/Input/s_fk.bin          >> $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/${NOMBRE[k]}
      echo $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/Input/k_max.bin         >> $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/${NOMBRE[k]}
      echo $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/Input/k_min.bin         >> $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/${NOMBRE[k]}
      echo $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/Input/states_time_e.bin >> $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/${NOMBRE[k]}
      echo $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/Input/states_time_i.bin >> $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/${NOMBRE[k]}
      echo $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/Input/states_seq.bin    >> $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/${NOMBRE[k]}
    
      echo $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/Verif/trama_hanning.bin      >> $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/${NOMBRE[k]}
      echo $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/Verif/v_cfreq.bin            >> $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/${NOMBRE[k]}
      echo $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/Verif/v_distortionxState.bin >> $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/${NOMBRE[k]}
      echo $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/Verif/v_SxD.bin              >> $RAIZ/${TIPOS[$i]}/${CASOS[$j]}/${NOMBRE[k]}
      k=$(( $k + 1 ))
    done

    j=$(( $j + 1 ))
  done
  i=$(( $i + 1 ))
done

