#!/bin/bash

RAIZ=`pwd`

#Name of the file with parameters
NOMBRE=(ParamVerify_Lineal.dat ParamVerify_Circular.dat ParamVerify_GPU.dat)

#Directory for the double files and float files
TIPOS=(Datos Datosf)

#Normal data cases and their specific parameters
CASOS=(10min 150sec 15min 30min 5min 30sec)
BASES=(2114 713 2882 5241 1121 57)
ESTADOS=(5818 1528 9793 19496 2926 77)

# Put your trama number for verify
TRAMAS=1000

echo "Building parametersVerify file for testing"
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

