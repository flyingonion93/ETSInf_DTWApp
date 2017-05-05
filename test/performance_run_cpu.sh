#!/bin/bash

# Determine path
AppPATH=$(dirname `pwd`)

# Append data folders
BinPATH=$AppPATH/bin
DataPATH=$AppPATH/data
RawDataPATH=$DataPATH/Datos
ResDataPATH=$DataPATH/res

Time=$(date +"%Y-%m-%d-%H-%M-%S")

# Create a subfolder with the current timestamp
TestFolder=$ResDataPATH/perfomance_test_$Time
mkdir $TestFolder


# Define test parameters
BETA="0.0 1.0 1.5 1.7 2.0"
TIME="30sec 150sec 5min 10min 15min 30min"
N=20

# Execute the application.
# First using all the cores
# Then using only one core
until [ $N -eq 0 ];do
	echo "Iteration $N"
	echo "Iteration $N" >> $TestFolder/all_cores_$i-$j
	for i in $BETA
	do
		for j in $TIME
		do
			echo "(All cores) Beta $i with time $j"
			DTWApp_CPU $i $RawDataPATH/$j/parameters.dat >> $TestFolder/all_cores_$i-$j-execut
			line=$(tail -n1 $TestFolder/all_cores_$i-$j-execut)
			echo $line >> $TestFolder/all_cores_$i-$j
			rm $TestFolder/all_cores_$i-$j-execut
		done
	done
	let N-=1
done

export OMP_NUM_THREADS=1

N=20
until [ $N -eq 0 ];do
	echo "Iteration $N"
	echo "Iteration $N" >> $TestFolder/one_core_$i-$j
	for i in $BETA
	do
		for j in $TIME
		do
			echo "(One core) Beta $i with time $j"
			DTWApp_CPU $i $RawDataPATH/$j/parameters.dat >> $TestFolder/one_core_$i-$j-execut
			line=$(tail -n1 $TestFolder/one_core_$i-$j-execut)
			echo $line >> $TestFolder/one_core_$i-$j
			rm $TestFolder/one_core_$i-$j-execut
		done
	done
	let N-=1
done

echo "Performance test completed. Your files are stored in $TestFolder"
		
