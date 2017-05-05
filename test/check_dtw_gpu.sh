#!/bin/sh

# Determine path
AppPATH=$(dirname `pwd`)

# Append data folders
BinPATH=$AppPATH/bin
DataPATH=$AppPATH/data
RawDataPATH=$DataPATH/Datos
RawfDataPATH=$DataPATH/Datosf
ResDataPATH=$DataPATH/res

# Define test parameters
TestParam=$RawDataPATH/30sec/ParamVerify_Circular.dat
TestfParam=$RawfDataPATH/30sec/ParamVerify_Circular.dat
Beta=1.5
Threads=64
Memory=1
Cublas=1


Time=$(date +"%Y-%m-%d-%H:%M:%S")
FileName=execution_$Time
FileNamef=executionf_$Time

#Implementación original
$BinPATH/DTWAppGPU $Beta $TestParam $Threads $Memory $Cublas >> $ResDataPATH/$FileName
$BinPATH/DTWAppGPUf $Beta $TestfParam >> $ResDataPATH/$FileNamef

#
diff $ResDataPATH/$FileName $ResDataPATH/benchmark/execution_double >> $ResDataPATH/diff_res
diff $ResDataPATH/$FileNamef $ResDataPATH/benchmark/execution_float >> $ResDataPATH/diff_resf

lined=$(head -n 1 $ResDataPATH/diff_res)
linef=$(head -n 1 $ResDataPATH/diff_resf)

if [ "$lined" = "1001c1001" ]; then
	echo "(Double) Correct application execution"
else
	echo "(Double) Wrong application execution. There is some error with the processing code"
fi

if [ "$linef" = "1001c1001" ]; then
	echo "(Float) Correct application execution"
else
	echo "(Float) Wrong application execution. There is some error with the processing code"
fi

# Cleanup
rm $ResDataPATH/$FileName
rm $ResDataPATH/$FileNamef
rm $ResDataPATH/diff_res
rm $ResDataPATH/diff_resf

