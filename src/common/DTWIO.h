#pragma once

#ifndef DTW_IO_h
	#define DTW_IO_h

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>

#if FFTW
	#include <fftw3.h> //Tiene que irse fuera
#endif

#include "DTWCommon.h"

int Get_stats_output( DTW_const_st *DTW_const, double total_time );

int Config_DTW_data( DTW_const_st *DTW_const, DTW_files_st *DTW_files, DTW_validation_files_st *DTW_verify, const char *file_name );

int DTW_set_vectors( DTW_const_st *DTW_const, DTW_files_st *DTW_files );

void Read_frame( FILE *fp, precission_type *current_frame, const int FRAME_SIZE );

#endif

