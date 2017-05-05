#pragma once

#ifndef DTW_Malloc_h
	#define DTW_Malloc_h

#include <stdio.h>
#include <stdlib.h>

#if FFTW
	#include <fftw3.h>
#endif

#include "DTWCommon.h"

int DTW_alloc( const DTW_const_st DTW_const, const int threads );

int DTW_dealloc( DTW_const_st DTW_const, DTW_files_st DTW_files, DTW_validation_files_st DTW_verify );

int DTW_init( precission_type **pD, precission_type **pV, int **pS, precission_type **cst, const int NST, const int NCST, const int TB );
#endif

