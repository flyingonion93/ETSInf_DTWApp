#pragma once

#ifndef DTW_Preprocess_h
	#define DTW_Preprocess_H

#include <omp.h>
#include <math.h>
#include <string.h>

#if BLAS
	#include <cblas.h>
#endif

#if FFTW
	#include <fftw3.h>
#endif

#include "../common/DTWCommon.h"

int DTW_start_loop( DTW_const_st *DTW_const, const precission_type BETA );

int DTW_preprocess( DTW_const_st *DTW_const, const precission_type BETA );

#endif

