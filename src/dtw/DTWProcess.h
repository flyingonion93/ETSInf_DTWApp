#pragma once

#ifndef DTW_Process_h
	#define DTW_Process_h

#include <omp.h>
#include <string.h>
#include <float.h>

#if FFTW
	#include <fftw3.h> //Tiene que irse fuera
#endif

#include "../common/DTWCommon.h"

int DTW_process( DTW_const_st DTW_const );

#endif
