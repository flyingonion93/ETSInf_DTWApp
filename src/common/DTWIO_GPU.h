#pragma once
#ifndef DTW_IO_GPU_h
	#define DTW_IO_GPU_h

#include <stdio.h>
#include <stdlib.h>

#include "DTWCommon_GPU.h"

bool ReadParameters (DTWconst *, DTWfiles *, const char *);
bool ReadVector     (MyType   *, const int,  const char *);
bool ReadVectorInt64(int      *, const int,  const char *);
bool ReadS_fk       (MyType   *, const int, const char *);
bool ReadFrame      (MyType   *, const int,  FILE *);

bool CompareOutput(const MyType *, const int, const char *);

void FreeFiles(DTWfiles *);

bool escribirVectorEnFicheroTxtd(const MyType *, const int, const char *);
bool escribirVectorEnFicheroTxti(const int    *, const int, const char *);
#endif
