/*! \file DTWExec_CPU.c
 *	dtw
 *	
 *		\brief		Launcher for the DTW program for CPU architectures
 *		\author		Carlos Gómez Morillas
 *		\version	1.0
 *		\date		11/3/16
 *		\copyright	GNU Public License
 *		
 *		Last modification by Carlos Gómez Morillas on 22/8/16
 */

#include <stdio.h>
#include <stdarg.h>
#include <omp.h>

#include "DTWExec_CPU.h"
#include "../common/DTWIO.h"
#include "../common/DTWMalloc.h"
#include "../preproc/DTWPreprocess.h"
#include "../dtw/DTWProcess.h"

int DTWExec( int argc, char *argv[] )
{
	precission_type beta = (precission_type)1.5;
	int i;
	FILE *fp;	
	double total_time = 0.0;
#if TIMING
	double fft_time = 0.0;
	double distortion_time = 0.0;
	double vdi_time = 0.0;
	double dtw_time = 0.0;
	double acum_fft = 0.0;
	double acum_dis = 0.0;
	double acum_vdi = 0.0;
	double acum_dtw = 0.0;
#endif
	DTW_const_st  DTW_const;
	DTW_files_st  DTW_files;
	DTW_validation_files_st DTW_verify;

	if( 3 != argc )
	{
		printf( COLOR_YELLOW "Usage: %s <beta> <configuration file>. \nExample: main 1.5 parameters.dat\n" COLOR_RESET, argv[0] );
		return -EXIT_FAILURE;
	}
	printf( COLOR_YELLOW "DTW Testbench. Starting score alignment simulation\n" COLOR_RESET );

#if SIMPLE_PRECISSION
	beta = strtof( argv[1], NULL );
#else
	beta = strtod( argv[1], NULL );
#endif

	// DTWIO
	if( Config_DTW_data( &DTW_const, &DTW_files, &DTW_verify, argv[2] ) )
	{
		printf( COLOR_RED "Error starting the DTW structures\n" COLOR_RESET );
		return -EXIT_FAILURE;
	}
	printf( "Data configured\n" );

	// DTWMalloc
	if( DTW_alloc( DTW_const, omp_get_num_threads() ) )
	{
		printf( COLOR_RED "Error. There's an issue allocating memory for the DTW structures\n" COLOR_RESET );
		return -EXIT_FAILURE;
	}

	printf( "Memory allocated\n" );

	// DTWIO
	if( DTW_set_vectors( &DTW_const, &DTW_files ) )
	{
		printf( COLOR_RED "Error setting the values for the DTW vectors\n" COLOR_RESET );
		return -EXIT_FAILURE;
	}
	printf( "Vectors configured\n" );

	// DTWMalloc
	if( DTW_init( &pD, &pV, &pS, &costs, DTW_size - 1, DTW_const.nc, DTW_const.tb ) )
	{
		printf( COLOR_RED "Error. Some issues allocating memory for DTW structures\n" COLOR_RESET );
		return -EXIT_FAILURE;
	}
	printf( "DTW structures initialized\n" );


	frame = (precission_type *)calloc( DTW_const.frame_size, sizeof(precission_type ) );
	if( NULL == frame )
		return -EXIT_FAILURE;

	v_SxD = (precission_type *)calloc(DTW_size, sizeof(precission_type ) );
	if( NULL == v_SxD )
		return -EXIT_FAILURE;

	// DTWPreprocess
	if( DTW_start_loop( &DTW_const, beta ) )
	{
		printf( COLOR_RED "Error while initializing the process\n" COLOR_RESET );
		return -EXIT_FAILURE;
	}

	fp = fopen(DTW_files.file_frame, "rb" );
	if( NULL == fp )
	{
		printf( COLOR_RED "Error opening %s input file %s\n" , DTW_files.file_frame, COLOR_RESET );
		return -EXIT_FAILURE;
	}

	printf( "Applying alignment to the frames of the score. This could take a while...\n" );
	total_time=omp_get_wtime();
	for( i = 0; i < DTW_const.n_frames; i++ )
	{
		Read_frame( fp, frame, DTW_const.frame_size );

		// DTWPreprocess
		if( DTW_preprocess( &DTW_const, beta ) )
		{
			printf( COLOR_RED "Error while preprocessing the frame %d%s\n", i, COLOR_RESET );
			return -EXIT_FAILURE;
		}

		// DTWProcess
		if( DTW_process( DTW_const ) )
		{
			printf( COLOR_RED "Error while processing the frame %d%s\n", i, COLOR_RESET );
			return -EXIT_FAILURE;
		}
#if DEBUG
		printf("%d,%d\n", i, minimum_position);
#endif

	}
	total_time=omp_get_wtime() - total_time;

	// DTWIO
	if( Get_stats_output( &DTW_const, total_time ) )
		return -EXIT_FAILURE;

	fclose(fp);

	// DTWMalloc
	printf( "Deallocating\n" );
	if( DTW_dealloc( DTW_const, DTW_files, DTW_verify ) )
	{
		printf( COLOR_RED "Error. Some issues deallocating the DTW structures\n" COLOR_RESET );
		return -EXIT_FAILURE;
	}

	return 0;
}

