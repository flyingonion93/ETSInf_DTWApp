/*! \file DTWCommon.h
 *	dtw
 *	
 *		\brief 		Common structures and variables for the DTW program
 *		\author		Carlos Gómez Morillas
 * 		\version	1.0
 *		\date		11/3/16
 *		\copyright	GNU Public License
 * 
 *		Last modification by Carlos Gómez Morillas on 19/8/16
 */
#pragma once

#ifndef DTW_Common_h
	#define DTW_Common_h
#if SIMPLE_PRECISSION
	#define precission_type float
#if FFTW
	#define FFTW_plan_type fftwf_plan
#endif
#else
	#define precission_type	double
#if FFTW
	#define FFTW_plan_type fftw_plan
#endif
#endif

#define COLOR_RED     "\x1b[31m"
#define COLOR_GREEN   "\x1b[32m"
#define COLOR_YELLOW  "\x1b[33m"
#define COLOR_BLUE    "\x1b[34m"
#define COLOR_RESET   "\x1b[0m"

#include <stdlib.h>
#include <float.h>
#include <stdio.h>

/*!
 *	\brief Struct to store common values for DTW from a configuration file.
 *	Each composition needs a file with values for
 *	these parameters. The structure of this text file should be
 *	the same as this DTW_const_st, DTW_files_st, DTW_validation_files_st
 *	structs, one line per parameter. No empty lines and one empty line to the end
 *
 *	\par Example of this file:
 * 	114
 *	57
 *	77
 *	5700
 *	570
 *	1
 *	16384
 *	50
 *	4
 *	11.5129
 *	Datos/30sec/Input/hanning.bin
 *	....
 *	Datos/30sec/Verif/v_SxD.bin
 *	
 *	\param n_midi					Number of MIDI samples. The default value is 114
 *	\param n_bases					Number of bases in the score
 *	\param n_states					Number of states in the score
 *	\param n_frames			 		Number of frames in the score
 *	\param frame_size				Size of every frame. The default value is 5700
 *	\param sample_size				Size of every sample. The default value is 5700
 *	\param nfft					Size of FFT matrix. It must be a power of two for performance purposes. The default value is 16384
 *	\param tb					Extended buffers size
 *	\param nc					Number of costs per dimension. The default value is 4
 *	\param alpha					Program constant used to compute the distortions
 */
struct DTW_const_st
{
	int n_midi;
	int n_bases;
	int n_states;
	int n_frames;
	int frame_size;
	int sample_size;
	int nfft;
	int tb;
	int nc;
	precission_type alpha;
};

/*!
 *	\brief Struct to store the file names for each required step of the process
 *
 *	\param file_hanning				Name of the hanning vector file
 *	\param file_frame				Name of the frame file
 *	\param file_score				Name of the score file
 *	\param file_kmax				Parameter description
 *	\param file_kmin				Parameter description
 *	\param file_states_time_e			Parameter description
 *	\param file_states_time_i			Parameter description
 *	\param file_states_seq				Parameter description
 */
struct DTW_files_st
{
	char *file_hanning;
	char *file_frame;
	char *file_score;
	char *file_kmax;
	char *file_kmin;
	char *file_states_time_e;
	char *file_states_time_i;
	char *file_states_seq;
};

/*!
 *	\brief Struct to store the file names for validation
 *	
 *	\param file_validation_hanning	Parameter description
 *	\param file_validation_fft	Parameter description
 *	\param file_validation_vec_dist	Parameter description
 *	\param file_validation_distor	Parameter description	
 */
struct DTW_validation_files_st
{
	char *file_validation_hanning;
	char *file_validation_fft;
	char *file_validation_vec_dist;
	char *file_validation_distor;
};

//DTW
int *pS;				/*!< Vector that stores the number of jump of the score*/
int DTW_size;				/*!< */
int minimum_position;			/*!< Minimum position obtained by the DTW algorithm*/
int TBDes;				/*!< Actual position on buffers*/
precission_type *pD;			/*!< Vector that stores the distances between states*/
precission_type *pV;			/*!< Vector that stores the quotient of distance and jumps*/
precission_type *costs;			/*!< Vector that stores the cost of each state*/
//DTW

//FFT
precission_type	*X_fft;			/*!< */
precission_type *out_fft;		/*!< */
precission_type	*mod_fft;		/*!< */
int *kmax_fft;				/*!< */
int *kmin_fft;				/*!< */
//FFT

#if FFTW
FFTW_plan_type plan;			/*!< FFT solver for the library FFTW3*/
#endif

//COMMON
precission_type	*norms;			/*!< */
precission_type *s_fk;			/*!< */
precission_type *v_SxD;			/*!< */
precission_type *v_cfreq;		/*!< */
precission_type	*v_hanning;		/*!< */
precission_type	*v_distortion_x_states;	/*!< */
//COMMON
int *states_time_i;			/*!< */
int *states_time_e;			/*!< */
int *states_seq;			/*!< */
//COMMON

precission_type *frame;			/*!< Frame that the DTW processes*/
precission_type *aux_frame;		/*!< Auxiliar frame to use when moving the vector on memory*/
precission_type *ts_fk;			/*!< */

typedef struct DTW_const_st DTW_const_st;
typedef struct DTW_files_st DTW_files_st;
typedef struct DTW_validation_files_st DTW_validation_files_st;
#endif

