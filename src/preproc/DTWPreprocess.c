/*! \file DTWPreprocess.c
 *	dtw
 *
 *		\brief		Preprocessing functions of the application
 *		\author		Carlos Gómez Morillas
 *		\version	1.0
 *		\date		20/8/16
 *		\copyright	GNU Public License
 *
 *		Last modification by Carlos Gómez Morillas on 30/5/16
 */
#include "DTWPreprocess.h"

/*
 * Private module interface
 */
static int Init_auxi( const int N_MIDI, const int N_BASES, const precission_type BETA );

static void Compute_norms( const int N_MIDI, const int N_BASES, const precission_type BETA );

static void Apply_window( precission_type * __restrict frame, const int FRAME_SIZE );

static void FFT( const precission_type *trama, const int TAMTRAMA, const int NFFT, const int N_MIDI );

static void Compute_distortion_vector( precission_type * __restrict v_distortion_x_states, const precission_type BETA, const int N_BASES, const int N_MIDI );

static void Apply_distortion( precission_type *v_SxD, const int N_STATES, const precission_type ALPHA );

#if FFTW
static void FFTW_();
#endif
/*!
 *	\brief User function that starts the required vectors before the main application loop
 *	
 *	\param [in] DTW_const		Common values structure
 */
int DTW_start_loop( DTW_const_st *DTW_const, const precission_type BETA )
{
	if( Init_auxi( DTW_const->n_midi, DTW_const->n_bases, BETA ) )
		return -EXIT_FAILURE;

	Compute_norms( DTW_const->n_midi, DTW_const->n_bases, BETA );

	FFT( frame, DTW_const->frame_size, DTW_const->nfft, DTW_const->n_midi );

	return EXIT_SUCCESS;
}

/*!
 *	\brief User function that computes all of the required operations of the preprocess step
 *
 *	\param [in] *DTW_const		Common values structure
 *	\param [in] BETA			Parameter description
 */
int DTW_preprocess( DTW_const_st *DTW_const, const precission_type BETA )
{
	Apply_window( frame, DTW_const->frame_size );

	FFT( frame, DTW_const->frame_size, DTW_const->nfft, DTW_const->n_midi );

	Compute_distortion_vector( v_distortion_x_states, BETA, DTW_const->n_bases, DTW_const->n_midi );

	Apply_distortion( v_SxD, DTW_const->n_states, DTW_const->alpha );

	return EXIT_SUCCESS;
}

/*!
 *	\brief Starts the vector ts_fk to later use on the preprocessing
 *	
 *	\param [in] N_MIDI			Number of MIDI samples
 *	\param [in] N_BASES			Number of bases in the score
 *	\param [in] BETA			Parameter description
 */
static int Init_auxi( const int N_MIDI, const int N_BASES, const precission_type BETA )
{
	int i;
	if( (precission_type)1.5 >= BETA && (precission_type)1.5 <= BETA )
	{
		#pragma omp parallel for
		for( i = 0; i < N_MIDI * N_BASES; i++ )
#if SIMPLE_PRECISSION
			ts_fk[i] = sqrtf( s_fk[i] );
#else
			ts_fk[i] = sqrt( s_fk[i] );
#endif
	}else if( (precission_type)0.0 != BETA && (precission_type)1.0 != BETA )
	{
		#pragma omp parallel for
		for( i = 0; i < N_MIDI * N_BASES; i++ )
#if SIMPLE_PRECISSION
			ts_fk[i] = powf( s_fk[i], BETA - 1.0f );
#else
			ts_fk[i] = pow( s_fk[i], BETA - 1.0 );
#endif
	}
	return 0;
}

static void Compute_norms( const int N_MIDI, const int N_BASES, const precission_type BETA )
{
	int i;
	
	if( BETA >= (precission_type)0.0 && BETA <= (precission_type)0.0 )
	{
		for( i = 0; i < N_BASES; i++ )
			norms[i] = (precission_type) N_MIDI;
	}else if( BETA >= (precission_type)1.0 && BETA <= (precission_type)1.0 )
	{
		#pragma omp parallel for
		for( i = 0; i < N_BASES; i++ )
		{
			int k = i * N_MIDI;
#if BLAS
#if SIMPLE_PRECISSION
			norms[i] = cblas_sasum( N_MIDI, &s_fk[k], 1 );
#else
			norms[i] = cblas_dasum( N_MIDI, &s_fk[k], 1 );
#endif
#endif
		}
	}else if( BETA >= (precission_type)2.0 && BETA <= (precission_type)2.0 )
	{
		#pragma omp parallel for
		for( i = 0; i < N_BASES; i++ )
		{
			int k = i * N_MIDI;
#if BLAS
#if SIMPLE_PRECISSION
			norms[i] = cblas_sasum( N_MIDI, &s_fk[k], 1 );
#else
			norms[i] = cblas_dasum( N_MIDI, &s_fk[k], 1 );
#endif
#endif
		}
	}else if( BETA >= (precission_type)1.5 && (precission_type)BETA <= 1.5 )
	{
		#pragma omp parallel for
		for( i = 0; i < N_BASES; i++ )
		{
			int j;
			int k;
			precission_type data;
			k = i * N_MIDI;

			data = (precission_type)0.0;
			for( j = 0; j < N_MIDI; j++ )
#if SIMPLE_PRECISSION
			data += s_fk[k+j] * sqrtf( s_fk[k+j] );
#else
			data += s_fk[k+j] * sqrt( s_fk[k+j] );
#endif
			norms[i] = data;
		}
	}else
	{
		#pragma omp parallel for
		for( i = 0; i < N_BASES; i++ )
		{
			int j;
			int k;
			precission_type data;
			k = i * N_MIDI;
	
			data = (precission_type)0.0;
			for( j = 0; j < N_MIDI; j++ )
#if SIMPLE_PRECISSION
			data += powf( s_fk[k+j], BETA );
#else
			data += pow( s_fk[k+j], BETA );
#endif
			norms[i] = data;
		}
	}
}

static void Apply_window( precission_type * __restrict frame, const int FRAME_SIZE )
{
	int i;

	#pragma omp parallel for if(FRAME_SIZE > 1000)
	for( i = 0; i < FRAME_SIZE; i++ )
		frame[i] = frame[i] * v_hanning[i];

}

/*!
 *	\brief Computes the FFT for the current frame. It depends on some third-party libraries
 *
 *	\param [in] FRAME_SIZE			Size of the current frame
 *	\param [in] NFFT			Size of FFT matrix
 *	\param [in] N_MIDI			Number of MIDI samples
 */
static void FFT( const precission_type *frame, const int FRAME_SIZE, const int NFFT, const int N_MIDI )
{
	int i;

	memcpy( X_fft, frame, sizeof(precission_type) * FRAME_SIZE );
	memset( &X_fft[FRAME_SIZE], 0, sizeof(precission_type) * ( NFFT - FRAME_SIZE ) );

#if FFTW
	FFTW_();
#endif	
	mod_fft[0] = out_fft[0] * out_fft[0];
	mod_fft[NFFT / 2] = out_fft[NFFT / 2] *out_fft[NFFT / 2];

	#pragma omp parallel for
	for( i = 1; i < NFFT / 2; i++ )
		mod_fft[i] = out_fft[i] * out_fft[i] + out_fft[NFFT - i] * out_fft[NFFT - i];

	#pragma omp parallel for
	for( i =0; i < N_MIDI; i++ )
	{
#if BLAS
#if SIMPLE_PRECISSION
		v_cfreq[i] = sqrtf( cblas_sasum( kmax_fft[i] - kmin_fft[i] + 1, &mod_fft[kmin_fft[i]], 1) );
#else
		v_cfreq[i] = sqrt( cblas_dasum( kmax_fft[i] - kmin_fft[i] + 1, &mod_fft[ kmin_fft[i]], 1 ) );
#endif
#endif
	}
}

/*!
 *	\brief Computes the distortion vector for the problem
 *
 *	\param [out] *v_distortion_x_states		Parameter description
 *	\param [in] BETA						Parameter description
 *	\param [in] N_BASES						Number of bases in the score
 *	\param [in] N_MIDI						Number of MIDI samples
 */
static void Compute_distortion_vector( precission_type * __restrict v_distortion_x_states, const precission_type BETA, const int N_BASES, const int N_MIDI )
{
	int i;

	if( (precission_type)0.0 >= BETA && (precission_type)0.0 <= BETA )
	{
		#pragma omp parallel for
		for( i = 0; i < N_BASES; i++ )
		{
			int j;
			int itmp = N_MIDI * i;
			precission_type A_kt = 0.0;
			precission_type dsum = 0.0;
			precission_type dtmp;

			//BETA = 0 --> a^(BETA - 1) = a^-1 = 1 / a
			for( j = 0; j < N_MIDI; j++ )
				A_kt += v_cfreq[j] / s_fk[itmp+j];
			
			A_kt = A_kt / norms[i];
			for( j = 0; j < N_MIDI; j++ )
			{
				dtmp = v_cfreq[j] / ( s_fk[itmp + j] * A_kt );
#if SIMPLE_PRECISSION
				dsum += dtmp - logf( dtmp ) - 1.0f;
#else
				dsum += dtmp -  log( dtmp ) - 1.0;
#endif
			}
			v_distortion_x_states[i] = dsum;
		}
	}else if( (precission_type)1.0 >= BETA && (precission_type)1.0 <= BETA )
	{
		precission_type A_kt = (precission_type)0.0;
			
		/*	BETA = 1 --> a^(BETA - 1) = a^0 = 1 due to a >= 0. So next inner-loop
			is moved here (out) because it doesn't depend on index i, is always the
			same result/operation/value
		*/
#if BLAS
#if SIMPLE_PRECISSION
		A_kt = cblas_sasum( N_MIDI, v_cfreq, 1 );
#else
		A_kt = cblas_dasum( N_MIDI, v_cfreq, 1 );
#endif
#endif
      		#pragma omp parallel for
		for( i = 0; i < N_BASES; i++ )
		{
			int j;
			int itmp  = N_MIDI * i;
			precission_type dsum = 0.0;
			precission_type dtmp;
			precission_type dtmp2 = A_kt / norms[i];
      
			for( j = 0; j < N_MIDI; j++ )
			{
				dtmp = s_fk[itmp + j] * dtmp2;
#if SIMPLE_PRECISSION
			dsum += v_cfreq[j] * logf( v_cfreq[j] / dtmp ) + dtmp - v_cfreq[j];
#else
			dsum += v_cfreq[j] * log( v_cfreq[j] / dtmp ) + dtmp - v_cfreq[j];
#endif
			}
			v_distortion_x_states[i] = dsum;
		}
	}else if( (precission_type)2.0 <= BETA && (precission_type)2.0 >= BETA )
	{
		for( i = 0; i < N_MIDI; i++ )
			aux_frame[i] = v_cfreq[i] * v_cfreq[i];

		#pragma omp parallel for
		for( i = 0; i < N_BASES; i++ )
		{
			int j;
			int itmp = N_MIDI*i;
			precission_type A_kt = 0.0;
			precission_type dsum = 0.0;	
			precission_type dtmp;
			//BETA = 2 --> a^(BETA - 1) = a^1 = a			
#if BLAS
#if SIMPLE_PRECISSION
			A_kt = cblas_sdot( N_MIDI, v_cfreq, 1,  &s_fk[itmp], 1 );
#else
			A_kt = cblas_ddot( N_MIDI, v_cfreq, 1,  &s_fk[itmp], 1 );
#endif
#endif
			A_kt = A_kt / norms[i];   
			for( j = 0; j < N_MIDI; j++ )
			{
				dtmp = s_fk[itmp + j] * A_kt;
				dsum += (aux_frame[j] + dtmp * dtmp - ((precission_type)2.0 * v_cfreq[j] * dtmp)) * (precission_type)0.5;
			}
			v_distortion_x_states[i] = dsum;
		}
	}else if( (precission_type)1.5 >= BETA && (precission_type)1.5 <= BETA )
	{
		precission_type d_const = (precission_type)1.0 / ( (precission_type)1.5 *(precission_type)0.5 );
		for( i = 0; i < N_MIDI; i++ )
#if SIMPLE_PRECISSION
			aux_frame[i] = v_cfreq[i] * sqrtf( v_cfreq[i] );
#else
			aux_frame[i] = v_cfreq[i] * sqrt( v_cfreq[i] );
#endif
		#pragma omp parallel for
		for( i = 0; i < N_BASES; i++ )
		{
			int j;
			int itmp = N_MIDI*i;
			precission_type A_kt = 0.0;
			precission_type dsum = 0.0;
			precission_type dtmp;
			precission_type dtmp2;
			precission_type dtmp3;
#if BLAS
#if SIMPLE_PRECISSION
			A_kt = cblas_sdot( N_MIDI, v_cfreq, 1,  &s_fk[itmp], 1 ) / norms[i];
			dtmp3 = sqrtf( A_kt );
#else
			A_kt = cblas_ddot( N_MIDI, v_cfreq, 1,  &s_fk[itmp], 1 ) / norms[i];
			dtmp3 = sqrt(A_kt);
#endif
#ifdef SIMPLE_PRECISSION
			dtmp3 = sqrtf( A_kt );
#else
			dtmp3 = sqrt( A_kt );
#endif
#endif
			for( j = 0; j < N_MIDI; j++ )
			{
				dtmp  = s_fk[itmp + j] * A_kt;
				dtmp2 = ts_fk[itmp + j] * dtmp3;
				dsum += ( aux_frame[j] + (precission_type)0.5 * dtmp * dtmp2 - (precission_type)1.5 * v_cfreq[j] * dtmp2) * d_const;
			}
			v_distortion_x_states[i] = dsum;
		}
	}else
	{
		precission_type beta_minus_one = BETA - (precission_type)1.0;
		precission_type	d_const = (precission_type)1.0 / ( BETA * ( beta_minus_one ) );
		for( i = 0; i < N_MIDI; i++ )
#if SIMPLE_PRECISSION
			aux_frame[i] = powf( v_cfreq[i], BETA );
#else
			aux_frame[i] = pow( v_cfreq[i], BETA );
#endif
		#pragma omp parallel for
		for( i = 0; i < N_BASES; i++ )
		{
			int j;
			int itmp = N_MIDI*i;
			precission_type A_kt = (precission_type)0.0;
			precission_type dsum = (precission_type)0.0;
			precission_type dtmp;
			precission_type dtmp2;
			precission_type dtmp3;

#if BLAS
#if SIMPLE_PRECISSION
			A_kt = cblas_sdot( N_MIDI, v_cfreq, 1,  &s_fk[itmp], 1 ) / norms[i];
			dtmp3 = powf( A_kt, beta_minus_one );
#else
			A_kt = cblas_ddot( N_MIDI, v_cfreq, 1,  &s_fk[itmp], 1 ) / norms[i];
			dtmp3 = pow( A_kt, beta_minus_one );
#endif
#if SIMPLE_PRECISSION
			dtmp3 = powf( A_kt, beta_minus_one );
#else
			dtmp3 = pow( A_kt, beta_minus_one );
#endif
#endif

			for( j = 0; j < N_MIDI; j++ )
			{
				dtmp  = s_fk[itmp + j] * A_kt;
				dtmp2 = ts_fk[itmp + j] * dtmp3;
				dsum += ( aux_frame[j] + beta_minus_one * dtmp2 * dtmp - BETA * v_cfreq[j] * dtmp2 ) * d_const;
			}
			v_distortion_x_states[i] = dsum;
		}
	}
}

/*!
 *	\brief Applies a distortion 
 *
 *	\param [in] N_STATES				Parameter description
 *	\param [in] ALPHA				Parameter description
 */
static void Apply_distortion( precission_type *v_SxD, const int N_STATES, const precission_type ALPHA )
{
	int i;
	int position;
	int size = states_time_e[N_STATES - 1] + 1;

	precission_type vnorm = (precission_type)0.0;

	#pragma omp parallel for private(position) reduction(+: vnorm)
	for( i = 0; i < N_STATES; i++ )
	{
		for( position = states_time_i[i]; position <= states_time_e[i]; position++ )
		{
			v_SxD[position] = v_distortion_x_states[states_seq[i]];
			vnorm += v_SxD[position] * v_SxD[position];
		}
	}
	//ALPHA is defined as (-1.0) * ALPHA in function Read_structure
#if SIMPLE_PRECISSION
	vnorm = 1.0f / ( sqrtf( vnorm ) + FLT_EPSILON );
	
	#pragma omp parallel for
	for( i = 0; i < size; i++ )
		v_SxD[i] = 1.0f - expf( ALPHA * fabsf( v_SxD[i] * vnorm ) );
#else
	
	vnorm = 1.0 /  ( sqrt( vnorm ) + DBL_EPSILON );

	#pragma omp parallel for
	for( i = 0; i < size; i++ )
		v_SxD[i] = 1.0 - exp( ALPHA * fabs( v_SxD[i] * vnorm ) );
#endif

}

#if FFTW
/*!
 *	\brief Library-dependent function that obtains the FFT for the score based on an specified FFTW plan
 *
 */
static void FFTW_()
{
#if SIMPLE_PRECISSION
	fftwf_execute( plan );
#else
	fftw_execute( plan );
#endif
}
#endif

