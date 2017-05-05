/*! \file DTWMalloc.c
 *	dtw
 *
 *		\brief		Memory allocation operations for the DTW structures
 *		\author		Carlos Gómez Morillas
 *		\version	1.0
 *		\date		11/8/16
 *		\copyright	GNU Public License
 *
 *		Last modification by Carlos Gómez Morillas on 30/4/16
 */

#include "DTWMalloc.h"

/*
 * Private module interface
 */
static int Init_data( const DTW_const_st DTW_const );

static int End_data( DTW_const_st *DTW_const, DTW_files_st *DTW_files, DTW_validation_files_st *DTW_verify );

static int DTW_end();

static int FFT_init( const int NFFT, const int N_MIDI, const int threads );

static int FFT_end();

static int FFT_alloc( const int NFFT, const int N_MIDI );

static int FFT_dealloc();
#if FFTW
static int FFTW_init( const int NFFT, const int threads );

static int FFTW_end();
#endif

/*!
 *	\brief User function that allocates memory for the vectors required for preprocessing and DTW
 *
 *	\param [in] DTW_const		Structure that contains the parameters of the appllication
 *	\param [in] threads		Number of threads that the application is using
 */
int DTW_alloc( const DTW_const_st DTW_const, const int threads )
{
	if( Init_data( DTW_const ) )
		return -EXIT_FAILURE;		

	if( FFT_init( DTW_const.nfft, DTW_const.n_midi, threads ) )
		return -EXIT_FAILURE;

	return EXIT_SUCCESS;
}

/*!
 *	\brief User function that frees the memory for the vectors required for preprocessing and DTW
 *
 *	\param [in] DTW_const		Structure that contains the parameters of the application
 *	\param [in] DTW_files		Structure that contains the names of the files used by the application
 */
int DTW_dealloc( DTW_const_st DTW_const, DTW_files_st DTW_files, DTW_validation_files_st DTW_verify )
{
	if( End_data( &DTW_const, &DTW_files, &DTW_verify ) )
		return -EXIT_FAILURE;

	if( DTW_end() )
		return -EXIT_FAILURE;

	if( FFT_end() )
		return -EXIT_FAILURE;

	return EXIT_SUCCESS;
}

/*!
 *	\brief Allocates memory for the DTW-related structures and starts them with a fixed value to 
 *	adjust to the alignment process
 *	
 *	\param [in,out] pD		Description goes here
 *	\param [in,out] pV		Description goes here
 *	\param [in,out] pS		Description goes here
 *	\param [in,out] costs		Description goes here
 *	\param [in] NST			Number of states of a frame
 *	\param [in] NCST		Number of DTW dependencies
 *	\param [in] TB			Extended buffer size
 */
int DTW_init( precission_type **pD, precission_type **pV, int **pS, precission_type **costs, const int NST, const int NCST, const int TB )
{
	int i;
	int size = ( NST + NCST ) * ( NCST + TB );
	(*pD) = (precission_type *)calloc( sizeof(precission_type), size );
	if( NULL == (*pD) )
		return -EXIT_FAILURE;

   	(*pV) = (precission_type *)calloc( sizeof(precission_type), size );
	if( NULL == (*pV) )
		return -EXIT_FAILURE;

	(*pS) = (int *)calloc( size, sizeof(int) );
	if( NULL == (*pS) )
		return -EXIT_FAILURE;

	(*costs) = (precission_type *)calloc( 2 * NCST - 1, sizeof(precission_type) );
	if( NULL == (*costs) )
		return -EXIT_FAILURE;

	for( i = 0; i < size; i++ )
	{
#if SIMPLE_PRECISSION
		(*pD)[i] = FLT_MAX;
		(*pV)[i] = FLT_MAX;
#else
		(*pD)[i] = DBL_MAX;
		(*pV)[i] = DBL_MAX;
#endif
	}

	(*pD)[( NST + NCST ) * NCST - NST ] = (precission_type)0.00f;
	(*pV)[( NST + NCST ) * NCST - NST ] = (precission_type)0.00f;

	(*costs)[0] = 1.00f;
	(*costs)[1] = 1.05f;
	(*costs)[2] = 1.10f;
	(*costs)[3] = 1.20f;
	(*costs)[4] = 1.05f;
	(*costs)[5] = 1.10f;
	(*costs)[6] = 1.20f;

	return EXIT_SUCCESS;
}


/*!
 *	\brief Allocates memory for the common structures and vectors used
 *
 *	\param [in] DTW_const		Common values structure
 */
static int Init_data( const DTW_const_st DTW_const )
{
	v_hanning = (precission_type *)calloc( DTW_const.frame_size, sizeof( precission_type ) );
	if( NULL == v_hanning )
		return -EXIT_FAILURE;
	
	states_time_i = (int *)calloc( DTW_const.n_states, sizeof(int) );
	if( NULL == states_time_i )
		return -EXIT_FAILURE;

	states_time_e = (int *)calloc( DTW_const.n_states, sizeof(int) );
	if( NULL == states_time_e )
		return -EXIT_FAILURE;
	
	states_seq = (int *)calloc( DTW_const.n_states, sizeof(int) );
	if( NULL == states_seq )
		return -EXIT_FAILURE;

	s_fk = (precission_type *)calloc( DTW_const.n_midi * DTW_const.n_bases, sizeof(precission_type) );
	if( NULL == s_fk )
		return -EXIT_FAILURE;

	norms = (precission_type *)calloc( DTW_const.n_bases, sizeof(precission_type ) );
	if( NULL == norms )
		return -EXIT_FAILURE;

	v_distortion_x_states=(precission_type *)calloc( DTW_const.n_bases, sizeof(precission_type ) );
	if( NULL == v_distortion_x_states )
		return -EXIT_FAILURE;

	v_cfreq = (precission_type *)calloc( DTW_const.n_midi, sizeof(precission_type ) );  
	if( NULL == v_cfreq )
		return -EXIT_FAILURE;

      	aux_frame = (precission_type *)calloc( DTW_const.n_midi, sizeof(precission_type ) );
	if( NULL == aux_frame )
		return -EXIT_FAILURE;

	ts_fk = (precission_type *)calloc( DTW_const.n_midi * DTW_const.n_bases, sizeof(precission_type) );
	if( NULL == ts_fk )
		return -EXIT_FAILURE;

	return EXIT_SUCCESS;
}

/*!
 *	\brief Deallocates the memory used by the DTW structures
 *
 */
static int DTW_end()
{
	free( pD );
	free( pV );
	free( pS );
	free( costs );

	return EXIT_SUCCESS;
}

/*!
 *	\brief Deallocates the memory used by common structures and vectors
 *	
 *	\param [in] *DTW_const		Common values structure
 *	\param [in] *DTW_files		Name files structure
 *	\param [in] *DTW_verofy		Description goes here
 */
static int End_data( DTW_const_st *DTW_const, DTW_files_st *DTW_files, DTW_validation_files_st *DTW_verify )
{
	free( v_hanning );
	free( states_time_i );
	free( states_time_e );
	free( states_seq );
	free( s_fk );
	free( norms );
	free( v_distortion_x_states );
	free( v_SxD );
	free( v_cfreq );
	free( frame );
	free( aux_frame );
	free( ts_fk );

	free( DTW_files->file_hanning );
	free( DTW_files->file_frame );
	free( DTW_files->file_score );
	free( DTW_files->file_kmax );
	free( DTW_files->file_kmin );
	free( DTW_files->file_states_time_e );
	free( DTW_files->file_states_time_i );
	free( DTW_files->file_states_seq );

#if VERIFY
	free( DTW_verify->FileVerifyHanning );
	free( DTW_verify->FileVerifyFFT );
	free( DTW_verify->FileVerifyVecDist );
	free( DTW_verify->FileVerifyDistor );
#endif
	
	return EXIT_SUCCESS;
}

/*!
 *	\brief Starts the required structures and vectors for the FFT computation
 *
 *	\param [in] NFFT		Size of the FFT matrix
 *	\param [in] N_MIDI		Number of MIDI samples
 *	\param [in] threads		Number of threads that the application is using 
 */
static int FFT_init( const int NFFT, const int N_MIDI, const int threads )
{
	if( FFT_alloc( NFFT, N_MIDI ) )
		return -EXIT_FAILURE;

#if FFTW
	if( FFTW_init( NFFT, threads ) )
		return -EXIT_FAILURE;
#endif

	return EXIT_SUCCESS;
}

/*!
 *	\brief Frees all of the structures used by the FFT computation
 *
 */
int FFT_end()
{
	if( FFT_dealloc() )
		return -EXIT_FAILURE;

#if FFTW
	if( FFTW_end() )
		return -EXIT_FAILURE;
#endif

	return EXIT_SUCCESS;
}

/*!
 *	\brief Allocates memory for the FFT vectors
 *	
 *	\param [in] NNFT		Size of the FFT matrix
 *	\paran [in] N_MIDI		Number of MIDI samples 
 */
int FFT_alloc( const int NFFT, const int N_MIDI )
{
	X_fft = (precission_type *)calloc( 2 * NFFT + 1, sizeof( precission_type ) );
	if( NULL == X_fft )
		return -EXIT_FAILURE;

	mod_fft = (precission_type *)calloc( NFFT, sizeof( precission_type ) );
	if( NULL == mod_fft )
		return -EXIT_FAILURE;

   	kmax_fft = (int *)calloc(N_MIDI, sizeof(int));
	if( NULL == kmax_fft )
		return -EXIT_FAILURE;

	kmin_fft=(int *)calloc(N_MIDI, sizeof(int));
	if( NULL == kmin_fft )
		return -EXIT_FAILURE;

	return EXIT_SUCCESS;
}

/*!
 *	\brief Frees the used memory of the common FFT vectors
 *
 */
int FFT_dealloc()
{
	free( X_fft );
	free( mod_fft );
	free( kmin_fft );
	free( kmax_fft );

	return EXIT_SUCCESS;
}

#if FFTW
/*!
 *	\brief Library-dependent function. Allocates memory for the structures required by FFTW
 *
 *	\param [in] NFFT		Size of the FFT matrix
 *	\param [in] threads		Number of threads that the application is using 
 */
static int FFTW_init( const int NFFT, const int threads )
{
#if SIMPLE_PRECISSION
	fftwf_init_threads();
	out_fft = (precission_type *)fftwf_malloc( sizeof(precission_type) *NFFT );
	if( NULL == out_fft )
		return -EXIT_FAILURE;

	fftwf_plan_with_nthreads( threads );
	plan = fftwf_plan_r2r_1d( NFFT, X_fft, out_fft, FFTW_R2HC, FFTW_MEASURE );
#else
	fftw_init_threads();
	out_fft  = (precission_type *)fftw_malloc( sizeof(precission_type) * NFFT );
	if( NULL == out_fft )
		return -EXIT_FAILURE;

	fftw_plan_with_nthreads( threads );
	plan = fftw_plan_r2r_1d( NFFT, X_fft, out_fft, FFTW_R2HC, FFTW_MEASURE );
#endif
	if( NULL == plan )
	{
		printf( "Error creating the plan\n" );
		return -EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/*!
 *	\brief Library-dependent function. Frees memory for the structures used by FFTW
 *	
 */
static int FFTW_end()
{
#if SIMPLE_PRECISSION
	fftwf_free( out_fft );
	fftwf_destroy_plan( plan );
	fftwf_cleanup_threads();
#else
	fftw_free(out_fft);
	fftw_destroy_plan(plan);
	fftw_cleanup_threads();	
#endif
	return EXIT_SUCCESS;
}
#endif

