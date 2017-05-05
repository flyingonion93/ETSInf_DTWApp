/*! \file DTWIO.c
 *	dtw
 *
 *		\brief		I/O operations to write on structs or text files and read the required DTW vectors
 *		\author		Carlos Gómez Morillas
 *		\version	1.0
 *		\date		12/3/16
 *		\copyright	GNU Public License
 *
 *		Last modification by Carlos Gómez Morillas on 30/5/16
 */
#include "DTWIO.h"

/*
 * Private module interface
 */
static int Allocate_file_data_names( DTW_files_st *DTW_files, DTW_validation_files_st *DTW_verify );

static int Read_structure( DTW_const_st *DTW_const, FILE *fp );

static int Read_file_data( DTW_files_st *DTW_files, FILE *fp );

static int Read_vector( const char *file_name, precission_type *vector, const int size );

static int Read_vector_int_64( const char *file_name, int *vector, const int size );

static int Read_s_fk( const int N_MIDI, const int N_BASES, const char *file_name );

static int Read_int_64_vectors( DTW_const_st *DTW_const, DTW_files_st *DTW_files );

#if VERIFY
static int Read_verify_data( DTW_validation_files_st *DTW_verify, FILE *fp );
#endif

/*!
 *	\brief User function that reads the configuration file and sets the structures containing
 *	the application constants and vector file names
 *
 *	\param [out] *DTW_const		Pointer to common values structure
 *	\param [out] *DTW_files		Pointer to file names structures
 *	\param [out] *DTW_verify	Pointer to validation file names structure
 *	\param [in]  *file_name		Name of the configuration file
 */
int Config_DTW_data( DTW_const_st *DTW_const, DTW_files_st *DTW_files, DTW_validation_files_st *DTW_verify, const char *file_name )
{
	FILE *fp;
	fp = fopen( file_name, "r" );

	if( NULL == fp )
		return -EXIT_FAILURE;

	if( Allocate_file_data_names( DTW_files, DTW_verify ) )
		return EXIT_FAILURE;

	if( Read_structure( DTW_const, fp ) )
		return -EXIT_FAILURE;

	if( Read_file_data( DTW_files, fp ) )
		return -EXIT_FAILURE;

#if VERIFY
	if( Read_verify_data( DTW_verify, fp ) )
		return -EXIT_FAILURE;
#endif

	fclose( fp );

	return EXIT_SUCCESS;
}

/*!
 *	\brief User function to set the value of the required vectors based on their files
 *	
 *	\param [in] *DTW_const		Pointer to common values structure
 *	\param [in] *DTW_files		Pointer to file names structure
 */
int DTW_set_vectors( DTW_const_st *DTW_const, DTW_files_st *DTW_files )
{
	if( Read_vector( DTW_files->file_hanning, v_hanning, DTW_const->frame_size ) )
		return -EXIT_FAILURE;

	if( Read_int_64_vectors( DTW_const, DTW_files ) )
		return -EXIT_FAILURE;

	if( Read_s_fk( DTW_const->n_midi, DTW_const->n_bases, DTW_files->file_score ) )
		return -EXIT_FAILURE;
	
	DTW_size = states_time_e[DTW_const->n_states - 1] + 1;

	return EXIT_SUCCESS;
}

/*!
 *	\brief Reads the current frame of the score
 *	
 *	\param [in] *fp				File descriptor
 *	\param [out] *current_frame		Pointer to the vector that stores frames
 *	\param [in] FRAME_SIZE			Size of the frame that we are reading
 */
void Read_frame( FILE *fp, precission_type *current_frame, const int FRAME_SIZE )
{
	fread( current_frame, sizeof(precission_type), FRAME_SIZE, fp );
}

int Get_stats_output( DTW_const_st *DTW_const, double total_time )
{
	char* colored_time;
	double time_per_frame = total_time * 1000.0 / DTW_const->n_frames;
	if( time_per_frame > 1.28E+01 )
		colored_time = COLOR_RED;
	
	else if( time_per_frame <= 1.28E+01 && time_per_frame >= 1E+01 )
		colored_time = COLOR_YELLOW;

	else
		colored_time = COLOR_GREEN;

	printf( COLOR_BLUE "------------ STATS ------------\n" COLOR_RESET );
	printf( COLOR_BLUE "# of frames: %s\t %d\n", COLOR_RESET, DTW_const->n_frames );
	printf( COLOR_BLUE "Total time: %s\t %1.7E sec.\n", COLOR_RESET, total_time );
	printf( COLOR_BLUE "Time per frame:\t %s %1.7E ms.  %s\n", colored_time, time_per_frame, COLOR_RESET );
	printf( COLOR_BLUE "-------------------------------\n" COLOR_RESET );

	return EXIT_SUCCESS;
}

/*!
 *	\brief Reads the configuration file and store all of their values on a common structure
 *
 *	\param [out] *DTW_const		Pointer to common values structure
 *	\param [in] *fp			File descriptor
 */
static int Read_structure( DTW_const_st *DTW_const, FILE *fp )
{
	fscanf( fp, "%d\n",  &DTW_const->n_midi );
	fscanf( fp, "%d\n",  &DTW_const->n_bases );
	fscanf( fp, "%d\n",  &DTW_const->n_states );
	fscanf( fp, "%d\n",  &DTW_const->frame_size );	
	fscanf( fp, "%d\n",  &DTW_const->sample_size );
	fscanf( fp, "%d\n",  &DTW_const->n_frames );
	fscanf( fp, "%d\n",  &DTW_const->nfft );
	fscanf( fp, "%d\n",  &DTW_const->tb );
	fscanf( fp, "%d\n",  &DTW_const->nc );
#if SIMPLE_PRECISSION
	fscanf( fp, "%f\n",  &DTW_const->alpha );
#else
	fscanf( fp, "%lf\n", &DTW_const->alpha );
#endif
	DTW_const->alpha = - ( DTW_const->alpha );

	return EXIT_SUCCESS;
}

/*!
 *	\brief Sets the file name data on the structure
 *
 *	\param [out] *DTW_files		Pointer to the structure that stores the file names
 *	\param [in] *fp			File descriptor
 */
static int Read_file_data( DTW_files_st *DTW_files, FILE *fp )
{
	fscanf(fp, "%s\n", DTW_files->file_hanning);
	fscanf(fp, "%s\n", DTW_files->file_frame);
	fscanf(fp, "%s\n", DTW_files->file_score);
	fscanf(fp, "%s\n", DTW_files->file_kmax);
	fscanf(fp, "%s\n", DTW_files->file_kmin);
	fscanf(fp, "%s\n", DTW_files->file_states_time_e);
	fscanf(fp, "%s\n", DTW_files->file_states_time_i);
	fscanf(fp, "%s\n", DTW_files->file_states_seq);

	return EXIT_SUCCESS;
}

#if VERIFY
/*!
 *	\brief Description goes here
 *
 *	\param [out] *DTW_verify	Pointer to the structure that stores the names of the test files 
 *	\param [in] *fp				File descriptor
 */
static int Read_verify_data( DTW_validation_files_st *DTW_verify, FILE *fp )
{
	fscanf( fp, "%s\n", DTW_verify->file_validation_hanning );
	fscanf( fp, "%s\n", DTW_verify->file_validation_fft );
	fscanf( fp, "%s\n", DTW_verify->file_validation_vec_dist );
	fscanf( fp, "%s\n", DTW_verify->file_validation_distor );

	return EXIT_SUCCESS;
}
#endif

/*!
 *	\brief Reads a specified file containing the vector values and sets them on a pointer
 *
 *	\param [in] *file_name	Name of the file which contains the values of the vector
 *	\param [out] *vector	Pointer to the vector that we are writing
 *	\param [in] size		Size of the vector
 */
static int Read_vector( const char *file_name, precission_type *vector, const int size )
{
	FILE *fp;
	precission_type value;
	int line_count;
	int read;

	fp = fopen( file_name, "rb" );
	if( NULL == fp )
	{
		printf( "ERROR. File \"%s\" not found \n", file_name );
		return -EXIT_FAILURE;
	}

	line_count = 0;
	read = fread( &value, sizeof(precission_type), 1, fp );
	if( 1 != read )
		return -EXIT_FAILURE;

	while( !feof( fp ) )
	{
		if( line_count < size )
			vector[line_count] = value;

		line_count++;
		read = fread( &value, sizeof( precission_type ), 1 , fp );
		if( ( 1 != read ) && ( !feof( fp ) ) )
			return -EXIT_FAILURE;
	}
	fclose( fp );
      
	if( line_count  == size )
		return EXIT_SUCCESS;
	
	else
		return line_count;
}

/*!
 *	\brief Reads a specified file containing the values of a long integer and sets them on a pointer
 *
 *	\param [in] *file_name		Name of the configuration file
 *	\param [out] *vector		Pointer to the vector that we are reading
 *	\param [in] size		Size of the vector
 */
static int Read_vector_int_64( const char *file_name, int *vector, const int size )
{
	FILE *fp;
	int line_count;
	int read;

#if ARM32
	long long int long_value;
	int n_bytes = sizeof(long long int);
#else
	long int long_value;
	int n_bytes = sizeof(long int);
#endif
	fp = fopen( file_name, "rb" );
	if( NULL == fp )
	{
		printf( "ERROR. File \"%s\" not found \n", file_name );
		return -EXIT_FAILURE;
	}

	line_count = 0;
	read = fread( &long_value, n_bytes, 1, fp );
	if( 1 != read )
		return -EXIT_FAILURE;

	while( !feof( fp ) )
	{
		if( line_count < size )
			vector[line_count] = (int)long_value;
	
		line_count++;
		read = fread( &long_value, n_bytes, 1, fp );
		if( ( 1 != read ) && ( !feof( fp ) ) )
			return -EXIT_FAILURE;
	}
	fclose( fp );

	if( line_count == size )
		return EXIT_SUCCESS;
	else
		return line_count;
}

/*!
 *	\brief Reads the <vector> from its file
 *
 *	\param [out] *s_fk			Pointer to <vector>
 *	\param [in] N_MIDI			Number of MIDI samples
 *	\param [in] N_BASES			Number of bases on the score
 *	\param [in] *file_name		Name of the file which contains the values of the vector
 */
static int Read_s_fk( const int N_MIDI, const int N_BASES, const char *file_name )
{
	FILE *fp;
	int i = 0;
	int k;
	precission_type data;
	
	fp = fopen( file_name, "rb" );
	if( NULL == fp )
		return -EXIT_FAILURE;
	
	k = fread( &data, sizeof(precission_type), 1, fp );
	if( 1 != k )
		return -EXIT_FAILURE;

	while( !feof( fp ) ) 
	{
		if( i < ( N_BASES * N_MIDI ) )
			s_fk[i] = data;
		
		i++;
		k = fread( &data, sizeof(precission_type), 1, fp );
		if( ( 1 != k ) && ( !feof( fp ) ) )
			return -EXIT_FAILURE;
	}
	fclose(fp);
	if( N_BASES * N_MIDI == i )
		return EXIT_SUCCESS;

	else
		return i; 
}

/*!
 *	\brief Writes every long vector depending on their file name
 *
 *	\param [in] *DTW_const		Pointer to common values structure
 *	\param [in] *DTW_files		Pointer to file names structure
 */
static int Read_int_64_vectors( DTW_const_st *DTW_const, DTW_files_st *DTW_files )
{
	if( Read_vector_int_64( DTW_files->file_states_seq, states_seq ,DTW_const->n_states ) )
		return EXIT_FAILURE;

	if( Read_vector_int_64( DTW_files->file_states_time_i, states_time_i, DTW_const->n_states ) )
		return EXIT_FAILURE;

	if( Read_vector_int_64( DTW_files->file_states_time_e, states_time_e, DTW_const->n_states ) )
		return EXIT_FAILURE;

	if( Read_vector_int_64( DTW_files->file_kmax, kmax_fft, DTW_const->n_midi ) )
		return EXIT_FAILURE;

	if( Read_vector_int_64( DTW_files->file_kmin, kmin_fft, DTW_const->n_midi ) )
		return EXIT_FAILURE;

	return EXIT_SUCCESS;
}

/*!
 *	\brief Allocates memory for the parameters of the file names structures.
 *
 *	\param [in] *DTW_files		Pointer to the structure that stores the file names
 *	\param [in] *DTW_verify		Pointer to the structure that stores the names of the test files 
 */
static int Allocate_file_data_names( DTW_files_st *DTW_files, DTW_validation_files_st *DTW_verify )
{
	DTW_files->file_hanning = (char *) malloc( 1024 );
	DTW_files->file_frame = (char *) malloc( 1024 );
	DTW_files->file_score = (char *) malloc( 1024 );
	DTW_files->file_kmax = (char *) malloc( 1024 );
	DTW_files->file_kmin = (char *) malloc( 1024 );
	DTW_files->file_states_time_e = (char *) malloc( 1024 );
	DTW_files->file_states_time_i = (char *) malloc( 1024 );
	DTW_files->file_states_seq = (char *) malloc( 1024 );

#if VERIFY
	DTW_verify->file_validation_hanning = (char *) malloc( 1024 );
	DTW_verify->file_validation_fft = (char *) malloc( 1024 );
	DTW_verify->file_validation_vec_dist = (char *) malloc( 1024 );
	DTW_verify->file_validation_distor = (char *) malloc( 1024 );
#endif

	return EXIT_SUCCESS;
}

