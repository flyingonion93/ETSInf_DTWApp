/*! \file DTWProcess.c
 *	dtw
 *
 *		\brief		Processing functions of the application
 *		\author		Carlos Gómez Morillas
 *		\version	1.0
 *		\date		25/7/16
 *		\copyright	GNU Public License
 *
 *		Last modification by Carlos Gómez Morillas on 31/5/16
 */
#include "DTWProcess.h"


/*
 * Private module interface
 */
static int Compute_minimum_value( const int N, const precission_type *v );

static int DTW_compute( const precission_type *sequence, const precission_type *costs, precission_type * __restrict pD, precission_type * __restrict pV, int * __restrict pS, const int NST, const int NCST, const int TB, int *TBDes );

/*!
 *	\brief User function that computes the DTW algorithm to obtain the minimum position
 *
 *	\param [in] DTW_const		Common values structure
 */	
int DTW_process( DTW_const_st DTW_const )
{
	minimum_position = Compute_minimum_value( DTW_size - 1, &v_SxD[1] );
	if( 0 != minimum_position )
		minimum_position = DTW_compute( &v_SxD[1], costs, pD, pV, pS, DTW_size - 1, DTW_const.nc, DTW_const.tb, &TBDes );

	return EXIT_SUCCESS;
}

/*!
 *	\brief Returns the position of the minimum value in a vector
 *
 *	\param [in] N			Size of the vector
 *	\param [in] v			Vector that we need to loop
 */
static int Compute_minimum_value( const int N, const precission_type *v )
{
	int i;
	int minimum_position = 0;

	precission_type minimum_value = v[0];

	#pragma omp parallel
	{
#if SIMPLE_PRECISSION
		precission_type value = FLT_MAX;      
#else
		precission_type value = DBL_MAX;
#endif     
		int position = 1;
		
		#pragma omp for nowait
		for( i = 1; i < N; i++ )
		{
			if( v[i] < value )
			{
				value = v[i];
				position = i;
			}
		}
		if( value < minimum_value )
		{
			#pragma omp critical (ZonaCritica1)
			{
				if( value < minimum_value )
				{
					minimum_value = value;
					minimum_position = position;
				}
			}
		}
	}
	return minimum_position;
}

/*!
 *	\brief DTW algorithm function.
 *
 *	\param [in] sequence		Parameter description
 *	\param [out] pD			Vector that stores the distance between frames
 *	\param [out] pV			Vector that stores the quotient between frame distance and number of jumps
 *	\param [out] pS			Vector that stores the number of jumps
 *	\param [in] NST			Number of states of a frame
 *	\param [in] NCST		Number of DTW dependencies
 *	\param [in] TB			Extended buffer size
 *	\param [out] TBDes		Parameter description
 */
static int DTW_compute(const precission_type *sequence, const precission_type *costs, precission_type * __restrict pD, precission_type * __restrict pV,
            int * __restrict pS, const int NST, const int NCST, const int TB, int *TBDes)
{
	int NST_plus_NC = NST + NCST;
	int back_step = ( NST_plus_NC * NCST ) + ( (*TBDes) * NST_plus_NC ) - NST - 1;
	int relative_position = back_step + NST_plus_NC + 1;
	int j;

	#pragma omp parallel for
	for( j = ( NST - 1 ); j >= 0; j-- )
	{
		precission_type d;
		precission_type d2;
		precission_type v;
		precission_type v2;
		
		int s = 0;
		int s2;
		int k;
		int position;

#if SIMPLE_PRECISSION
		d = FLT_MAX;
		v = FLT_MAX; 
#else
		d = DBL_MAX;
		v = DBL_MAX;
#endif
	
		position = back_step + j;

		for( k = 0; k < NCST; k++ )
		{
			d2 = sequence[j] * costs[k] + pD[position-k];
			s2 = 1  + pS[position-k];
			v2 = d2 / s2;
			
			if( v2 < v )
			{
				d = d2;
				v = v2;
				s = s2;
			}
		}

		for( k = NCST; k < 2 * NCST - 1; k++ )
		{
			position = position - NST_plus_NC;
			d2 = sequence[j] * costs[k] + pD[position];
			s2 = 1 + pS[position];
			v2 = d2 / s2;

			if( v2 < v )
			{
				d = d2;
				v = v2;
				s = s2;
			}
		}
		
		pD[relative_position + j] = d;
		pS[relative_position + j] = s;
		pV[relative_position + j] = v;
	}

	j = Compute_minimum_value( NST, &pV[relative_position] );

	(*TBDes)++;
	if( (*TBDes) == TB)
	{
		memmove( pD, &pD[NST_plus_NC * TB], sizeof(precission_type) * NST_plus_NC * NCST );
		memmove( pV, &pV[NST_plus_NC * TB], sizeof(precission_type) * NST_plus_NC * NCST );
     		memmove( pS, &pS[NST_plus_NC * TB], sizeof(int) * NST_plus_NC * NCST );
		(*TBDes) = 0;
	}

	return j;
}
