    //
    //  3mm.h
    //  3MatrixMultiplex
    //
    //  Created by Alexander Makhov on 20/10/20.
    //


/*!
 \file 3mm.h
 \brief Header file with function definitions
 
 */

#ifndef _mm_h
#define _mm_h

#define THREAD_NUM    16

#ifndef _3MM_H
#define _3MM_H
# if !defined ( MINI_DATASET ) && !defined ( SMALL_DATASET ) && !defined ( MEDIUM_DATASET ) && !defined ( LARGE_DATASET ) && !defined ( EXTRALARGE_DATASET )
#define LARGE_DATASET
# endif
# if !defined ( NI ) && !defined ( NJ ) && !defined ( NK ) && !defined ( NL ) && !defined ( NM )
# ifdef MINI_DATASET
#define NI    16
#define NJ    18
#define NK    20
#define NL    22
#define NM    24
# endif
# ifdef SMALL_DATASET
#define NI    40
#define NJ    50
#define NK    60
#define NL    70
#define NM    80
# endif
# ifdef MEDIUM_DATASET
#define NI    180
#define NJ    190
#define NK    200
#define NL    210
#define NM    220
# endif
# ifdef LARGE_DATASET
#define NI    800
#define NJ    900
#define NK    1000
#define NL    1100
#define NM    1200
# endif
# ifdef EXTRALARGE_DATASET
#define NI    1600
#define NJ    1800
#define NK    2000
#define NL    2200
#define NM    2400
# endif
#endif
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>

#endif

double bench_t_start, bench_t_end;

double rtclock( void )
{
    struct timeval Tp;
    int            stat;
    stat = gettimeofday( &Tp, NULL );
    if ( stat != 0 )
        printf( "Error return from gettimeofday: %d", stat );
    return ( Tp.tv_sec + Tp.tv_usec * 1.0e-6 );
}

    //! Benchmark timer starter
void bench_timer_start( void );

    //! Benchmark timer stopper
void bench_timer_stop( void );

    //! Benchmark timer printer
void bench_timer_print( void );

/*!
     Initializes arrays for matrixes
     @author Ubsefor
     @version 1.0.1
     @param ni Lenght of matrix A
     @param nk Height of matrix A, Length of matrix B
     @param nj Height of matrix B, Length of matrix C
     @param nm Height of matrix C, Length of matrix D
     @param nl Height of matrix D
     @param A A Matrix to initialize
     @param B B Matrix to initialize
     @param C C Matrix to initialize
     @param D D Matrix to initialize
*/
void init_array( int ni, int nj, int nk, int nl, int nm,
                       float A[ni][nk], float B[nk][nj],
                       float C[nj][nm], float D[nm][nl] );
/*!
     Dumps arrays into STDIN
     @author Ubsefor
     @version 1.0.0
     @warning For debug only
     @param ni Matrix length
     @param nl Matrix height
     @param G Matrix itself
*/
void print_array( int ni, int nl, float G[ni][nl] );

/*!
     Kernel for multiplication of matrixes
     @author Ubsefor
     @version 1.0.1
     @param ni Length of matrixes E, A, G
     @param nj Height of matrixes E, B; Length of matrixes F, C
     @param nk Height of matrix A; Lenght of matrix B
     @param nl Height of matrixes F, D, G
     @param nm Height of matrix C; Length of matrix D
     @param A Matrix A pre inited
     @param B Matrix B pre inited
     @param C Matrix C pre inited
     @param D Matrix D pre inited
     @param E Resultig matrix of A*B
     @param F Resulting matrix of C*D
     @param G Resulting matrix of E*F
*/
void kernel_3mm( int ni, int nj, int nk, int nl, int nm,
                       float E[ni][nj], float A[ni][nk], float B[nk][nj],
                       float F[nj][nl], float C[nj][nm], float D[nm][nl],
                       float G[ni][nl] );

#endif /* _mm_h */
