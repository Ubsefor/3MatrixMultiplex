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
#define TIMEOUT       10000000

#define MASTER 0               /* task id of first task  */
#define FROM_MASTER 1          /* setting a message type */
#define FROM_WORKER 2          /* setting a message type */

#ifndef _3MM_H
#define _3MM_H

# if !defined ( MINI_DATASET ) && !defined ( SMALL_DATASET ) && !defined ( MEDIUM_DATASET ) && !defined ( LARGE_DATASET ) && !defined ( EXTRALARGE_DATASET )
#define MEDIUM_DATASET
# endif

# if !defined ( NI ) && !defined ( NJ ) && !defined ( NK ) && !defined ( NL ) && !defined ( NM )
# ifdef MINI_DATASET
#define NI    32
#define NJ    36
#define NK    40
#define NL    44
#define NM    48
# endif
# ifdef SMALL_DATASET
#define NI    160
#define NJ    200
#define NK    240
#define NL    280
#define NM    320
# endif
# ifdef MEDIUM_DATASET
#define NI    720
#define NJ    760
#define NK    800
#define NL    840
#define NM    880
# endif
# ifdef LARGE_DATASET
#define NI    1600
#define NJ    1800
#define NK    2000
#define NL    2200
#define NM    2400
# endif
# ifdef EXTRALARGE_DATASET
#define NI    3200
#define NJ    3600
#define NK    4000
#define NL    4400
#define NM    4800
# endif

#endif

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
#include <errno.h>
#include <time.h>
#include <sys/time.h>
#include <signal.h>
#include <mpi.h>
#include <mpi-ext.h>
#include <stdarg.h>
#include <stdbool.h>

#endif

double
    bench_t_start, bench_t_end,
    mpi_t_start, mpi_t_end;

int rc; /* error code from MPI functions */
int victim;
int spare;

bool setdie = true;

char estr[MPI_MAX_ERROR_STRING]=""; int strl; /* error messages */

MPI_Comm world; /* a world comm for the work, w/o the spares */
MPI_Comm rworld; /* and a temporary handle to store the repaired copy */
MPI_Status status;
MPI_Request request;

int stage = 0;

char** gargv;

int ping = 0;

int ni = NI;
int nj = NJ;
int nk = NK;
int nl = NL;
int nm = NM;

double A[NI][NK];
double B[NK][NJ];
double E[NI][NJ] = {0};

double C[NJ][NM];
double D[NM][NL];
double F[NJ][NL] = {0};

double G[NI][NL] = {0};

int
    np,              /* number of tasks in partition */
    rank,                /* a task identifier */
    numworkers,            /* number of worker tasks */
    source,                /* task id of message source */
    dest,                  /* task id of message destination */
    mtype,                 /* message type */
    rows,                  /* rows of matrix A sent to each worker */
    averow, extra, offset, /* used to determine rows sent to each worker */
    i, j, k, rc;           /* misc */

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
void print_array( int ni, int nl, double G[ni][nl] );

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
                double E[ni][nj], double A[ni][nk], double B[nk][nj],
                double F[nj][nl], double C[nj][nm], double D[nm][nl],
                double G[ni][nl] );

int validate_param( char *a );

#endif /* _mm_h */
