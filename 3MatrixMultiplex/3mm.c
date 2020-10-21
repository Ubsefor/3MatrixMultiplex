    //
    //  3mm.c
    //  3MatrixMultiplex
    //
    //  Created by Alexander Makhov on 20/10/20.
    //

#include "3mm.h"

void bench_timer_start( void )
{
    bench_t_start = rtclock();
}

void bench_timer_stop( void )
{
    bench_t_end = rtclock();
}

void bench_timer_print( void )
{
    printf( "Time in seconds = %0.6lf\n", bench_t_end - bench_t_start );
}


void init_array( int ni, int nj, int nk, int nl, int nm,
                       float A[ni][nk], float B[nk][nj],
                       float C[nj][nm], float D[nm][nl] )
{
    int i, j;
    omp_set_num_threads( THREAD_NUM );
#pragma omp parallel shared(A,B,C,D) private(i,j)
    {
#pragma omp for  schedule(static)
        for ( i = 0; i < ni; i++ )
        for ( j = 0; j < nk; j++ )
        A[i][j] = (float) ( ( i * j + 1 ) % ni ) / ( 5 * ni );
        
#pragma omp for  schedule(static)
        for ( i = 0; i < nk; i++ )
        for ( j = 0; j < nj; j++ )
        B[i][j] = (float) ( ( i * ( j + 1 ) + 2 ) % nj ) / ( 5 * nj );
        
#pragma omp for  schedule(static)
        for ( i = 0; i < nj; i++ )
        for ( j = 0; j < nm; j++ )
        C[i][j] = (float) ( i * ( j + 3 ) % nl ) / ( 5 * nl );
        
#pragma omp for  schedule(static)
        for ( i = 0; i < nm; i++ )
        for ( j = 0; j < nl; j++ )
        D[i][j] = (float) ( ( i * ( j + 2 ) + 2 ) % nk ) / ( 5 * nk );
    }
}

void print_array( int ni, int nl, float G[ni][nl] )
{
    int i, j;
    
    fprintf( stderr, "==BEGIN DUMP_ARRAYS==\n" );
    fprintf( stderr, "begin dump: %s", "G" );
    for ( i = 0; i < ni; i++ )
    for ( j = 0; j < nl; j++ )
    {
        if ( ( i * ni + j ) % 20 == 0 )
            fprintf( stderr, "\n" );
        fprintf( stderr, "%0.2f ", G[i][j] );
    }
    fprintf( stderr, "\nend   dump: %s\n", "G" );
    fprintf( stderr, "==END   DUMP_ARRAYS==\n" );
}


void kernel_3mm( int ni, int nj, int nk, int nl, int nm,
                       float E[ni][nj], float A[ni][nk], float B[nk][nj],
                       float F[nj][nl], float C[nj][nm], float D[nm][nl], float G[ni][nl] )
{
    int i = 0, j = 0, k = 0;
    omp_set_num_threads( THREAD_NUM );
#pragma omp parallel shared(E,A,B,F,C,D,G) private(i,j,k)
    {
#pragma omp for  schedule(static)
        for ( i = 0; i < ni; i++ )
        for ( j = 0; j < nj; j++ )
        {
            E[i][j] = 0.0f;
            for ( k = 0; k < nk; ++k )
            E[i][j] += A[i][k] * B[k][j];
        }
#pragma omp for  schedule(static)
        for ( i = 0; i < nj; i++ )
        for ( j = 0; j < nl; j++ )
        {
            F[i][j] = 0.0f;
            for ( k = 0; k < nm; ++k )
            F[i][j] += C[i][k] * D[k][j];
        }
        
#pragma omp for  schedule(static)
        for ( i = 0; i < ni; i++ )
        for ( j = 0; j < nl; j++ )
        {
            G[i][j] = 0.0f;
            for ( k = 0; k < nj; ++k )
            G[i][j] += E[i][k] * F[k][j];
        }
    }
}

int main( int argc, char** argv )
{
    int ni = NI;
    int nj = NJ;
    int nk = NK;
    int nl = NL;
    int nm = NM;
    
    float (*E)[ni][nj]; E = ( float(*)[ni][nj] )malloc( ( ni ) * ( nj ) * sizeof ( float ) );
    float (*A)[ni][nk]; A = ( float(*)[ni][nk] )malloc( ( ni ) * ( nk ) * sizeof ( float ) );
    float (*B)[nk][nj]; B = ( float(*)[nk][nj] )malloc( ( nk ) * ( nj ) * sizeof ( float ) );
    float (*F)[nj][nl]; F = ( float(*)[nj][nl] )malloc( ( nj ) * ( nl ) * sizeof ( float ) );
    float (*C)[nj][nm]; C = ( float(*)[nj][nm] )malloc( ( nj ) * ( nm ) * sizeof ( float ) );
    float (*D)[nm][nl]; D = ( float(*)[nm][nl] )malloc( ( nm ) * ( nl ) * sizeof ( float ) );
    float (*G)[ni][nl]; G = ( float(*)[ni][nl] )malloc( ( ni ) * ( nl ) * sizeof ( float ) );
    
    init_array( ni, nj, nk, nl, nm,
               *A,
               *B,
               *C,
               *D );
    
    bench_timer_start();
    
    kernel_3mm( ni, nj, nk, nl, nm,
               *E,
               *A,
               *B,
               *F,
               *C,
               *D,
               *G );
    
    bench_timer_stop();
    bench_timer_print();
    
    if ( argc > 42 && !strcmp( argv[0], "" ) )
        print_array( ni, nl, *G );
    
    free( (void*) E );
    free( (void*) A );
    free( (void*) B );
    free( (void*) F );
    free( (void*) C );
    free( (void*) D );
    free( (void*) G );
    
    return 0;
}
