    //
    //  3mm.c
    //  3MatrixMultiplex
    //
    //  Created by Alexander Makhov on 20/10/20.
    //

#include "3mm.h"

int validate_param( char *a )
{
    unsigned x;
    for ( x = 0; x < strlen( a ); x++ )
    if ( !isdigit( a[x] ) )
        return 1;
    return 0;
}

void MPI_timer_start( void ){
    mpi_t_start = MPI_Wtime();
}

void MPI_timer_end( void ){
    mpi_t_end = MPI_Wtime();
}

void MPI_timer_print( void ){
#ifndef BENCH
    printf("time =  %lf seconds\n", mpi_t_end - mpi_t_start);
#else
    printf("%0.6lf\n", mpi_t_end - mpi_t_start);
#endif
}

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
#ifndef BENCH
    printf( "Time in seconds = %0.6lf\n", bench_t_end - bench_t_start );
#else
    printf( "%0.6lf\n", bench_t_end - bench_t_start );
#endif
}

void init_array( int ni, int nj, int nk, int nl, int nm,
                float A[ni][nk], float B[nk][nj],
                float C[nj][nm], float D[nm][nl] )
{
    for ( int i = 0; i < ni; i++ )
        for ( int j = 0; j < nk; j++ )
            A[i][j] = (float) ( ( i * j + 1 ) % ni ) / ( 5 * ni );

    for ( int i = 0; i < nk; i++ )
        for ( int j = 0; j < nj; j++ )
            B[i][j] = (float) ( ( i * ( j + 1 ) + 2 ) % nj ) / ( 5 * nj );

    for ( int i = 0; i < nj; i++ )
        for ( int j = 0; j < nm; j++ )
            C[i][j] = (float) ( i * ( j + 3 ) % nl ) / ( 5 * nl );
    
    for ( int i = 0; i < nm; i++ )
        for ( int j = 0; j < nl; j++ )
            D[i][j] = (float) ( ( i * ( j + 2 ) + 2 ) % nk ) / ( 5 * nk );
}

void print_array( int ni, int nl, double G[ni][nl] )
{
    fprintf( stderr, "==BEGIN DUMP_ARRAYS==\n" );
    fprintf( stderr, "begin dump: %s", "G" );
    for ( int i = 0; i < ni; i++ )
        for ( int j = 0; j < nl; j++ )
        {
            if ( ( i * ni + j ) % 20 == 0 )
                fprintf( stderr, "\n" );
            fprintf( stderr, "%0.2lf ", G[i][j] );
        }
    fprintf( stderr, "\nend   dump: %s\n", "G" );
    fprintf( stderr, "==END   DUMP_ARRAYS==\n" );
}

void kernel_3mm( int ni, int nj, int nk, int nl, int nm,
                double E[ni][nj], double A[ni][nk], double B[nk][nj],
                double F[nj][nl], double C[nj][nm], double D[nm][nl], double G[ni][nl] )
{
    for ( int i = 0; i < ni; i++ )
        for ( int j = 0; j < nj; j++ )
        {
            E[i][j] = 0.0f;
            for ( int k = 0; k < nk; ++k )
                E[i][j] += A[i][k] * B[k][j];
        }
    for ( int i = 0; i < nj; i++ )
        for ( int j = 0; j < nl; j++ )
        {
            F[i][j] = 0.0f;
            for ( int k = 0; k < nm; ++k )
                F[i][j] += C[i][k] * D[k][j];
        }
    for ( int i = 0; i < ni; i++ )
        for ( int j = 0; j < nl; j++ )
        {
            G[i][j] = 0.0f;
            for ( int k = 0; k < nj; ++k )
                G[i][j] += E[i][k] * F[k][j];
        }
}

static void MPI_init_matrixes() {
    if (taskid == MASTER){
        for (int i=0; i< ni; i++)
        for (int j=0; j < nj; j++)
        A[i][j] = (double) ( ( i * j + 1 ) % ni ) / ( 5 * ni );
        
        for (int i=0; i< ni; i++)
        for (int j=0; j < nj; j++)
        B[i][j] = (double) ( ( i * ( j + 1 ) + 2 ) % nj ) / ( 5 * nj );
        
        for ( int i = 0; i < nj; i++ )
        for ( int j = 0; j < nm; j++ )
        C[i][j] = (double) ( i * ( j + 3 ) % nl ) / ( 5 * nl );
        
        for ( int i = 0; i < nm; i++ )
        for ( int j = 0; j < nl; j++ )
        D[i][j] = (double) ( ( i * ( j + 2 ) + 2 ) % nk ) / ( 5 * nk );
        
#ifndef BENCH
        printf("Init successful\n");
#endif
    }
}

/*
 float A[ni][nk] *
 float B[nk][nj] =
 float E[ni][nj]
 
 float C[nj][nm] *
 float D[nm][nl] =
 float F[nj][nl]
 
 float E[ni][nj] *
 float F[nj][nl] =
 float G[ni][nl]
 */
static void MPI_kernel_3mm() {
    if (taskid == MASTER)
    {
#ifndef BENCH
        printf("MPI_kernel_3mm has started with %d tasks.\n",numtasks);
#endif
        
        /* Measure start time */
        MPI_timer_start();
        
        /* Send A, B matrix data to the worker tasks */
        averow = NI / numworkers;
        extra = NI % numworkers;
        offset = 0;
        mtype = FROM_MASTER;
        
        for (int dest=1; dest <= numworkers; dest++)
        {
                // Master sending data
            rows = (dest <= extra) ? averow+1 : averow;
#ifndef BENCH
            printf("Sending %d rows of A and B to task %d offset=%d\n",rows,dest,offset);
#endif
            MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
#ifndef BENCH
            printf("Offset sent to %d\n",dest);
#endif
            MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
#ifndef BENCH
            printf("Rows sent to %d\n",dest);
#endif
            MPI_Send(&A[offset][0], rows * NK, MPI_DOUBLE, dest, mtype,
                     MPI_COMM_WORLD);
#ifndef BENCH
            printf("A sent to %d\n",dest);
#endif
            MPI_Send(&B, NK * NJ, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
#ifndef BENCH
            printf("B sent to %d\n",dest);
#endif
            offset = offset + rows;
        }
        
        mtype = FROM_WORKER;
        for (int i = 1; i <= numworkers; i++)
        {
                // Master collecting data
            source = i;
#ifndef BENCH
            printf("%d: Requesting data from %d task\n", taskid, i);
#endif
            MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
#ifndef BENCH
            printf("%d: Received offset from task %d for matrix E\n", taskid, source);
#endif
            MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
#ifndef BENCH
            printf("%d: Received rows from task %d for matrix E\n", taskid, source);
#endif
            MPI_Recv(&E[offset][0], rows * NJ, MPI_DOUBLE, source, mtype,
                     MPI_COMM_WORLD, &status);
#ifndef BENCH
            printf("Received results from task %d for matrix E, waiting for %d tasks\n",source, numworkers-i);
#endif
        }
        
#ifndef BENCH
        printf ("Awaiting sync...\n");
#endif
        MPI_Barrier(MPI_COMM_WORLD);
        
        /* Send C, D matrix data to the worker tasks */
        averow = NJ / numworkers;
        extra = NJ % numworkers;
        offset = 0;
        mtype = FROM_MASTER;
        
        for (int dest=1; dest <= numworkers; dest++)
        {
                // Master sending data
            rows = (dest <= extra) ? averow+1 : averow;
#ifndef BENCH
            printf("Sending %d rows of C and D to task %d offset=%d\n", rows, dest, offset);
#endif
            MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
#ifndef BENCH
            printf("Offset sent to %d\n", dest);
#endif
            MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
#ifndef BENCH
            printf("Rows sent to %d\n",dest);
#endif
            MPI_Send(&C[offset][0], rows * NM, MPI_DOUBLE, dest, mtype,
                     MPI_COMM_WORLD);
#ifndef BENCH
            printf("C sent to %d\n",dest);
#endif
            MPI_Send(&D, NM * NL, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
#ifndef BENCH
            printf("D sent to %d\n",dest);
#endif
            offset = offset + rows;
        }
        
        mtype = FROM_WORKER;
        for (int i = 1; i <= numworkers; i++)
        {
                // Master collecting data
            source = i;
#ifndef BENCH
            printf("%d: Requesting data from %d task\n", taskid, i);
#endif
            MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
#ifndef BENCH
            printf("%d: Received offset from task %d for matrix E\n", taskid, source);
#endif
            MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
#ifndef BENCH
            printf("%d: Received rows from task %d for matrix E\n", taskid, source);
#endif
            MPI_Recv(&F[offset][0], rows * NL, MPI_DOUBLE, source, mtype,
                     MPI_COMM_WORLD, &status);
#ifndef BENCH
            printf("Received results from task %d for matrix F\n",source);
#endif
        }
        
#ifndef BENCH
        printf ("Awaiting sync...\n");
#endif
        MPI_Barrier(MPI_COMM_WORLD);
        
        /* Send E, F matrix data to the worker tasks */
        averow = NI / numworkers;
        extra = NI % numworkers;
        offset = 0;
        mtype = FROM_MASTER;
        
        for (int dest=1; dest <= numworkers; dest++)
        {
                // Master sending data
            rows = (dest <= extra) ? averow+1 : averow;
#ifndef BENCH
            printf("Sending %d rows of E and F to task %d offset=%d\n",rows,dest,offset);
#endif
            MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
#ifndef BENCH
            printf("Offset sent to %d\n", dest);
#endif
            MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
#ifndef BENCH
            printf("Rows sent to %d\n",dest);
#endif
            MPI_Send(&E[offset][0], rows * NJ, MPI_DOUBLE, dest, mtype,
                     MPI_COMM_WORLD);
#ifndef BENCH
            printf("E sent to %d\n",dest);
#endif
            MPI_Send(&F, NJ * NL, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
#ifndef BENCH
            printf("F sent to %d\n",dest);
#endif
            offset = offset + rows;
            
        }
        
        mtype = FROM_WORKER;
        for (int i = 1; i <= numworkers; i++)
        {
                // Master collecting data
            source = i;
#ifndef BENCH
            printf("%d: Requesting data from %d task\n", taskid, i);
#endif
            MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
#ifndef BENCH
            printf("%d: Received offset from task %d for matrix E\n", taskid, source);
#endif
            MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
#ifndef BENCH
            printf("%d: Received rows from task %d for matrix E\n", taskid, source);
#endif
            MPI_Recv(&G[offset][0], rows * NL, MPI_DOUBLE, source, mtype,
                     MPI_COMM_WORLD, &status);
#ifndef BENCH
            printf("Received results from task %d\n for matrix G\n",source);
#endif
        }
        
#ifndef BENCH
        printf ("Awaiting sync...\n");
#endif
        
        MPI_Barrier(MPI_COMM_WORLD);
#ifndef BENCH
        printf("%d: All workers exited\n", taskid);
#endif
        
        MPI_timer_end();
        MPI_timer_print();
        
    }
    
    // Workers
    
    if (taskid > MASTER)
    {
            // Workers calculating E
        {
#ifndef BENCH
            printf("E calculation in %d task\n", taskid);
#endif
            
            mtype = FROM_MASTER;
#ifndef BENCH
            printf("%d: Requesting data from MASTER\n", taskid);
#endif
            MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
#ifndef BENCH
            printf("%d: Got offset from MASTER\n", taskid);
#endif
            MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
#ifndef BENCH
            printf("%d: Got rows from MASTER\n", taskid);
#endif
            MPI_Recv(&A, rows * NK, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
#ifndef BENCH
            printf("%d: Got first matrix from MASTER\n", taskid);
#endif
            MPI_Recv(&B, NK * NJ, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
#ifndef BENCH
            printf("%d: Got second matrix from MASTER\n", taskid);
            printf("Received all data for task %d\n", taskid);
#endif
            
            for (int k = 0; k < NJ; k++)
                for (int i = 0; i < rows; i++)
                {
                    E[i][k] = 0.0;
                    for (int j = 0; j < NI; j++)
                        E[i][k] = E[i][k] + A[i][j] * B[j][k];
                }
            
            mtype = FROM_WORKER;
            
#ifndef BENCH
            printf("Result: sending data from %d task to MASTER\n", taskid);
#endif
            MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
#ifndef BENCH
            printf("Result: send offset from %d task to MASTER\n", taskid);
#endif
            MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
#ifndef BENCH
            printf("Result: sent rows from %d task to MASTER\n", taskid);
#endif
            MPI_Send(&E, rows * NJ, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
#ifndef BENCH
            printf("Successfully sent all data from %d task to MASTER\n", taskid);
            printf("Awaiting sync...\n");
#endif
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
            // Workers calculating F
        
        {
            mtype = FROM_MASTER;
#ifndef BENCH
            printf("%d: Requesting data from MASTER\n", taskid);
#endif
            MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
#ifndef BENCH
            printf("%d: Got offset from MASTER\n", taskid);
#endif
            MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
#ifndef BENCH
            printf("%d: Got rows from MASTER\n", taskid);
#endif
            MPI_Recv(&C, rows * NM, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
#ifndef BENCH
            printf("%d: Got first matrix from MASTER\n", taskid);
#endif
            MPI_Recv(&D, NM * NL, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
#ifndef BENCH
            printf("%d: Got second matrix from MASTER\n", taskid);
            printf("Received all data for task %d\n", taskid);
#endif
            
            for (int k = 0; k < NL; k++)
                for (int i = 0; i < rows; i++)
                {
                    F[i][k] = 0.0;
                    for (int j = 0; j < NJ; j++)
                        F[i][k] = F[i][k] + C[i][j] * D[j][k];
                }
            mtype = FROM_WORKER;
#ifndef BENCH
            printf("Result: sending data from %d task to MASTER\n", taskid);
#endif
            MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
#ifndef BENCH
            printf("Result: sent offset from %d task to MASTER\n", taskid);
#endif
            MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
#ifndef BENCH
            printf("Result: sent rows from %d task to MASTER\n", taskid);
#endif
            MPI_Send(&F, rows * NL, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
#ifndef BENCH
            printf("Successfully sent all data from %d task to MASTER\n", taskid);
            printf("Awaiting sync...\n");
#endif
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
            // Workers calculating G
        {
            mtype = FROM_MASTER;
#ifndef BENCH
            printf("%d: Requesting data from MASTER\n", taskid);
#endif
            MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
#ifndef BENCH
            printf("%d: Got offset from MASTER\n", taskid);
#endif
            MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
#ifndef BENCH
            printf("%d: Got rows from MASTER\n", taskid);
#endif
            MPI_Recv(&E, rows * NJ, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
#ifndef BENCH
            printf("%d: Got first matrix from MASTER\n", taskid);
#endif
            MPI_Recv(&F, NJ * NL, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
#ifndef BENCH
            printf("%d: Got second matrix from MASTER\n", taskid);
            printf("Received all data for task %d\n", taskid);
#endif
            
            for (k=0; k < NL; k++)
                for (i=0; i < rows; i++)
                {
                    G[i][k] = 0.0;
                    for (j=0; j < NI; j++)
                        G[i][k] = G[i][k] + E[i][j] * F[j][k];
                }
            
            mtype = FROM_WORKER;
#ifndef BENCH
            printf("Result: sending data from %d task to MASTER\n", taskid);
#endif
            MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
#ifndef BENCH
            printf("Result: sent offset from %d task to MASTER\n", taskid);
#endif
            MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
#ifndef BENCH
            printf("Result: sent rows from %d task to MASTER\n", taskid);
#endif
            MPI_Send(&G, rows * NL, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
        }
#ifndef BENCH
        printf("Successfully sent all data from %d task to MASTER\n", taskid);
        printf("Awaiting sync...\n");
#endif
        
        MPI_Barrier(MPI_COMM_WORLD);
        
#ifndef BENCH
        printf("%d: Done working.\n", taskid);
#endif
        
    }
}

int main( int argc, char** argv )
{
        
#ifndef BENCH
#ifdef MINI_DATASET
    printf("Mini dataset selected...\n");
#endif
#ifdef SMALL_DATASET
    printf("Small dataset selected...\n");
#endif
#ifdef MEDIUM_DATASET
    printf("Medium dataset selected...\n");
#endif
#ifdef LARGE_DATASET
    printf("Large dataset selected...\n");
#endif
#ifdef EXTRALARGE_DATASET
    printf("Extra dataset selected...\n");
#endif
#endif
   
#ifndef BENCH
    printf("Setting defines successfull\n");
#endif

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
    
    if (numtasks < 2 ) {
        bench_timer_start();
        kernel_3mm( ni, nj, nk, nl, nm,
                   E,
                   A,
                   B,
                   F,
                   C,
                   D,
                   G );
            
        bench_timer_stop();
        bench_timer_print();
        MPI_Finalize();
        return 0;
    }
    numworkers = numtasks-1;
    
    MPI_init_matrixes();
    
    
#ifndef BENCH
    if (taskid == MASTER)
        printf("Tasks comm port set successful\n");
#endif
    
    MPI_kernel_3mm();
    
    if ( (taskid == MASTER) && argc > 42 && !strcmp( argv[0], "" ) )
        print_array( ni, nl, G );
    
#ifndef BENCH
    if (taskid == MASTER)
        printf( "Just finished!\n" );
#endif
    MPI_Finalize();
    return 0;
}

