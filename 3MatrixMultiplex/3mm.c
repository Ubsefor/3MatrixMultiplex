    //
    //  3mm.c
    //  3MatrixMultiplex
    //
    //  Created by Alexander Makhov on 20/10/20.
    //

#include "3mm.h"

// microsleep in nanoseconds
int msleep(long msec)
{
    struct timespec ts;
    int res;

    if (msec < 0)
    {
        errno = EINVAL;
        return -1;
    }

    ts.tv_sec = msec / 1000;
    ts.tv_nsec = (msec % 1000) * 1000000;

    do {
        res = nanosleep(&ts, &ts);
    } while (res && errno == EINTR);

    return res;
}

// waits for timeout in ms, then calls abort on operations
static int MPIX_Wait_timeout(MPI_Request *request, long timeout){
    long tries = 0;
    int flag = 0;
    while (!flag) {
        if (tries > timeout) {
            MPIX_Comm_revoke(world);
            return MPI_ERR_PROC_FAILED_PENDING;
        }
        MPI_Test(request, &flag, &status);
        msleep(10);
        tries++;
    }
    return MPI_SUCCESS;
}

int MPIX_Comm_replace(MPI_Comm comm, MPI_Comm *newcomm) {
    MPI_Comm icomm, /* the intercomm between the spawnees and the old (shrinked) world */
             scomm, /* the local comm for each sides of icomm */
             mcomm; /* the intracomm, merged from icomm */
    MPI_Group cgrp, sgrp, dgrp;
    int rc, flag, rflag, i, nc, ns, nd, crank, srank, drank;

redo:
    if( comm == MPI_COMM_NULL ) { /* am I a new process? */
        /* I am a new spawnee, waiting for my new rank assignment
         * it will be sent by rank 0 in the old world */
        MPI_Comm_get_parent(&icomm);
        scomm = MPI_COMM_WORLD;
        // receiving rank
        MPI_Recv(&crank, 1, MPI_INT, 0, 1, icomm, MPI_STATUS_IGNORE);
        // receiving previous stage
        MPI_Recv(&stage, sizeof(int), MPI_INT, 0, 1, icomm, MPI_STATUS_IGNORE);
        MPI_Comm_rank(scomm, &srank);
        printf("Spawnee %d: crank=%d\n", srank, crank);
    }
    else {
        /* I am a survivor: Spawn the appropriate number
         * of replacement processes (we check that this operation worked
         * before we procees further) */
        /* First: remove dead processes */
        MPIX_Comm_shrink(comm, &scomm);
        MPI_Comm_size(scomm, &ns);
        MPI_Comm_size(comm, &nc);
        nd = nc-ns; /* number of deads */
        printf("Report survivor: %d\n", rank);
        if( 0 == nd ) {
            /* Nobody was dead to start with. We are done here */
            MPI_Comm_free(&scomm);
            *newcomm = comm;
            return MPI_SUCCESS;
        }
        /* We handle failures during this function ourselves... */
        MPI_Comm_set_errhandler( scomm, MPI_ERRORS_RETURN );
        
        
        rc = MPI_Comm_spawn(gargv[0], &gargv[1], nd, MPI_INFO_NULL,
                            0, scomm, &icomm, MPI_ERRCODES_IGNORE);
        printf("%d: Tried spawning a spare process\n", rank);
        
        flag = (MPI_SUCCESS == rc);
        MPIX_Comm_agree(scomm, &flag);
        if( !flag ) {
            if( MPI_SUCCESS == rc ) {
                MPIX_Comm_revoke(icomm);
                MPI_Comm_free(&icomm);
            }
            MPI_Comm_free(&scomm);
            fprintf(stderr, "%d: comm_spawn failed, redo\n", rank);
            goto redo;
        }

        
        /* remembering the former rank: we will reassign the same
         * ranks in the new world. */
        MPI_Comm_rank(comm, &crank);
        MPI_Comm_rank(scomm, &srank);
        printf("%d Got former rank\n", srank);
        /* the rank 0 in the scomm comm is going to determine the
         * ranks at which the spares need to be inserted. */
        if(0 == srank) {
            /* getting the group of dead processes:
             *   those in comm, but not in scomm are the deads */
            printf("Root reporting!\n");
            MPI_Comm_group(comm, &cgrp);
            MPI_Comm_group(scomm, &sgrp);
            MPI_Group_difference(cgrp, sgrp, &dgrp);
            /* Computing the rank assignment for the newly inserted spares */
            for(i=0; i<nd; i++) {
                MPI_Group_translate_ranks(dgrp, 1, &i, cgrp, &drank);
                /* sending their new assignment to all new procs */
                MPI_Send(&drank, 1, MPI_INT, i, 1, icomm);
                printf("%d: Sent new rank assignment to spare\n", srank);
                /* sending previous stage */
                MPI_Send(&stage, sizeof(int), MPI_INT, i, 1, icomm);
                printf("%d: Sent stage to spare\n", srank);
            }
            MPI_Group_free(&cgrp); MPI_Group_free(&sgrp); MPI_Group_free(&dgrp);
            printf("%d: Freed old group\n", srank);
        }
        
    }

    /* Merge the intercomm, to reconstruct an intracomm (we check
     * that this operation worked before we proceed further) */
    
    printf("%d: Trying to merge intercomm\n", srank);
    rc = MPI_Intercomm_merge(icomm, 1, &mcomm);
    rflag = flag = (MPI_SUCCESS==rc);
    MPIX_Comm_agree(scomm, &flag);
    if( MPI_COMM_WORLD != scomm ) MPI_Comm_free(&scomm);
    MPIX_Comm_agree(icomm, &rflag);
    MPI_Comm_free(&icomm);
    if( !(flag && rflag) ) {
        if( MPI_SUCCESS == rc ) {
            MPI_Comm_free(&mcomm);
        }
        fprintf(stderr, "%d: Intercomm_merge failed, redo\n", rank);
        goto redo;
    }

    /* Now, reorder mcomm according to original rank ordering in comm
     * Split does the magic: removing spare processes and reordering ranks
     * so that all surviving processes remain at their former place */
    printf("%d: Splitting comm\n", crank);
    rc = MPI_Comm_split(mcomm, 1, crank, newcomm);

    /* Split or some of the communications above may have failed if
     * new failures have disrupted the process: we need to
     * make sure we succeeded at all ranks, or retry until it works. */
    printf("%d: Ensure operation succeeds\n", crank);
    flag = (MPI_SUCCESS==rc);
    MPIX_Comm_agree(mcomm, &flag);
    printf("%d: Operation succeeded\n", crank);
    MPI_Comm_free(&mcomm);
    printf("%d: Freed mcomm\n", crank);
    if( !flag ) {
        if( MPI_SUCCESS == rc ) {
            printf("%d: freeing new comm\n", crank);
            MPI_Comm_free( newcomm );
            printf("%d: freed newcomm\n", crank);
        }
        fprintf(stderr, "%04d: comm_split failed, redo\n", rank);
        goto redo;
    }

    /* restore the error handler */
    if( MPI_COMM_NULL != comm ) {
        printf("%d: restoring error handler\n", crank);
        MPI_Errhandler errh;
        MPI_Comm_get_errhandler( comm, &errh );
        MPI_Comm_set_errhandler( *newcomm, errh );
    }
    printf("%d: Successfully managed error\n",crank);
    return MPI_SUCCESS;
}

static void MPIX_Error_handler(MPI_Comm *pcomm, int *perr, ...)
{
    MPIX_Comm_replace( world, &rworld );
    printf("%d: Replaced comm\n", rank);
    MPI_Comm_free( &world );
    world = rworld;
    // set stage for everyone
    MPI_Bcast(&stage, 1, MPI_INT, 0, world);
}

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
    if (rank == MASTER){
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
    
kern_start:
    printf("%d: Stage selected: %d\n",rank,stage);
    int err = 0;
    // restore stage if failed
    
    
    switch (stage) {
        case 2:
            if (rank == MASTER) goto axb; else goto axbw;
            
        case 3:
            if (rank == MASTER) goto cxd; else goto cxdw;
            
        case 4:
            if (rank == MASTER) goto exf; else goto exfw;
            
        default:
            break;
    }
    
    
    if (rank == MASTER)
    {
#ifndef BENCH
        printf("MPI_kernel_3mm has started with %d tasks.\n",np);
#endif
        
        /* Measure start time */
        MPI_timer_start();
    axb:
        stage = 2;
        MPI_Bcast(&stage, 1, MPI_INT, 0, world);
        /* Send A, B matrix data to the worker tasks */
        MPI_Comm_size(world, &np);
        numworkers = np - 1;
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
            rc = MPI_Isend(&offset, 1, MPI_INT, dest, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error sending data to %d, retrying stage\n", rank, dest);
                goto kern_start;
            }
            printf("Offset sent to %d\n",dest);
#endif
            rc = MPI_Isend(&rows, 1, MPI_INT, dest, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error sending data to %d, retrying stage\n", rank, dest);
                goto kern_start;
            }
            printf("Rows sent to %d\n",dest);
#endif
            rc = MPI_Isend(&A[offset][0], rows * NK, MPI_DOUBLE, dest, mtype,
                     world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error sending data to %d, retrying stage\n", rank, dest);
                goto kern_start;
            }
            printf("A sent to %d\n",dest);
#endif
            rc = MPI_Isend(&B, NK * NJ, MPI_DOUBLE, dest, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error sending data to %d, retrying stage\n", rank, dest);
                goto kern_start;
            }
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
            printf("%d: Requesting data from %d task\n", rank, i);
#endif
            rc = MPI_Irecv(&offset, 1, MPI_INT, source, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error getting data from %d, retrying stage\n", rank, source);
                goto kern_start;
            }
            printf("%d: Received offset from task %d for matrix E\n", rank, source);
#endif
            rc = MPI_Irecv(&rows, 1, MPI_INT, source, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error getting data from %d, retrying stage\n", rank, source);
                goto kern_start;
            }
            printf("%d: Received rows from task %d for matrix E\n", rank, source);
#endif
            rc = MPI_Irecv(&E[offset][0], rows * NJ, MPI_DOUBLE, source, mtype,
                     world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error getting data from %d, retrying stage\n", rank, source);
                goto kern_start;
            }
            printf("Received results from task %d for matrix E, waiting for %d tasks\n",source, numworkers-i);
#endif
        }
        
#ifndef BENCH
        printf ("Awaiting sync...\n");
#endif
        
        // STAGE 3
    cxd:
        stage = 3;
        MPI_Bcast(&stage, 1, MPI_INT, 0, world);
        MPI_Comm_size(world, &np);
        numworkers = np - 1;
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
            rc = MPI_Isend(&offset, 1, MPI_INT, dest, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error sending data to %d, retrying stage\n", rank, dest);
                goto kern_start;
            }
            printf("Offset sent to %d\n", dest);
#endif
            rc = MPI_Isend(&rows, 1, MPI_INT, dest, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error sending data to %d, retrying stage\n", rank, dest);
                goto kern_start;
            }
            printf("Rows sent to %d\n",dest);
#endif
            rc = MPI_Isend(&C[offset][0], rows * NM, MPI_DOUBLE, dest, mtype,
                     world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error sending data to %d, retrying stage\n", rank, dest);
                goto kern_start;
            }
            printf("C sent to %d\n",dest);
#endif
            rc = MPI_Isend(&D, NM * NL, MPI_DOUBLE, dest, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error sending data to %d, retrying stage\n", rank, dest);
                goto kern_start;
            }
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
            printf("%d: Requesting data from %d task\n", rank, i);
#endif
            rc = MPI_Irecv(&offset, 1, MPI_INT, source, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error getting data from %d, retrying stage\n", rank, source);
                goto kern_start;
            }
            printf("%d: Received offset from task %d for matrix E\n", rank, source);
#endif
            rc = MPI_Irecv(&rows, 1, MPI_INT, source, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error getting data from %d, retrying stage\n", rank, source);
                goto kern_start;
            }
            printf("%d: Received rows from task %d for matrix E\n", rank, source);
#endif
            rc = MPI_Irecv(&F[offset][0], rows * NL, MPI_DOUBLE, source, mtype,
                     world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error getting data from %d, retrying stage\n", rank, source);
                goto kern_start;
            }
            printf("Received results from task %d for matrix F\n",source);
#endif
        }
        
#ifndef BENCH
        printf ("Awaiting sync...\n");
#endif
        
        
        // STAGE 4
    exf:
        stage = 4;
        MPI_Comm_size(world, &np);
        numworkers = np - 1;
        MPI_Bcast(&stage, 1, MPI_INT, 0, world);
        
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
            rc = MPI_Isend(&offset, 1, MPI_INT, dest, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error sending data to %d, retrying stage\n", rank, dest);
                goto kern_start;
            }
            printf("Offset sent to %d\n", dest);
#endif
            rc = MPI_Isend(&rows, 1, MPI_INT, dest, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error sending data to %d, retrying stage\n", rank, dest);
                goto kern_start;
            }
            printf("Rows sent to %d\n",dest);
#endif
            rc = MPI_Isend(&E[offset][0], rows * NJ, MPI_DOUBLE, dest, mtype,
                     world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error sending data to %d, retrying stage\n", rank, dest);
                goto kern_start;
            }
            printf("E sent to %d\n",dest);
#endif
            rc = MPI_Isend(&F, NJ * NL, MPI_DOUBLE, dest, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error sending data to %dest, retrying stage\n", rank, dest);
                goto kern_start;
            }
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
            printf("%d: Requesting data from %d task\n", rank, i);
#endif
            rc = MPI_Irecv(&offset, 1, MPI_INT, source, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error getting data from %d, retrying stage\n", rank, source);
                goto kern_start;
            }
            printf("%d: Received offset from task %d for matrix E\n", rank, source);
#endif
            rc = MPI_Irecv(&rows, 1, MPI_INT, source, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error getting data from %d, retrying stage\n", rank, source);
                goto kern_start;
            }
            printf("%d: Received rows from task %d for matrix E\n", rank, source);
#endif
            rc = MPI_Irecv(&G[offset][0], rows * NL, MPI_DOUBLE, source, mtype,
                     world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error getting data from %d, retrying stage\n", rank, source);
                goto kern_start;
            }
            printf("Received results from task %d\n for matrix G\n",source);
#endif
        }
        
#ifndef BENCH
        printf ("Awaiting sync...\n");
#endif
        
        
#ifndef BENCH
        printf("%d: All workers exited\n", rank);
#endif
        
        MPI_timer_end();
        //MPI_timer_print();
        
    }
    
    // Workers
    
    if (rank > MASTER)
    {
            // Workers calculating E
        
        // STAGE 2
    axbw:
        stage = 2;
    
        {
#ifndef BENCH
            printf("E calculation in %d task\n", rank);
#endif
            
            mtype = FROM_MASTER;
#ifndef BENCH
            printf("%d: Requesting data from MASTER\n", rank);
#endif
            rc = MPI_Irecv(&offset, 1, MPI_INT, MASTER, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error getting data from MASTER, retrying stage\n", rank);
                goto kern_start;
            }
            printf("%d: Got offset from MASTER\n", rank);
#endif
            rc = MPI_Irecv(&rows, 1, MPI_INT, MASTER, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error getting data from MASTER, retrying stage\n", rank);
                goto kern_start;
            }
            printf("%d: Got rows from MASTER\n", rank);
#endif
            rc = MPI_Irecv(&A, rows * NK, MPI_DOUBLE, MASTER, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error getting data from MASTER, retrying stage\n", rank);
                goto kern_start;
            }
            printf("%d: Got first matrix from MASTER\n", rank);
#endif
            rc = MPI_Irecv(&B, NK * NJ, MPI_DOUBLE, MASTER, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error getting data from MASTER, retrying stage\n", rank);
                goto kern_start;
            }
            printf("%d: Got second matrix from MASTER\n", rank);
            printf("Received all data for task %d\n", rank);
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
            printf("Result: sending data from %d task to MASTER\n", rank);
#endif
            rc = MPI_Isend(&offset, 1, MPI_INT, MASTER, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error sending data to MASTER, retrying stage\n", rank);
                goto kern_start;
            }
            printf("Result: send offset from %d task to MASTER\n", rank);
#endif
            rc = MPI_Isend(&rows, 1, MPI_INT, MASTER, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error sending data to MASTER, retrying stage\n", rank);
                goto kern_start;
            }
            printf("Result: sent rows from %d task to MASTER\n", rank);
#endif
            rc = MPI_Isend(&E, rows * NJ, MPI_DOUBLE, MASTER, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error sending data to MASTER, retrying stage\n", rank);
                goto kern_start;
            }
            printf("Successfully sent all data from %d task to MASTER\n", rank);
            printf("Awaiting sync...\n");
#endif
        }
        
        
            // Workers calculating F
        
        
        // STAGE 3
    cxdw:
        stage = 3;
        
        {
            mtype = FROM_MASTER;
#ifndef BENCH
            printf("%d: Requesting data from MASTER\n", rank);
#endif
            rc = MPI_Irecv(&offset, 1, MPI_INT, MASTER, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error getting data from MASTER, retrying stage\n", rank);
                goto kern_start;
            }
            printf("%d: Got offset from MASTER\n", rank);
#endif
            rc = MPI_Irecv(&rows, 1, MPI_INT, MASTER, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error getting data from MASTER, retrying stage\n", rank);
                goto cxdw;
            }
            printf("%d: Got rows from MASTER\n", rank);
#endif
            rc = MPI_Irecv(&C, rows * NM, MPI_DOUBLE, MASTER, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error getting data from MASTER, retrying stage\n", rank);
                goto kern_start;
            }
            printf("%d: Got first matrix from MASTER\n", rank);
#endif
            rc = MPI_Irecv(&D, NM * NL, MPI_DOUBLE, MASTER, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error getting data from MASTER, retrying stage\n", rank);
                goto kern_start;
            }
            printf("%d: Got second matrix from MASTER\n", rank);
            printf("Received all data for task %d\n", rank);
#endif
            
            for (int k = 0; k < NL; k++)
                for (int i = 0; i < rows; i++)
                {
                    F[i][k] = 0.0;
                    
                    // ERROR HERE
                    if (setdie && rank == 2 && k == 20){
                        printf("Rank %d: Goodbye, cruel world!\n", rank);
                        raise( SIGKILL );
                    }
                    
                    for (int j = 0; j < NJ; j++)
                        F[i][k] = F[i][k] + C[i][j] * D[j][k];
                }
            mtype = FROM_WORKER;
#ifndef BENCH
            printf("Result: sending data from %d task to MASTER\n", rank);
#endif
            rc = MPI_Isend(&offset, 1, MPI_INT, MASTER, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error sending data to MASTER, retrying stage\n", rank);
                goto kern_start;
            }
            printf("Result: sent offset from %d task to MASTER\n", rank);
#endif
            rc = MPI_Isend(&rows, 1, MPI_INT, MASTER, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error sending data to MASTER, retrying stage\n", rank);
                goto kern_start;
            }
            printf("Result: sent rows from %d task to MASTER\n", rank);
#endif
            rc = MPI_Isend(&F, rows * NL, MPI_DOUBLE, MASTER, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error sending data to MASTER, retrying stage\n", rank);
                goto kern_start;
            }
            printf("Successfully sent all data from %d task to MASTER\n", rank);
            printf("Awaiting sync...\n");
#endif
        }
        
            // Workers calculating G
        // STAGE 4
    exfw:
        stage = 4;

        {
            mtype = FROM_MASTER;
#ifndef BENCH
            printf("%d: Requesting data from MASTER\n", rank);
#endif
            rc = MPI_Irecv(&offset, 1, MPI_INT, MASTER, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error getting data from MASTER, retrying stage\n", rank);
                goto kern_start;
            }
            printf("%d: Got offset from MASTER\n", rank);
#endif
            rc = MPI_Irecv(&rows, 1, MPI_INT, MASTER, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error getting data from MASTER, retrying stage\n", rank);
                goto kern_start;
            }
            printf("%d: Got rows from MASTER\n", rank);
#endif
            rc = MPI_Irecv(&E, rows * NJ, MPI_DOUBLE, MASTER, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error getting data from MASTER, retrying stage\n", rank);
                goto kern_start;
            }
            printf("%d: Got first matrix from MASTER\n", rank);
#endif
            rc = MPI_Irecv(&F, NJ * NL, MPI_DOUBLE, MASTER, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error getting data from MASTER, retrying stage\n", rank);
                goto kern_start;
            }
            printf("%d: Got second matrix from MASTER\n", rank);
            printf("Received all data for task %d\n", rank);
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
            printf("Result: sending data from %d task to MASTER\n", rank);
#endif
            rc = MPI_Isend(&offset, 1, MPI_INT, MASTER, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error sending data to MASTER, retrying stage\n", rank);
                goto kern_start;
            }
            printf("Result: sent offset from %d task to MASTER\n", rank);
#endif
            rc = MPI_Isend(&rows, 1, MPI_INT, MASTER, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
#ifndef BENCH
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error sending data to MASTER, retrying stage\n", rank);
                goto kern_start;
            }
            printf("Result: sent rows from %d task to MASTER\n", rank);
#endif
            rc = MPI_Isend(&G, rows * NL, MPI_DOUBLE, MASTER, mtype, world, &request);
            err = MPIX_Wait_timeout(&request, 3000);
            if (err != MPI_SUCCESS){
                printf("rank %d: There was an error sending data to MASTER, retrying stage\n", rank);
                goto kern_start;
            }
        }
#ifndef BENCH
        printf("Successfully sent all data from %d task to MASTER\n", rank);
        printf("Awaiting sync...\n");
#endif
        
        MPI_Barrier(world);
        
#ifndef BENCH
        printf("%d: Done working.\n", rank);
#endif
        
    }
}

int main( int argc, char* argv[] )
{
    gargv = argv;
  
    MPI_Init( &argc, &argv );
    MPI_Errhandler errh;
    
    MPI_Comm_get_parent( &world );
    if( MPI_COMM_NULL == world ) {
        /* First run: Let's create an initial world,
         * a copy of MPI_COMM_WORLD */
        MPI_Comm_dup( MPI_COMM_WORLD, &world );
        MPI_Comm_size( world, &np );
        MPI_Comm_rank( world, &rank );
        
        MPI_Comm_create_errhandler(MPIX_Error_handler, &errh);
        MPI_Comm_set_errhandler(world, errh);
    } else {
        /* I am a spare, lets get the repaired world */
        printf("Spare process created!\n");
        MPIX_Comm_replace( MPI_COMM_NULL, &world );
        MPI_Comm_size( world, &np );
        MPI_Comm_rank( world, &rank );
        printf("Spare rank: %d, seen number of processes: %d\n", rank,np);
        MPI_Comm_create_errhandler(MPIX_Error_handler, &errh);
        MPI_Comm_set_errhandler(world, errh);
        setdie = false;
        goto joinwork;
    }
  
    if (rank == MASTER){
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
    }
    
    if (np < 2 ) {
        printf("Guess I'll die ¯\\_(ツ)_/¯ \n");
        raise ( SIGKILL );
        // Basically, there's no point doing this;
        // This is left here as legacy for future references
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
    
    numworkers = np-1;
    if (rank == MASTER && !spare){
        MPI_init_matrixes();
    }
    
    
    /* The victim is now the median process (for simplicity) */
    printf("Choosing victim from %d processes\n", np);
    victim = (rank == np / 2)? 1 : 0;
    
    if (rank == MASTER){
#ifndef BENCH
        printf("Tasks comm port set successful\n");
#endif
        bench_timer_start();
    }
    
joinwork:
    
    MPI_kernel_3mm();
    
    
    if (rank == MASTER){
        bench_timer_stop();
        bench_timer_print();
    }
        
    
    if ( (rank == MASTER) && argc > 42 && !strcmp( argv[0], "" ) )
        print_array( ni, nl, G );
    
#ifndef BENCH
    if (rank == MASTER)
        printf( "Just finished!\n" );
#endif
    MPI_Finalize();
    return 0;
}

