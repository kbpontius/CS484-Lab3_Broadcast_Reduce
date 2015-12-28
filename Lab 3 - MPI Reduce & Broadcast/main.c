//
//  main.c
//  Lab 3 - MPI Reduce & Broadcast
//
//  Created by Kyle on 12/15/15.
//  Copyright © 2015 Kyle Pontius. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sys/time.h>
#include "mpi.h"

#define VECSIZE 100000
#define ITERATIONS 100

// I decided to just use an int array, then reduce that. It saved extra time
// in figuring out how to send a struct through MPI.
typedef struct {
    double val;
    int rank;
} element;

void printArray(double vector[VECSIZE]) {
    int i;
    for (i = 0; i < VECSIZE; i++) {
        fprintf(stderr, "\t----> received vector[%i] = %f\n", i, vector[i]);
    }
}

double When()
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return ((double) tp.tv_sec + (double) tp.tv_usec * 1e-6);
}

void getMaxVectorValues(double myVector[VECSIZE], double receivedVector[VECSIZE]) {
    int i;
    
    for (i = 0; i < VECSIZE; i++) {
        if (myVector[i] < receivedVector[i]) {
//            printf("FOUND LARGER VALUE: %f > %f\n", receivedVector[i], myVector[i]);
            myVector[i] = receivedVector[i];
        }
    }
}

void my_MPI_Reduce_Max(double *vector, int size, int numDim, int rank, MPI_Status status, char hostName[]) {
    int notParticipating = 0;
    int bitmask = 1;
    int curDim;
    
    for (curDim = 0; curDim < numDim; curDim++) {
        if ((rank & notParticipating) == 0) {
            if ((rank & bitmask) != 0) {
                int msgDestination = rank ^ bitmask;
                MPI_Send(vector, VECSIZE, MPI_DOUBLE, msgDestination, 0, MPI_COMM_WORLD);
//                printf("\nHOSTNAME SENT: --------%s-------- || RANK: %i || DESTINATION: %i\n", hostName, rank, msgDestination);
            } else {
                int msgSource = rank ^ bitmask;
                double receivedArray[VECSIZE];
                MPI_Recv(receivedArray, VECSIZE, MPI_DOUBLE, msgSource, 0, MPI_COMM_WORLD, &status);
//                printf("\nHOSTNAME RECEIVED: --------%s-------- || RANK: %i || SOURCE: %i\n", hostName, rank, msgSource);
                
                getMaxVectorValues(vector, receivedArray);
            }
        }
        
        notParticipating = notParticipating ^ bitmask;
        bitmask <<= 1;
    }
}

void my_MPI_Broadcast_Max(double vector[VECSIZE], int size, int numDim, int rank, MPI_Status status) {
    int notParticipating = pow(2, numDim - 1) - 1;
    int bitmask = pow(2, numDim - 1);
    int curDim;
    double receivedArray[VECSIZE];
    
    for(curDim = 0; curDim < numDim; curDim++) {
        if ((rank & notParticipating) == 0) {
            if ((rank & bitmask) == 0) {
                int msgDestination = rank ^ bitmask;
                MPI_Send(vector, VECSIZE, MPI_DOUBLE, msgDestination, 0, MPI_COMM_WORLD);
            } else {
                int msgSource = rank ^ bitmask;
                MPI_Recv(receivedArray, VECSIZE, MPI_DOUBLE, msgSource, 0, MPI_COMM_WORLD, &status);
                
                int i;
                for (i = 0; i < VECSIZE; i++) {
                    vector[i] = receivedArray[i];
                }
            }
        }
        
        notParticipating >>= 1;
        bitmask >>=1;
    }
}

int main(int argc, char *argv[])
{
    int nproc, j, iter;
    int myRank, root = 0;
    char hostName[255];
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Status status;
    
    int numDim = (int)(log2(nproc));
    double vector[VECSIZE];
    
    gethostname(hostName, 253);
    
    /* EXECUTE VECTOR MAX BCAST & REDUCE */
    
    srand(myRank + 5);
    double startMaxVector = When();
    
    for (j = 0; j < VECSIZE; j++) {
        vector[j] = rand();
//        fprintf(stderr, "%i: vector[%i] = %f\n", myRank, j, vector[j]);
    }
    
    // Reduce the vector the root node using the hypercube.
    my_MPI_Reduce_Max(vector, VECSIZE, numDim, myRank, status, hostName);
    
    // Uncomment for error checking.
//    if (myRank == root) {
//        for (iter = 0; iter < VECSIZE; iter++) {
//            fprintf(stderr, "root vector[%d] = %f\n", iter, vector[iter]);
//        }
//    }
    
    // Now broadcast this max vector to everyone else.
    my_MPI_Broadcast_Max(vector, VECSIZE, numDim, myRank, status);
    
    // Uncomment for error checking.
//    for (iter = 0; iter < VECSIZE; iter++) {
//        fprintf(stderr, "final proc %d [%d] = %f\n", myRank, iter, vector[iter]);
//    }
    
    double endMaxVector = When();
    
    /* MAXIMUM BANDWIDTH CALCULATION */
    
    double startBandwidth = 0.0, endBandwidth = 0.0;
    
    if (myRank == root || myRank == 1) {
        
        startBandwidth = When();
        
        for (iter = 0; iter < ITERATIONS; iter++) {
            if (myRank == root) {
                MPI_Send(vector, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
            } else if (myRank == 1) {
                MPI_Recv(vector, 1, MPI_DOUBLE, root, 0, MPI_COMM_WORLD, &status);
            }
        }
        
        endBandwidth = When();
    }
    
    MPI_Finalize();
    
    
    if (myRank == root) {
        double bandwidthTime = endBandwidth - startBandwidth;
        double bytesPerSecond = sizeof(double) * ITERATIONS / bandwidthTime;
        fprintf(stderr, "Time %f\n", endMaxVector - startMaxVector);
        fprintf(stderr, "MAX BANDWIDTH: %f bytes/second\n", bytesPerSecond);
    }
    
    return 0;
}
