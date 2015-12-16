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

#define VECSIZE 8
#define ITERATIONS 10000

typedef struct {
    double val;
    int rank;
} element;

double When()
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return ((double) tp.tv_sec + (double) tp.tv_usec * 1e-6);
}
//
//// Reduce values to one node
//float ReduceSum(int numdim, int rank, float value)
//{
//    int notparticipating = 0;
//    int bitmask = 1;
//    float sum = value;
//    float newvalue;
//    for(int i = 0; i < numdim; i++) {
//        if ((rank & notparticipating) == 0) {
//            if ((rank & bitmask) != 0) {
//                int msg_dest = rank ^ bitmask;
//                send(sum, msg_dest);
//            } else {
//                int msg_src = rank ^ bitmask;
//                recv(&newvalue, msg_src);
//                sum += newvalue;
//            }
//        }
//        notparticipating = notparticipating ^ bitmask;
//        bitmask <<=1;
//    }
//    return(sum);
//}

void my_MPI_Reduce_Max(element vector[VECSIZE], int size, int numDim, int myRank) {
    int notParticipating = 0;
    int bitmask = 1;
    float sum = // value;
    float newValue;
    
    for (int curDim = 0; curDim < numDim; curDim++) {
        if ((myRank & notParticipating) == 0) {
            if ((rank & bitmask) != 0) {
                int msgDestination = rank ^ bitmask;
                MPI_Send(<#const void *buf#>, <#int count#>, <#MPI_Datatype datatype#>, <#int dest#>, <#int tag#>, <#MPI_Comm comm#>)
            } else {
                int msgSource = rank ^ bitmask;
                MPI_Recv(<#void *buf#>, <#int count#>, <#MPI_Datatype datatype#>, <#int source#>, <#int tag#>, <#MPI_Comm comm#>, <#MPI_Status *status#>)
            }
        }
    }
}

void my_MPI_Broadcast_Max(element vector[VECSIZE], int size, int numDim, int myRank) {
    
}

int main(int argc, char *argv[])
{
    int nproc, i, iter;
    int myRank, root = 0;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    
    int numDim = (int)(log2(nproc));
    element vector[VECSIZE];
    
    // Start time here
    srand(myRank + 5);
    double start = When();
    
    for (iter = 0; iter < ITERATIONS; iter++) {
        for (i = 0; i < VECSIZE; i++) {
            vector[i].val = rand();
            vector[i].rank = myRank;
//             printf("init proc %d [%d]=%f\n", myRank, i, ain[i]);
        }
        
        my_MPI_Reduce_Max(vector, VECSIZE, numDim, myRank);
        
        // At this point, the answer resides on process root
        if (myRank == root) {
            for (i = 0; i < VECSIZE; i++) {
                printf("------------ root vector[%d] = %f from %d\n", i, vector[i].val, vector[i].rank);
            }
        }
        
        // Now broadcast this max vector to everyone else.
        my_MPI_Broadcast_Max(vector, VECSIZE, numDim, myRank);
        
         for (i = 0; i < VECSIZE; i++) {
             printf("final proc %d [%d]=%f from %d\n", myRank, i, vector[i].val, vector[i].rank);
         }
    }
    
    MPI_Finalize();
    
    double end = When();
    if (myRank == root) {
        printf("Time %f\n", end - start);
    }
    
    return 0;
}
