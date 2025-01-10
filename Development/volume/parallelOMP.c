#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <stddef.h>

#define ROOT_SELF (-1)
#define MAX_ITER 100

typedef struct {
    int src, dst, weight;
} Edge;

/* 
File format should be:
    line 1: <num_vertices> <num_edges>
    next lines: <src> <dst> <weight>
*/
void readGraph(const char *filename, Edge **edges, int *numVertices, int *numEdges, int world_rank) {
    FILE *f = fopen(filename, "r");
    if (!f) {
        fprintf(stderr, "Error: Cannot open file %s\n", filename);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    fscanf(f, "%d %d", numVertices, numEdges);

    // Reads and store Edges on the defined structure
    *edges = (Edge *)malloc((*numEdges) * sizeof(Edge));
    int i;
    for (i = 0; i < *numEdges; i++) {
        fscanf(f, "%d %d %d", &(*edges)[i].src, &(*edges)[i].dst, &(*edges)[i].weight);
    }
    fclose(f);
}

void distributeEdges(Edge *globalEdges, int totalEdges, Edge **localEdges, int *localCount, int world_rank, int world_size, MPI_Datatype MPI_EDGE_TYPE) {
    // Equal distribution of edges
    int baseCount = totalEdges / world_size;
    int remainder = totalEdges % world_size;

    int *counts = (int *)malloc(world_size * sizeof(int));
    int *displs = (int *)malloc(world_size * sizeof(int));

    int i;
    for (i = 0; i < world_size; i++) {
        // First ranks get 1 extra edge
        counts[i] = baseCount + (i < remainder ? 1 : 0);
        displs[i] = (i == 0) ? 0 : displs[i - 1] + counts[i - 1];
    }

    *localCount = counts[world_rank];
    *localEdges = (Edge *)malloc((*localCount) * sizeof(Edge));

    MPI_Scatterv(globalEdges, counts, displs, MPI_EDGE_TYPE, *localEdges, *localCount, MPI_EDGE_TYPE, 0, MPI_COMM_WORLD);

    free(counts);
    free(displs);
}

int findRoot(int v, int *parent) {
    if (parent[v] == ROOT_SELF || parent[v] == v)
        return v;

    // Path compression
    parent[v] = findRoot(parent[v], parent);
    return parent[v];
}

void unionRoots(int r1, int r2, int *parent) {
    r1 = findRoot(r1, parent);
    r2 = findRoot(r2, parent);
    if (r1 != r2)
        parent[r2] = r1;
}

void minEdgeReduce(void *in, void *inout, int *len, MPI_Datatype *datatype) {
    Edge *input = (Edge *)in;
    Edge *inoutResult = (Edge *)inout;

    // Reduction
    int i;
    for (i = 0; i < *len; i++)
        if (input[i].weight < inoutResult[i].weight)
            inoutResult[i] = input[i];
}

void boruvkaMST(int numVertices, Edge *localEdges, int localCount, int *parent, int world_rank, int world_size, Edge **mstEdges, int *mstCount, MPI_Datatype MPI_EDGE_TYPE) {
    int *bestWeight = (int *)malloc(numVertices * sizeof(int));
    Edge *bestEdge = (Edge *)malloc(numVertices * sizeof(Edge));

    // Initialize parent array
    int v;
    for (v = 0; v < numVertices; v++)
        parent[v] = ROOT_SELF;

    // Initialize MST edges
    if (world_rank == 0)
        *mstEdges = (Edge *)malloc((numVertices - 1) * sizeof(Edge));

    int done = 0;

    // Create custom MPI_Op for Edge reduction
    MPI_Op MPI_MIN_EDGE_OP;
    MPI_Op_create((MPI_User_function *)minEdgeReduce, 1, &MPI_MIN_EDGE_OP);

    while (!done) {
        // Initialize best edges and weights
        for (v = 0; v < numVertices; v++) {
            bestWeight[v] = INT_MAX;
            bestEdge[v].src = -1;
            bestEdge[v].dst = -1;
            bestEdge[v].weight = INT_MAX;
        }

        // Find best edge for each component
        int i;
        #pragma omp parallel for private(i)
        for (i = 0; i < localCount; i++) {
            int src = localEdges[i].src;
            int dst = localEdges[i].dst;
            int weight = localEdges[i].weight;

            int r1 = findRoot(src, parent);
            int r2 = findRoot(dst, parent);

            if (r1 != r2) {
                #pragma omp critical
                {
                    if (weight < bestWeight[r1]) {
                        bestWeight[r1] = weight;
                        bestEdge[r1] = localEdges[i];
                    }
                    if (weight < bestWeight[r2]) {
                        bestWeight[r2] = weight;
                        bestEdge[r2] = localEdges[i];
                    }
                }
            }
        }

        // Find global best edges using custom MPI_Op
        Edge *globalBestEdge = (Edge *)malloc(numVertices * sizeof(Edge));
        MPI_Allreduce(bestEdge, globalBestEdge, numVertices, MPI_EDGE_TYPE, MPI_MIN_EDGE_OP, MPI_COMM_WORLD);

        // Add global best edges to MST
        if (world_rank == 0) {
            for (v = 0; v < numVertices; v++) {
                if (globalBestEdge[v].weight != INT_MAX) {
                    int src = globalBestEdge[v].src;
                    int dst = globalBestEdge[v].dst;

                    int rootSrc = findRoot(src, parent);
                    int rootDst = findRoot(dst, parent);

                    if (rootSrc != rootDst) {
                        unionRoots(rootSrc, rootDst, parent);
                        (*mstEdges)[(*mstCount)++] = globalBestEdge[v];
                    }
                }
            }
        }
        
        // Count the number of components
        int rootCount = 0;
        if (world_rank == 0)
            for (v = 0; v < numVertices; v++) 
                if (parent[v] == ROOT_SELF && findRoot(v, parent) == v) 
                    rootCount++;

        // Broadcast the globalRootCount from rank 0 to all processes
        MPI_Bcast(&rootCount, 1, MPI_INT, 0, MPI_COMM_WORLD);


        // Each process checks the termination condition
        if (rootCount <= 1) {
            done = 1;
        } else {
            // Broadcast the updated parent array to all ranks if not terminated
            MPI_Bcast(parent, numVertices, MPI_INT, 0, MPI_COMM_WORLD);
        }
        free(globalBestEdge);
    }

    // Cleanup
    MPI_Op_free(&MPI_MIN_EDGE_OP);
    free(bestWeight);
    free(bestEdge);
}

void writeMST(Edge *mstEdges, int mstCount) {
    // Calculate total cost and write MST edges and total cost to file
    long long totalCost = 0;
    int i;
    for (i = 0; i < mstCount; i++)
        totalCost += mstEdges[i].weight;

    printf("Total cost = %lld\n", totalCost);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    // Synchronize all processes before starting the timer
    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Number of OpenMP threads per MPI process
    omp_set_num_threads(atoi(argv[1]));

    int numVertices = 0, numEdges = 0;
    Edge *globalEdges = NULL;

    // Define MPI_Datatype for Edge
    MPI_Datatype MPI_EDGE_TYPE;
    {
        int blocklengths[3] = {1, 1, 1};
        MPI_Aint offsets[3] = {offsetof(Edge, src), offsetof(Edge, dst), offsetof(Edge, weight)};
        MPI_Datatype types[3] = {MPI_INT, MPI_INT, MPI_INT};

        MPI_Type_create_struct(3, blocklengths, offsets, types, &MPI_EDGE_TYPE);
        MPI_Type_commit(&MPI_EDGE_TYPE);
    }

    // Read graph on rank 0
    if (world_rank == 0)
        readGraph("graph.txt", &globalEdges, &numVertices, &numEdges, world_rank);

    // Broadcast number of vertices and edges
    MPI_Bcast(&numVertices, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numEdges, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Distribute edges
    Edge *localEdges = NULL;
    int localCount = 0;
    distributeEdges(globalEdges, numEdges, &localEdges, &localCount, world_rank, world_size, MPI_EDGE_TYPE);

    // Rank 0 can free the global edges
    if (world_rank == 0)
        free(globalEdges);

    // Allocate parent array
    int *parent = (int *)malloc(numVertices * sizeof(int));

    // Initialize MST edges list
    Edge *mstEdges = NULL;
    int mstCount = 0;

    // Borůvka MST
    boruvkaMST(numVertices, localEdges, localCount, parent, world_rank, world_size, &mstEdges, &mstCount, MPI_EDGE_TYPE);

    // Stop the timer after MST computation completes
    double end_time = MPI_Wtime();
    double elapsed_time = end_time - start_time;

    // Gather the maximum elapsed_time across all processes
    double max_time;
    MPI_Reduce(&elapsed_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // Only rank 0 writes the MST and prints the timing
    if (world_rank == 0) {
        writeMST(mstEdges, mstCount);
        printf("MST Computation Time (excluding writing): %lf seconds\n", max_time);
    }

    // Free allocated memory
    free(localEdges);
    free(parent);
    if (mstEdges)
        free(mstEdges);

    // Free custom MPI datatype
    MPI_Type_free(&MPI_EDGE_TYPE);

    MPI_Finalize();
    return 0;
}