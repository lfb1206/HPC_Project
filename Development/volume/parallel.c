#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <stddef.h>

typedef struct {
    int src, dst, weight;
} Edge;

/*
File format:
    line 1: <num_vertices> <num_edges>
    next lines: <src> <dst> <weight>
*/
void readGraph(const char *filename, Edge **edges, int *numVertices, int *numEdges) {
    FILE *f = fopen(filename, "r");
    if (!f) {
        fprintf(stderr, "Error: Cannot open file %s\n", filename);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    fscanf(f, "%d %d", numVertices, numEdges);

    *edges = (Edge *)malloc((*numEdges) * sizeof(Edge));

    int i;
    for (i = 0; i < *numEdges; i++) {
        if (fscanf(f, "%d %d %d", &(*edges)[i].src, &(*edges)[i].dst, &(*edges)[i].weight) != 3) {
            fprintf(stderr, "Error: Invalid edge data at line %d\n", i + 2);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    fclose(f);
}

void distributeEdges(Edge *globalEdges, int totalEdges, Edge **localEdges, int *localCount, int world_rank, int world_size, MPI_Datatype MPI_EDGE_TYPE) {
    // Equal distribution of edges
    int baseCount = totalEdges / world_size;
    int remainder = totalEdges % world_size;

    int *counts = (int *)malloc(world_size * sizeof(int));
    int *displs = (int *)malloc(world_size * sizeof(int));

    if (!counts || !displs) {
        fprintf(stderr, "Error: Not enough memory for distribution arrays\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Dislplacement for each array in globalEdges
    displs[0] = 0;
    int i;
    for (i = 0; i < world_size; i++) {
        counts[i] = baseCount + ((i < remainder) ? 1 : 0);
        if (i > 0)
            displs[i] = displs[i - 1] + counts[i - 1];
    }

    *localCount = counts[world_rank];
    *localEdges = (Edge *)malloc((*localCount) * sizeof(Edge));
    if (*localEdges == NULL && *localCount > 0) {
        fprintf(stderr, "Error: Not enough memory to store local edges\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MPI_Scatterv(globalEdges, counts, displs, MPI_EDGE_TYPE,
                 *localEdges, *localCount, MPI_EDGE_TYPE,
                 0, MPI_COMM_WORLD);

    free(counts);
    free(displs);
}

int findRoot(int v, int *parent, int *rank) {
    if (parent[v] != v) {
        parent[v] = findRoot(parent[v], parent, rank); // Path compression
    }
    return parent[v];
}

void unionRoots(int v1, int v2, int *parent, int *rank) {
    int r1 = findRoot(v1, parent, rank);
    int r2 = findRoot(v2, parent, rank);

    if (r1 != r2) {
        // Attach smaller rank tree under root of higher rank tree
        if (rank[r1] > rank[r2]) {
            parent[r2] = r1;
        } else if (rank[r1] < rank[r2]) {
            parent[r1] = r2;
        } else {
            // Same rank; choose one as root and increment its rank
            parent[r2] = r1;
            rank[r1]++;
        }
    }
}

void boruvkaMST(int numVertices, Edge *localEdges, int localCount, int *parent, int *rank, int world_rank, int world_size, long long *totalCost, MPI_Datatype MPI_EDGE_TYPE) {
    // Arrays to track best edge from each component
    int *bestWeight       = (int *)malloc(numVertices * sizeof(int));
    Edge *bestEdge        = (Edge *)malloc(numVertices * sizeof(Edge));
    Edge *globalBestEdge  = (Edge *)malloc(numVertices * sizeof(Edge));

    // For gathering best edges from all ranks (only rank 0 needs to allocate)
    Edge *gatherBuffer = NULL;
    if (world_rank == 0) {
        gatherBuffer = (Edge *)malloc(world_size * numVertices * sizeof(Edge));
        if (!gatherBuffer) {
            fprintf(stderr, "Error: Not enough memory for gatherBuffer\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    // Initialize parent and rank
    int v;
    for (v = 0; v < numVertices; v++) {
        parent[v] = v;
        rank[v]   = 0;
    }

    int done = 0;
    while (!done) {
        // Initialize "best" info for each vertex's component
        for (v = 0; v < numVertices; v++) {
            bestWeight[v]       = INT_MAX;
            bestEdge[v].src     = -1;
            bestEdge[v].dst     = -1;
            bestEdge[v].weight  = INT_MAX;
        }

        // Find best edge for each component
        int i;
        for (i = 0; i < localCount; i++) {
            int src  = localEdges[i].src;
            int dst  = localEdges[i].dst;
            int w    = localEdges[i].weight;

            int r1 = findRoot(src, parent, rank);
            int r2 = findRoot(dst, parent, rank);

            if (r1 != r2) {
                if (w < bestWeight[r1]) {
                    bestWeight[r1] = w;
                    bestEdge[r1]   = localEdges[i];
                }
                if (w < bestWeight[r2]) {
                    bestWeight[r2] = w;
                    bestEdge[r2]   = localEdges[i];
                }
            }
        }

        // Gather all bestEdges on rank 0
        MPI_Gather(bestEdge, numVertices, MPI_EDGE_TYPE,
                   gatherBuffer, numVertices, MPI_EDGE_TYPE,
                   0, MPI_COMM_WORLD);

        // Rank 0 combines the results to find the global best edge per component
        if (world_rank == 0) {
            // Initialize global best edges
            for (v = 0; v < numVertices; v++) {
                globalBestEdge[v].src    = -1;
                globalBestEdge[v].dst    = -1;
                globalBestEdge[v].weight = INT_MAX;
            }

            // Combine from all ranks
            int r;
            for (r = 0; r < world_size; r++) {
                for (v = 0; v < numVertices; v++) {
                    Edge candidate = gatherBuffer[r * numVertices + v];
                    if (candidate.weight < globalBestEdge[v].weight) {
                        globalBestEdge[v] = candidate;
                    }
                }
            }
        }

        // On rank 0, union the components for all global best edges and add to MST if it merges two components.
        if (world_rank == 0) {
            for (v = 0; v < numVertices; v++) {
                if (globalBestEdge[v].weight != INT_MAX) {
                    int src = globalBestEdge[v].src;
                    int dst = globalBestEdge[v].dst;
                    if (src == -1 || dst == -1) 
                        continue;  // no valid edge

                    int r1 = findRoot(src, parent, rank);
                    int r2 = findRoot(dst, parent, rank);
                    if (r1 != r2) {
                        unionRoots(r1, r2, parent, rank);
                        totalCost += globalBestEdge[v].weight;
                    }
                }
            }
        }

        // Rank 0 counts number of distinct roots
        int rootCount = 0;
        if (world_rank == 0) {
            for (v = 0; v < numVertices; v++) {
                if (parent[v] == v) {
                    rootCount++;
                }
            }
        }

        // Broadcast rootCount and updated parent, rank arrays
        MPI_Bcast(&rootCount, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (rootCount <= 1) {
            done = 1;
        } else {
            // Only broadcast if not finished
            MPI_Bcast(parent, numVertices, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(rank, numVertices, MPI_INT, 0, MPI_COMM_WORLD);
        }
    }

    free(bestWeight);
    free(bestEdge);
    free(globalBestEdge);
    if (gatherBuffer)
        free(gatherBuffer);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    // Synchronize all processes before timing
    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int numVertices = 0, numEdges = 0;
    Edge *globalEdges = NULL;

    // Define MPI_Datatype for Edge (weight as int)
    MPI_Datatype MPI_EDGE_TYPE;
    {
        int blocklengths[3] = {1, 1, 1};
        MPI_Aint offsets[3];
        offsets[0] = offsetof(Edge, src);
        offsets[1] = offsetof(Edge, dst);
        offsets[2] = offsetof(Edge, weight);

        MPI_Datatype types[3] = {MPI_INT, MPI_INT, MPI_INT};

        MPI_Type_create_struct(3, blocklengths, offsets, types, &MPI_EDGE_TYPE);
        MPI_Type_commit(&MPI_EDGE_TYPE);
    }

    // Read the graph on rank 0 only
    if (world_rank == 0) {
        readGraph("graph.txt", &globalEdges, &numVertices, &numEdges);
        if (numVertices <= 1) {
            fprintf(stderr, "Graph must have at least 2 vertices\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    // Broadcast number of vertices and edges
    MPI_Bcast(&numVertices, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numEdges, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Distribute edges among processes
    Edge *localEdges = NULL;
    int localCount = 0;
    distributeEdges(globalEdges, numEdges, &localEdges, &localCount, world_rank, world_size, MPI_EDGE_TYPE);

    // Free global edges on rank 0
    if (world_rank == 0) {
        free(globalEdges);
        globalEdges = NULL;
    }

    // Allocate parent & rank arrays
    int *parent = (int *)malloc(numVertices * sizeof(int));
    int *rank   = (int *)malloc(numVertices * sizeof(int));
    if (!parent || !rank) {
        fprintf(stderr, "Error: Not enough memory for parent/rank arrays\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Prepare cost variable
    long long totalcost = 0;

    boruvkaMST(numVertices, localEdges, localCount,
               parent, rank,
               world_rank, world_size,
               &totalcost, MPI_EDGE_TYPE);

    double end_time = MPI_Wtime();
    double elapsed_time = end_time - start_time;

    // Gather the maximum elapsed time across ranks
    double max_time = 0.0;
    MPI_Reduce(&elapsed_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // Rank 0 writes MST result and timing
    if (world_rank == 0) {
        printf("MST has total cost = %lld\n", totalcost);
        printf("MST Computation Time: %lf seconds (max across ranks)\n", max_time);
    }

    // Cleanup
    free(localEdges);
    free(parent);
    free(rank);

    MPI_Type_free(&MPI_EDGE_TYPE);
    MPI_Finalize();
    return 0;
}
