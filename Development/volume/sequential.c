#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <time.h>

typedef struct Edge {
    int src, dest, weight;
    int covered; 
} Edge;

typedef struct Graph {
    int V, E;
    Edge* edges;
} Graph;

typedef struct Subset {
    int parent;
    int rank;
} Subset;

Graph* createGraph(int V, int E) {
    Graph* graph = (Graph*)malloc(sizeof(Graph));
    graph->V = V;
    graph->E = E;
    graph->edges = (Edge*)malloc(E * sizeof(Edge));
    int i;
    for (i = 0; i < E; i++) {
        graph->edges[i].covered = 0; 
    }
    return graph;
}

int find(Subset subsets[], int v) {
    if (subsets[v].parent != v) {
        subsets[v].parent = find(subsets, subsets[v].parent);
    }
    return subsets[v].parent;
}

void componentUnion(Subset subsets[], int x, int y) {
    int xroot = find(subsets, x);
    int yroot = find(subsets, y);

    // Union by rank to optimize find
    if (xroot != yroot) {
        if (subsets[xroot].rank < subsets[yroot].rank)
            subsets[xroot].parent = yroot;
        else if (subsets[xroot].rank > subsets[yroot].rank)
            subsets[yroot].parent = xroot;
        else {
            subsets[yroot].parent = xroot;
            subsets[xroot].rank++;
        }
    }
}

void boruvkaMST(Graph* graph) {
    int V = graph->V;
    int E = graph->E;
    Edge* edges = graph->edges;

    Subset* subsets = (Subset*)malloc(V * sizeof(Subset));
    int* cheapest = (int*)malloc(V * sizeof(int));

    // Initialize arrays
    int v;
    for (v = 0; v < V; ++v) {
        subsets[v].parent = v;
        subsets[v].rank = 0;
        cheapest[v] = -1;
    }

    int numTrees = V;
    int MSTweight = 0;

    int iteration = 0;

    while (numTrees > 1) {
        // Reset array
        memset(cheapest, -1, V * sizeof(int));

        // Find cheapest edge, if covered is skipped because it doesn't lead to a new component
        int i;
        for (i = 0; i < E; i++) {
            if (edges[i].covered) {
                continue; 
            }

            int u = edges[i].src;
            int v = edges[i].dest;
            int w = edges[i].weight;

            int set1 = find(subsets, u);
            int set2 = find(subsets, v);

            if (set1 == set2) {
                edges[i].covered = 1;
            } else {
                if (cheapest[set1] == -1 || edges[cheapest[set1]].weight > w)
                    cheapest[set1] = i;
                if (cheapest[set2] == -1 || edges[cheapest[set2]].weight > w)
                    cheapest[set2] = i;
            }
        }

        // Union of components with cheapest edge
        for (v = 0; v < V; v++) {
            if (cheapest[v] != -1) {
                int x = cheapest[v];
                int u = edges[x].src;
                int v = edges[x].dest;
                int w = edges[x].weight;

                int set1 = find(subsets, u);
                int set2 = find(subsets, v);

                if (set1 != set2) {
                    MSTweight += w;
                    componentUnion(subsets, set1, set2);
                    numTrees--;
                    edges[x].covered = 1; // Mark MST edge as covered
                }
            }
        }

        iteration++;
    }

    printf("Total weight of MST is %d\n", MSTweight);

    free(subsets);
    free(cheapest);
}

int main() {
    // Start timing
    clock_t start = clock();

    /* 
    File format should be:
        line 1: <num_vertices> <num_edges>
        next lines: <src> <dst> <weight>
    */
    FILE* file = fopen("graph.txt", "r");
    if (!file) {
        fprintf(stderr, "Error: Unable to open file\n");
        return 1;
    }

    int V, E;
    fscanf(file, "%d %d", &V, &E);
    Graph* graph = createGraph(V, E);

    // Reads and store Edges on the defined structure
    int i;
    for (i = 0; i < E; i++) {
        fscanf(file, "%d %d %d", &graph->edges[i].src, &graph->edges[i].dest, &graph->edges[i].weight);
    }

    fclose(file);

    boruvkaMST(graph);

    // End timing
    clock_t end = clock();

    // Calculate elapsed time in seconds
    double elapsed_time = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time taken: %.6f seconds\n", elapsed_time);

    free(graph->edges);
    free(graph);

    return 0;
}