#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Struct to store an edge
typedef struct {
    long long u;
    long long v;
    int weight;
} Edge;

void generate_weighted_graph(long long vertices, long long edges, const char *filename) {
    if (edges < vertices - 1) {
        fprintf(stderr, "Number of edges is too small to form a connected graph.\n");
        exit(EXIT_FAILURE);
    }
    if (edges > vertices * (vertices - 1) / 2) {
        fprintf(stderr, "Too many edges for the given number of vertices.\n");
        exit(EXIT_FAILURE);
    }

    // Allocate memory for the edge list
    Edge *edge_list = malloc(edges * sizeof(Edge));
    if (!edge_list) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(EXIT_FAILURE);
    }

    srand(time(NULL));

    // Step 1: Add edges to form a spanning tree (ensure connectivity)
    long long edge_count = 0;
    for (long long i = 1; i < vertices; i++) {
        long long u = i;
        long long v = rand() % i; 
        int weight = (rand() % 100) + 1;
        edge_list[edge_count++] = (Edge){u, v, weight};
    }

    // Step 2: Add remaining edges randomly
    while (edge_count < edges) {
        long long u = rand() % vertices;
        long long v = rand() % vertices;
        if (u != v) { // Ensure no self-loops
            int weight = (rand() % 100) + 1;
            edge_list[edge_count++] = (Edge){u, v, weight};
        }
    }

    // Write to file
    FILE *file = fopen(filename, "w");
    if (!file) {
        fprintf(stderr, "Error opening file.\n");
        free(edge_list);
        exit(EXIT_FAILURE);
    }

    fprintf(file, "%lld %lld\n", vertices, edges);
    for (long long i = 0; i < edge_count; i++) {
        fprintf(file, "%lld %lld %d\n", edge_list[i].u, edge_list[i].v, edge_list[i].weight);
    }

    fclose(file);
    free(edge_list);
}

int main() {
    long long num_vertices, degree, base_edges, variability, num_edges;
    printf("Enter the number of vertices: ");
    scanf("%lld", &num_vertices);
    printf("Enter the desired average vertex degree: ");
    scanf("%lld", &degree);

    base_edges = num_vertices * degree / 2;
    variability = (long long)(base_edges * 0.05);
    num_edges = base_edges > variability ? base_edges - variability : num_vertices - 1;
    num_edges = num_edges < base_edges + variability ? num_edges : base_edges + variability;

    double actual_avg_degree = (2.0 * num_edges) / num_vertices;
    const char *output_file = "graph.txt";

    generate_weighted_graph(num_vertices, num_edges, output_file);
    printf("Graph of %lld nodes and %lld edges written to %s\n", num_vertices, num_edges, output_file);
    printf("Expected average degree: %lld, Actual average degree: %.2f\n", degree, actual_avg_degree);

    return 0;
}
