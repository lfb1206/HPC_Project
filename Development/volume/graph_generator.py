import random

def generate_weighted_graph(vertices, edges, filename):
    if edges < vertices - 1:
        raise ValueError("Number of edges is too small to form a connected graph.")
    if edges > vertices * (vertices - 1) // 2:
        raise ValueError("Too many edges for the given number of vertices.")
    
    edge_set = set()
    
    # Step 1: Add edges to form a spanning tree (ensuring connectivity)
    for i in range(1, vertices):
        u = i
        v = random.randint(0, i - 1)  # Connect to a random previous node
        weight = random.randint(1, 100)
        edge_set.add((u, v, weight))
    
    # Step 2: Add remaining edges randomly
    while len(edge_set) < edges:
        u = random.randint(0, vertices - 1)
        v = random.randint(0, vertices - 1)
        if u != v:
            weight = random.randint(1, 100)
            edge_set.add((min(u, v), max(u, v), weight))
    
    # Write to the file
    with open(filename, 'w') as file:
        file.write(f"{vertices} {edges}\n")
        for u, v, weight in edge_set:
            file.write(f"{u} {v} {weight}\n")

if __name__ == "__main__":
    num_vertices = int(input("Enter the number of vertices: "))
    degree = int(input("Enter the desired average vertex degree: "))
    base_edges = num_vertices * degree // 2
    
    # Variability of ±5%
    variability = int(base_edges * 0.05)
    num_edges = max(base_edges - variability, num_vertices - 1)  # Ensure enough edges for connectivity
    num_edges = min(num_edges, base_edges + variability)
    
    # Calculate the actual average degree
    actual_avg_degree = (2 * num_edges) / num_vertices
    
    output_file = "graph.txt"

    try:
        generate_weighted_graph(num_vertices, num_edges, output_file)
        print(f"Graph of {num_vertices} nodes and {num_edges} edges written to {output_file}")
        print(f"Expected average degree: {degree}, Actual average degree: {actual_avg_degree:.2f}")
    except ValueError as e:
        print(f"Error: {e}")
