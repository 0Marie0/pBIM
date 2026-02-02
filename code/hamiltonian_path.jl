function create_graph(n)
    """
    This function create a graph that represents the links between each node of the lattice;
    Entry : 
    n (int) = the amount of node on a side
    Return : 
    Graph (Dict) = a graph where a key is a given node and the value is the list of the connected nodes
    """
    idx(i,j) = (i-1)*n + j
    graph = Dict{Int, Vector{Int}}()

    for i in 1:n, j in 1:n
        neighbors = Int[]

        if i > 1
            push!(neighbors, idx(i-1,j)) # Up
        end
        if i < n
            push!(neighbors, idx(i+1,j)) # Down
        end
        if j > 1
            push!(neighbors, idx(i,j-1)) # Left
        end
        if j < n
            push!(neighbors, idx(i,j+1)) # Right
        end

        graph[idx(i,j)] = sort!(neighbors)
    end

    return graph
end


# Usage:
n = 3
graph = create_graph(n)
print(graph)