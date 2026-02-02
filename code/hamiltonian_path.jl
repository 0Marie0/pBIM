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

function hamiltonian_path(graph, n)
    """
    This function find a Hamiltonian path in a given graph using backtracking;
    Entry : 
    graph (Dict) = a graph where a key is a given node and the value is the list of the connected nodes
    n (int) = the amount of node on a side
    Return : 
    path (Vector) = a Hamiltonian path if it exists, otherwise an empty vector
    """
    total_nodes = n * n
    path = Int[]

    function backtrack(current_node)
        push!(path, current_node)

        if length(path) == total_nodes
            return true
        end

        for neighbor in graph[current_node]
            if neighbor ∉ path
                if backtrack(neighbor)
                    return true
                end
            end
        end

        pop!(path)
        return false
    end

    for start_node in keys(graph)
        if backtrack(start_node)
            return path
        end
    end

    return Int[]
end

# Example usage:
n = 3
graph = create_graph(n)
print(graph)