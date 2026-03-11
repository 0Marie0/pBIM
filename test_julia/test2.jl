function hamiltonian_path(graph, n, dimension, start_node=1)
    """
    This function find a Hamiltonian path in a given graph using backtracking;
    
    Entry : 
    graph (Dict) = a graph where a key is a given node and the value is the list of the connected nodes
    n (int) = the amount of node on a side
    dimension (int) = the dimension of the graph (2 or 3)
    
    Return : 
    path (Vector) = a Hamiltonian path if it exists, otherwise an empty vector
    """

    #dimension should be between 2 and 3
    if dimension == 2
        total_nodes = n * n
    elseif dimension == 3
        total_nodes = n * n * n
    else
        throw(ArgumentError("dimension should be 2 or 3, but $dim has been provided"))
    end

    paths = Vector{Vector{Int}}() # contains all the possible Hamiltonian path
    path = Int[] # contains the current Hamiltion path

    function backtrack(current_node)
        push!(path, current_node)

        if length(path) == total_nodes
            push!(paths, copy(path))
        end

        for neighbor in graph[current_node]
            if neighbor ∉ path
                backtrack(neighbor)
            end
        end

        pop!(path)
    end

    

    # for start_node in keys(graph)
    #     backtrack(start_node)
    # end


    push!(path, start_node)

    if n == 2
        imposed_first_neighbor = start_node + 1
        backtrack(imposed_first_neighbor)

    elseif n == 3
        imposed_first_neighbor  = start_node + 1
        imposed_second_neighbor = start_node + 2

        push!(path, imposed_first_neighbor)
        backtrack(imposed_second_neighbor)
        pop!(path)
    end

    
    return paths
end

#One simple 2*2 graph

"""
1--2
|  |
3--4

graph = Dict(
    1 => [2, 3],
    2 => [1, 4],
    3 => [1, 4],
    4 => [2, 3]
)

n = 2

path = hamiltonian_path(graph, n, 2)
println("Hamiltonian paths found: ", path)
"""

function create_graph2D(n)
    """
    This function create a graph that represents the links between each node of a 2D lattice;
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

"""
#test for 2D graph 3*3
n = 3
graph = create_graph2D(n)
print(graph)
print("\n")
path = hamiltonian_path(graph, n, 2)
println("Hamiltonian paths found: ", path)
print("Number of Hamiltonian paths found: ", length(path)) 
# On ne devrait avoir que 3 chemins. 
"""

function create_graph3D(n)
    """
    Create a graph representing the links between each node of a 3D lattice.
    Entry:
    n (int) = the amount of node on a side
    Return : 
    Graph (Dict) = a graph where a key is a given node and the value is the list of the connected nodes
    """
    #linear numerotation of nodes 
    idx(i,j,k) = (i-1)*n^2 + (j-1)*n + k
    graph = Dict{Int, Vector{Int}}()

    for i in 1:n, j in 1:n, k in 1:n
        neighbors = Int[]

        #neighbors on x-axis
        if i > 1
            push!(neighbors, idx(i-1,j,k))
        end
        if i < n
            push!(neighbors, idx(i+1,j,k))
        end

        #neighbors on y-axis
        if j > 1
            push!(neighbors, idx(i,j-1,k))
        end
        if j < n
            push!(neighbors, idx(i,j+1,k))
        end

        #neighbors on z-axis
        if k > 1
            push!(neighbors, idx(i,j,k-1))
        end
        if k < n
            push!(neighbors, idx(i,j,k+1))
        end

        graph[idx(i,j,k)] = sort!(neighbors)
    end

    return graph
end


#3D graph 2*2*2
n = 2
graph3D = create_graph3D(n)
println(graph3D)
#test for 3D graph 2*2*2
path = hamiltonian_path(graph3D, n, 3)
println("Hamiltonian paths found: ", path)
print("Number of Hamiltonian paths found: ", length(path)) 
# On ne devrait avoir que 5 chemins, [1, 2, 4, 3, 7, 8, 6, 5] et [1, 2, 6, 5, 7, 8, 4, 3] sont identiques par rotation

#Doesn't work for now, i think they are too many paths, 
n = 3
graph3D = create_graph3D(n)
println(graph3D)

#test for 3D graph 3*3*3
path = hamiltonian_path(graph3D, n, 3)
#println("Hamiltonian paths found: ", path)
print("Number of Hamiltonian paths found: ", length(path))

path

#bis function just to know the different way to start in the medium layer of a 3*3*3
#So we return every start of path from a specific node of a 3*3 graph
function hamiltonian_path_bis(graph, n, dimension, first_node, second_node)

    #dimension should be between 2 and 3
    if dimension == 2
        total_nodes = n * n
    elseif dimension == 3
        total_nodes = n * n * n
    else
        throw(ArgumentError("dimension should be 2 or 3, but $dim has been provided"))
    end

    paths = Vector{Vector{Int}}() #contains all the possible Hamiltonian path
    path = Int[] #contains the current Hamiltion path

    function backtrack(current_node)
        push!(path, current_node)

        push!(paths, copy(path)) #ici on s'en fiche de ne pas avoir des paths finis

        # if length(path) == total_nodes
        #     push!(paths, copy(path))
        # end

        for neighbor in graph[current_node]
            if neighbor ∉ path
                backtrack(neighbor)
            end
        end

        pop!(path)
    end


    # for start_node in keys(graph)
    #     backtrack(start_node)
    # end


    imposed_first_neighbor = first_node
    imposed_second_neighbor = second_node

    push!(path, imposed_first_neighbor)
    backtrack(imposed_second_neighbor)

    return paths
end

  


function create_graph2D_bis()

    #manual mapping of positions
    grid = [
        2   11   20;
        5   14   23;
        8   17   26
    ]

    graph = Dict{Int, Vector{Int}}()

    for i in 1:3, j in 1:3
        node = grid[i, j]
        neighbors = Int[]

        # Up
        if i > 1
            push!(neighbors, grid[i-1, j])
        end

        # Down
        if i < 3
            push!(neighbors, grid[i+1, j])
        end

        # Left
        if j > 1
            push!(neighbors, grid[i, j-1])
        end

        # Right
        if j < 3
            push!(neighbors, grid[i, j+1])
        end

        graph[node] = sort!(neighbors)
    end

    return graph
end


n = 3
graph = create_graph2D_bis()
println(graph)

paths_1 = hamiltonian_path_bis(graph, n, 2, 5, 2)
println("Hamiltonian paths found: ", paths_1)
println("Number of Hamiltonian paths found: ", length(paths_1)) 

paths_2 = hamiltonian_path_bis(graph, n, 2, 5, 14)
println("Hamiltonian paths found: ", paths_2)
println("Number of Hamiltonian paths found: ", length(paths_2)) 

#liste exaustive des paths de départ à définir

for path in paths_1
    push!(path, path[end]-1)
end

for path in paths_2
    push!(path, path[end]-1)
end

paths_corner = [[1,2,3,6],[1,2,5]]

paths_final = vcat(paths_1, paths_2, paths_corner)


println("Hamiltonian paths begginings found: ", paths_final)
println("Number of Hamiltonian paths begginings found: ", length(paths_final)) 


function hamiltonian_path_with_beginning(graph, n, dimension, beginning)
    #dimension should be between 2 and 3
    if dimension == 2
        total_nodes = n * n
    elseif dimension == 3
        total_nodes = n * n * n
    else
        throw(ArgumentError("dimension should be 2 or 3, but $dim has been provided"))
    end

    paths = Vector{Vector{Int}}() # contains all the possible Hamiltonian path
    path = Int[] # contains the current Hamiltion path

    function backtrack(current_node)
        push!(path, current_node)

        if length(path) == total_nodes
            push!(paths, copy(path))
        end

        for neighbor in graph[current_node]
            if neighbor ∉ path
                backtrack(neighbor)
            end
        end

        pop!(path)
    end


    append!(path, beginning[1:end-1])
    start_node = beginning[end]
    backtrack(start_node)
    
    return paths
end

"""
total_paths = []

for beginning in paths_final
    paths = hamiltonian_path_with_beginning(graph3D, n, 3, beginning)
    append!(total_paths,paths)
end

print("Number of Hamiltonian paths found: ", length(total_paths))
"""

#55020
#120458
#88875
"""
println(total_paths[1])
println(total_paths[end])
"""
# paths_final = [[1,2,3],[1,2,5],[5,2,3],[5,14,11],[5,2,1]]

# total_paths = []

# for beginning in paths_final
#     paths = hamiltonian_path_with_beginning(graph3D, n, 3, beginning)
#     append!(total_paths,paths)
# end

# print("Number of Hamiltonian paths found: ", length(total_paths))

using StaticArrays
using DataStructures

#cube 3*3*3
const N = 3
const TOTAL_NODES = N^3


#ID (1-27) to Coordinates (0-2, 0-2, 0-2)
function id_to_coords(id)
    z = div(id - 1, 9) #id - 1 because Julia indices start at 1
    rem = mod(id - 1, 9) #remainder of the entire division
    y = div(rem, 3)
    x = mod(rem, 3)
    return (x, y, z)
end


#Coord to ID
function coords_to_id(x, y, z)
    return x + 3y + 9z + 1
end


#Generating the table of 24 rotations (only det = 1)
#Reflections are excluded when counting chiral forms separately
function generate_rotation_table()
    #permutation table : table[id, r] = new_id -> under rotation r, node id becomes new_id
    table = zeros(Int, 27, 48) #27 -> each node ; 24 -> each rotation
    
    #6 possibles transposition : [1,2,3] -> (x,y,z)  ;  [1,3,2] -> (x,z,y)  ;  [2,1,3] -> (y,x,z)  etc.
    perms = [[1,2,3], [1,3,2], [2,1,3], [2,3,1], [3,1,2], [3,2,1]]
    
    #Each combination corresponds to a change of sign on the axes : total 2^3 = 8 sign combination
    signs = [[s1, s2, s3] for s1 in [-1,1], s2 in [-1,1], s3 in [-1,1]]
    
    #so we got 6*8 = 48 rotations
    #we only keep 24 rotations because we want to keep mirror inversions (so mirror are not considered as rotations)

    rot_idx = 1
    for p in perms, s in signs
        #Calculating the determinant to retain only the physical rotations
        
        #a transposition is even if the number of permutation is even, determinant of even=1, else -1
        p_parity = (p == [1,2,3] || p == [2,3,1] || p == [3,1,2]) ? 1 : -1
        
        #The transformation is a permutation of the axes followed by sign reversals
        #Its determinant is det(P)*det(D) = parity(permutation) × (s1*s2*s3)
        #This comes from the general property: det(AB) = det(A)*det(B)
        #Since the sign matrix is ​​diagonal, its determinant is the product of its diagonal
        det = p_parity * (s[1] * s[2] * s[3])

        #We only keep det = +1 to preserve the orientation (pure rotations)
        #if det == 1
            for id in 1:27 #rotation applied to all points
                c = id_to_coords(id) .- 1 #Centering on 0 for rotation, (0,1,2) becomes (-1,0,1)
                new_c = (c[p[1]]*s[1], c[p[2]]*s[2], c[p[3]]*s[3]) #application of the permutation
                table[id, rot_idx] = coords_to_id(new_c[1]+1, new_c[2]+1, new_c[3]+1)
            end
            rot_idx += 1 #next rotation
        #end
    end
    return table
end

const ROT_TABLE = generate_rotation_table()


#Canonical formatting function
function get_canonical(path)
    best_path = copy(path)
    #similar(path) create an object of same type and same lenght but empty
    current_variant = similar(path) #to temporarily store an applied rotation
    rev_variant = similar(path) #to temporarily store the reversed version of the path
    
    for s in 1:48 #each rotation
        #application of rotation
        for i in 1:27 #each node
            current_variant[i] = ROT_TABLE[path[i], s] #new id after rotation
        end
        
        #We are looking for the minimal form among all the rotations
        if current_variant < best_path #lexicographical comparison ; [1,3,9] < [1,4,2] because 1==1 then 3<4 
            best_path .= current_variant #copy element by element (.)
        end
        
        #We reverse the reading direction of the chain
        # for i in 1:27
        #     rev_variant[i] = current_variant[28-i]
        # end
        # if rev_variant < best_path #lexicographical comparison
        #     best_path .= rev_variant
        # end
    end
    return best_path
end


#Connectivity graph
function build_graph()
    adj = [Int[] for _ in 1:27] #init of 27 empty vectors
    for id in 1:27
        x, y, z = id_to_coords(id)
        #we consider the 6 directions in the cube (+x,-x,+y,-y,+z,-z)
        for (dx, dy, dz) in [(1,0,0), (-1,0,0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1)]
            nx, ny, nz = x+dx, y+dy, z+dz
            if 0 <= nx < 3 && 0 <= ny < 3 && 0 <= nz < 3 #if the neighbor exist
                push!(adj[id], coords_to_id(nx, ny, nz)) #add the neighbor to the vector
            end
        end
    end
    return adj
end

const ADJ = build_graph()


#Backtracking with filtering by Set
function find_all_unique_conformations()
    #the Set guarantees uniqueness
    unique_paths = OrderedSet{Vector{Int8}}()
    visited = zeros(Bool, 27)
    current_path = zeros(Int8, 27)

    #Recursive function to explore all Hamiltonian paths
    function backtrack(node, depth)
        visited[node] = true #indicates that this node is explored
        current_path[depth] = node #we add the node to the current path
        
        if depth == 27
            #we transform the path into its unique reference version
            push!(unique_paths, get_canonical(current_path)) #if already in Set, not added
        else
            for neighbor in ADJ[node]
                if !visited[neighbor] #if this neighbor has not been visited
                    backtrack(neighbor, depth + 1) #backtrack with this new node
                end
            end
        end
        #after exploring all paths starting from this node, it is freed up for other paths
        visited[node] = false 
    end

    println("Calculation of all conformations")
    
    println("Step 1: Exploration from the corner (node 1)")
    backtrack(1, 1) 
    
    println("Step 2: Exploration from the center of the face (node 5)")
    backtrack(5, 1) 
    
    return unique_paths,length(unique_paths)
end

#script launch
paths, count = find_all_unique_conformations()
println("Final result : $count")

println(first(paths))
println(typeof(paths))

function paths_to_file(filename, paths)
    open(filename, "w") do file
        for path in paths
            println(file, join(path, ",")) 
        end
    end
end

paths_to_file("paths_bis.csv", paths)

function all_edges_3x3x3()
    """
    We retreive all the edges from a square : we have 54 tuples. 
    """
    edges = Set{Tuple{Int,Int}}()

    for id in 1:27
        for neighbor in ADJ[id]
            edge = (min(id, neighbor), max(id, neighbor)) #to avoid duplicates (1,2) and (2,1)
            push!(edges, edge)
        end
    end

    return edges, length(edges)
end

#script launch
edges, count = all_edges_3x3x3()
const ALL_EDGES = edges

println("Final result : $count edges")
println(edges)

function contact_map(path)
    """
    Retreive all the contact for a given path
    """
    in_path = Set{Tuple{Int,Int}}() # Set with all the edges used in the path
    
    for k in 1:length(path)-1
        a, b = path[k], path[k+1]
        a, b = a < b ? (a,b) : (b,a) # avoid duplicates (such as (1,2) and (2,1))
        push!(in_path, (a, b))
    end

    return setdiff(ALL_EDGES, in_path) # We only keep "non-used" eedges

end

contact_map(first(paths))

const aa = ["CYS", "MET", "PHE", "ILE", "LEU", "VAL", "TRP", "TYR", "ALA", "GLY", "THR", "SER", "GLN", "ASN", "GLU", "ASP", "HIS", "ARG", "LYS"]
const aa_idx = Dict(i => aa[i] for i in 1:length(aa)) #i changed the logic of the dictionnary, now key=1:20 and value=aa

aa_idx[1] # -> CYS

using CSV
using DataFrames

# Reading the file and converting into a matrix 
df = CSV.read("MJ.csv", DataFrame; delim='\t', header=false)
MJ = Matrix{Float64}(df)

# Reading the value for MET and ILE in the MJ matrix
MJ[2, 20]

using Random

#Generate a random 27 a.a sequence
function generate_seq()
    return rand(1:20, 27)
end

#Compute the total energy of one sequence for one structure
function energy_one_struct(ctc_map, seq)
    score=0
    for contact in ctc_map
        pos1,pos2 = contact
        aa1 = seq[pos1]
        aa2 = seq[pos2]
        score+=MJ[aa1, aa2]
    end
    return score
end

function all_contact_maps(paths)
    ctc_maps=[]
    for path in paths
        push!(ctc_maps, contact_map(path))
    end
    return ctc_maps
end

#compute the denominator = sum of all the exponential of minus the energy 
function total_energy(ctc_maps, seq)
    total=0
    for ctc_map in ctc_maps
        total+=exp(-energy_one_struct(ctc_map,seq))
    end
    return total
end


#Compute the score/probability of one struture for one sequence among all the other structures for all the structures
function probability_each_structure(ctc_maps, seq, total_energy)
    proba=[]
    for ctc_map in ctc_maps
        score=energy_one_struct(ctc_map,seq)
        p = exp(-score)/total_energy
        push!(proba,p)
    end
    return proba
end

function probability_target_structure(ctc_maps, seq, total_energy,target_number)
    s=energy_one_struct(ctc_maps[target_number], seq)
    return exp(-s)/total_energy
end

#random seq
seq = generate_seq()
println(seq)

#choose a random structure
target_number = 75005
target_structure = paths[target_number]
#println(target_structure)
#for target_number=75005
#[5, 2, 1, 4, 7, 16, 25, 26, 27, 18, 15, 24, 21, 12, 3, 6, 9, 8, 17, 14, 13, 10, 19, 22, 23, 20, 11]

all_ctc_maps = all_contact_maps(paths)

ttl_energy=total_energy(all_ctc_maps,seq)
println(ttl_energy)

probas=probability_each_structure(all_ctc_maps, seq, ttl_energy)
prob=probas[target_number]
println(prob)

function mutation(seq)
    pos = rand(1:27)
    current_aa = seq[pos]

    #new aa different from actual
    new_aa = rand(1:19)
    if new_aa >= current_aa #we skip if the same aa
        new_aa += 1
    end

    return (pos, new_aa)
end

"""

function metropolis_non_fonctionnelle(mut, seq, target_number, all_ctc_maps, beta=1000)
    old_score=energy_one_struct(all_ctc_maps[target_number], seq)
    new_seq = copy(seq)
    pos,new_aa = mut
    new_seq[pos]=new_aa
    new_score=energy_one_struct(all_ctc_maps[target_number], new_seq)

    #println(old_score, new_score)

    #because smaller score = better sequence (p = exp(-score)/total_energy)
    if new_score<=old_score
        #println("Mutation accepted")
        return new_seq
    else
        #not exactly the reel probs because we should recompute ttl_energy for the new seq but same result/order of magnitude
        delta = old_score-new_score #new_score>old_score so delta is <0
        limit = exp(beta * delta) #because its a difference between score and not probability, the gap can be big so we use a little delta (10 for exemple)
        #println(limit)
        val=rand()
        #println(val)
        if val <= limit
            println("test, moins bonne mutation acceptée")
            return new_seq
        end
    end
    #println("No mutation")
    return seq
end

"""

function metropolis(mut, seq, target_number, all_ctc_maps, ttl_energy, old_prob, i, beta=1000)
    old_log_proba=log(old_prob)
    
    new_seq = copy(seq)
    pos,new_aa = mut
    new_seq[pos]=new_aa

    new_ttl_energy=total_energy(all_ctc_maps, new_seq)
    new_prob=probability_target_structure(all_ctc_maps, new_seq, new_ttl_energy, target_number)

    new_log_proba=log(new_prob)

    # println("new log proba $new_log_proba")
    # println("old log proba $old_log_proba")

    if new_log_proba>=old_log_proba
        #println("meilleure sequence acceptée, epoch $i")
        return new_seq, new_ttl_energy, new_prob
    else
        delta = new_log_proba-old_log_proba
        limit = exp(beta * delta)
        val=rand()
        #println(val)
        if val <= limit
            #println("test, moins bonne mutation acceptée, epoch $i")
            return new_seq, new_ttl_energy, new_prob
        end
    end
    #println("moins bonnne séquence rejetée, pas de modif, epoch $i")
    return seq, ttl_energy, old_prob
end

"""

function algo_non_fonctionnel(seq, target_number, all_ctc_maps, epochs=100)
    actual_seq=copy(seq)
    mut_last_iter=false

    epochs=epochs
    p_struct_time_t=[]
    ttl_energy=total_energy(all_ctc_maps,actual_seq)
    probas=probability_each_structure(all_ctc_maps, actual_seq, ttl_energy)

    for i in 1:epochs
        if mut_last_iter
            ttl_energy=total_energy(all_ctc_maps,actual_seq)
            probas=probability_each_structure(all_ctc_maps, actual_seq, ttl_energy)
        end
        push!(p_struct_time_t, probas[target_number])

        mut=mutation(actual_seq)
        seq_metropolis = metropolis(mut, actual_seq, target_number, all_ctc_maps, 10)
        if (seq_metropolis!=actual_seq)
            mut_last_iter=true
            actual_seq=seq_metropolis
        else
            mut_last_iter=false
        end
    end
    return actual_seq, p_struct_time_t
end    

"""

function algo(seq, target_number, all_ctc_maps, epochs=100, beta=100)
    actual_seq=copy(seq)
    mut_last_iter=false

    epochs=epochs
    p_struct_time_t=[]
    new_ttl_energy=total_energy(all_ctc_maps,actual_seq)
    new_prob=probability_target_structure(all_ctc_maps, actual_seq, new_ttl_energy,target_number)

    for i in 1:epochs
        if i%100==0
            println("epoch number $i")
        end
        # if mut_last_iter
        #     ttl_energy=total_energy(all_ctc_maps,actual_seq)
        #     probas=probability_each_structure(all_ctc_maps, actual_seq, ttl_energy)
        # end
        push!(p_struct_time_t, new_prob)

        mut=mutation(actual_seq)
        seq_metropolis, new_ttl_energy, new_prob = metropolis(mut, actual_seq, target_number, all_ctc_maps, new_ttl_energy, new_prob, i, beta)
        if (seq_metropolis!=actual_seq)
            #println("test, changement de sequence a l'iter $i")
            mut_last_iter=true
            actual_seq=seq_metropolis
        else
            mut_last_iter=false
        end
    end
    return actual_seq, p_struct_time_t
end    

#final_seq, p_struct_time_t = algo(seq, target_number, all_ctc_maps, epochs=1000, beta=10)

# using Plots

# t = 1:length(p_struct_time_t)

# plot(t, p_struct_time_t, xlabel="Time", ylabel="Probability", title="Evolution of the probability of the sequence over time")

final_seq, p_struct_time_t = algo(seq, target_number, all_ctc_maps, 100, 100)

using Plots

t = 1:length(p_struct_time_t)

plot(t, p_struct_time_t, xlabel="Time", ylabel="Probability", title="Evolution of the probability of the sequence over time")

using Profile
using ProfileView
using PProf

#First JIT 
algo(seq, target_number, all_ctc_maps, 10, 100)

#reset profiler then real mesure
Profile.clear()
@profile algo(seq, target_number, all_ctc_maps, 500, 100)

ProfileView.view()
