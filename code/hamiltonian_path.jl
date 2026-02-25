# ------------------------------ Functions ------------------------------
using StaticArrays
using DataStructures
using CSV
using DataFrames
using Plots
plotlyjs()


# Functions to find all Hamiltonion paths in 3*3*3 cube
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

        for id in 1:27 #rotation applied to all points
            c = id_to_coords(id) .- 1 #Centering on 0 for rotation, (0,1,2) becomes (-1,0,1)
            new_c = (c[p[1]]*s[1], c[p[2]]*s[2], c[p[3]]*s[3]) #application of the permutation
            table[id, rot_idx] = coords_to_id(new_c[1]+1, new_c[2]+1, new_c[3]+1)
        end
        rot_idx += 1 #next rotation
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
paths, _ = find_all_unique_conformations()
const PATHS = paths


# Function for saving paths 
function paths_to_file(filename)
    """
    Save the paths to a CSV file, each path is a row and nodes are separated by commas.
    """
    open(filename, "w") do file
        for path in PATHS
            println(file, join(path, ",")) 
        end
    end
end

# Function for contact map
function all_edges_3x3x3()
    """
    Retreive all the edges of the 3x3x3 cube.
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

edges, _ = all_edges_3x3x3()
const ALL_EDGES = edges


function contact_map(path)
    """
    Retreive all the contact for a given path.
    """
    in_path = Set{Tuple{Int,Int}}() # Set with all the edges used in the path
    
    for k in 1:length(path)-1
        a, b = path[k], path[k+1]
        a, b = a < b ? (a,b) : (b,a) # avoid duplicates (such as (1,2) and (2,1))
        push!(in_path, (a, b))
    end
    return setdiff(ALL_EDGES, in_path) # We only keep "non-used" eedges

end

# Function for visualizing the contact map
function get_nodes(range = 0:2)
    # Iterators.product crée un itérateur paresseux (lazy)
    # vec() transforme le tout en vecteur sans copier inutilement
    points = vec(collect(Iterators.product(range, range, range)))
    
    # On sépare les composants proprement
    x = [p[1] for p in points]
    y = [p[2] for p in points]
    z = [p[3] for p in points]
    
    return x, y, z
end


function get_edges(range = 0:2)
    start, stop = first(range), last(range)
    n = length(range)
    # On pré-alloue pour éviter les redimensionnements dynamiques (3 points par arête * 3 axes * n^2)
    total_elements = 3 * 3 * n^2
    lx = Vector{Float64}(undef, total_elements)
    ly = Vector{Float64}(undef, total_elements)
    lz = Vector{Float64}(undef, total_elements)

    idx = 1
    for i in range, j in range
        # Axe X
        lx[idx:idx+2] .= [start, stop, NaN]; ly[idx:idx+2] .= [i, i, NaN]; lz[idx:idx+2] .= [j, j, NaN]
        idx += 3
        # Axe Y
        lx[idx:idx+2] .= [i, i, NaN]; ly[idx:idx+2] .= [start, stop, NaN]; lz[idx:idx+2] .= [j, j, NaN]
        idx += 3
        # Axe Z
        lx[idx:idx+2] .= [i, i, NaN]; ly[idx:idx+2] .= [j, j, NaN]; lz[idx:idx+2] .= [start, stop, NaN]
        idx += 3
    end
    return lx, ly, lz
end


const indices_hover = [(coords_to_id(x,y,z)) for (x, y, z) in zip(x, y, z)]

function plot_cube_ctc_map(x, y, z, lx, ly, lz, edges_to_highlight)
    # ----- skeleton -----
    p = plot(lx, ly, lz, color=:black, alpha=0.5, lw=1, label="")
    
    # ----- Highlighted_edges -----
    hx, hy, hz = Float64[], Float64[], Float64[]
    for (n1, n2) in edges_to_highlight
        c1 = id_to_coords(n1)
        c2 = id_to_coords(n2)
        
        append!(hx, [c1[1], c2[1], NaN])
        append!(hy, [c1[2], c2[2], NaN])
        append!(hz, [c1[3], c2[3], NaN])
    end
    
    # plotting ctcmap over skeleton
    plot!(p, hx, hy, hz, color=:red, lw=3, label="Contact")

    # ----- nodes -----
    scatter!(p, x, y, z, 
            markersize=4, 
            markercolor=:blue,
            xlims = (0, 2),
            ylims = (0, 2),
            zlims = (0, 2),
            hover = indices_hover,
            showaxis=false,
            grid=false,
            ticks=false,
            aspect_ratio=:equal)

    return p
end



# ------------------------------ Main ------------------------------

# Finding all paths 
paths, count = find_all_unique_conformations()
println("We find $count unique paths.")
println(first(paths))
println(typeof(paths))

# Saving paths to a file
paths_to_file("paths.csv")
println("Paths saved to paths.csv")

# Finding the contact map of a given path
println("The first path is : \n", first(paths))
println("Its contact map is : \n", contact_map(first(paths)))

# Visualization of the contact map
const x, y, z = get_nodes()
const lx, ly, lz = get_edges()
liste_contact = contact_map(first(paths))
p = plot_cube_ctc_map(x, y, z, lx, ly, lz, liste_contact)
display(p)