using TreeTools

# ------------------------------ Functions ------------------------------
function create_genealogy(root_seq, distance_btw_nodes, number_of_nodes, target_number, all_ctc_maps, beta=100)
    """
    Builds a binary tree of depth 'number_of_nodes' using Metropolis evolution.
    Each branch represents 'distance_btw_nodes' mutation steps.
    Returns:
        - newick_str : Newick string compatible with TreeTools.jl
        - seq_dict : Dict{String, Vector{Int}} mapping node name -> sequence
    """
    name_cpt = 0
    seq_dict = Dict{String, Vector{Int}}()  # name -> sequence

    function evolve(seq)
        evolved_seq,_,_ = algo_count_only_if_mut_accepted(copy(seq), target_number, all_ctc_maps,
                                  distance_btw_nodes, beta,
                                  distance_btw_nodes+1, distance_btw_nodes+1, false)
        return evolved_seq
    end

    function build_newick(seq, depth)
        name_cpt += 1
        name = "N$name_cpt"
        seq_dict[name] = copy(seq)

        if depth == 0
            return name
        else
            child1_seq = evolve(seq)
            child2_seq = evolve(seq)

            child1_str = build_newick(child1_seq, depth - 1)
            child2_str = build_newick(child2_seq, depth - 1)

            return "($(child1_str):$(distance_btw_nodes),$(child2_str):$(distance_btw_nodes))$(name)"
        end
    end

    newick_str = build_newick(root_seq, number_of_nodes) * ";"
    return newick_str, seq_dict
end


function distance_matrix(seq_dict)
    names = collect(keys(seq_dict))
    n = length(names)
    D = zeros(Int, n, n)

    for i in 1:n
        for j in i+1:n
            d = hamming_distance(seq_dict[names[i]], seq_dict[names[j]])
            D[i, j] = d
            D[j, i] = d
        end
    end

    return D, names
end


function simple_print_distance_matrix(seq_dict)
    D, names = distance_matrix(seq_dict)
    println("      ", join(rpad.(names, 6))) #6 characters padding
    for i in 1:length(names)
        println(rpad(names[i], 6), join(rpad.(D[i, :], 6)))
    end
end


function plot_distance_matrix(seq_dict)
    D, names = distance_matrix(seq_dict)
    n = length(names)

    p = heatmap(1:n, 1:n, D,
        color = :viridis,
        title = "Hamming Distance between nodes",
        xrotation = 45,
        aspect_ratio = :equal,
        xticks = (1:n, names),
        yticks = (1:n, names),
        yflip = true,
    )

    for i in 1:n, j in 1:n
        annotate!(p, j, i, text(string(D[i, j]), 9, :white))
    end

    return p
end
