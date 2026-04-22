# ------------------------------ Functions ------------------------------
function create_genealogy(root_seq, distance_btw_nodes, number_of_nodes, target_number, all_ctc_maps; beta=100)
    """
    Builds a binary tree of depth 'number_of_nodes' using Metropolis evolution.
    Each branch represents 'distance_btw_nodes' mutation steps.
    Returns:
        - newick_str : Newick string compatible with TreeTools.jl
        - seq_dict : Dict{String, Vector{Int}} mapping node name -> sequence
    """
    name_cpt = 0
    seq_dict = Dict{String,Vector{Int}}()  # name -> sequence

    function evolve(seq)
        evolved_seq, _, _ = algo_count_only_if_mut_accepted(copy(seq), target_number, all_ctc_maps,
            distance_btw_nodes, beta,
            distance_btw_nodes + 1, distance_btw_nodes + 1, false)
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
        color=:viridis,
        title="Hamming Distance between nodes",
        xrotation=45,
        aspect_ratio=:equal,
        xticks=(1:n, names),
        yticks=(1:n, names),
        yflip=true,
    )

    for i in 1:n, j in 1:n
        annotate!(p, j, i, text(string(D[i, j]), 9, :white))
    end

    return p
end



function fitch_reconstruction(tree, seq_dict)
    n_sites = length(first(values(seq_dict)))

    fitch_sets = Dict{String,Vector{Set{Int}}}() #up (Set for no duplicates)
    reconstructed = Dict{String,Vector{Int}}() #down

    #up part (leaves to root)
    for node in postorder_traversal(tree) #post_order = leaves to root
        name = label(node)

        if isleaf(node)
            #Leaf : Fitch set is the singleton {base}
            seq = seq_dict[name]
            fitch_sets[name] = [Set([aa]) for aa in seq]

        else
            #internal node : intersection if non-empty, else union
            children_sets = [fitch_sets[label(c)] for c in children(node)]
            fitch_sets[name] = Vector{Set{Int}}(undef, n_sites)

            for site in 1:n_sites
                #We start from the intersection with the first child's set
                inter = copy(children_sets[1][site])
                for cs in children_sets[2:end] #here it's binary tree, just in case
                    intersect!(inter, cs[site])
                end

                if isempty(inter)
                    #empty intersection -> we take the union of all the children
                    uni = copy(children_sets[1][site])
                    for cs in children_sets[2:end]
                        union!(uni, cs[site])
                    end
                    fitch_sets[name][site] = uni
                else
                    fitch_sets[name][site] = inter
                end
            end
        end
    end

    #down part
    root_node = root(tree)
    root_name = label(root_node)

    #The root arbitrarily chooses an element from its set
    reconstructed[root_name] = [first(s) for s in fitch_sets[root_name]]

    for node in preorder_traversal(tree)
        if isroot(node)
            continue
        end

        name = label(node)
        parent_name = label(ancestor(node))
        parent_seq = reconstructed[parent_name]

        reconstructed[name] = Vector{Int}(undef, n_sites)

        for site in 1:n_sites
            fs = fitch_sets[name][site]

            if parent_seq[site] in fs
                #the parent's aa is in the set -> we keep it (parsimony)
                reconstructed[name][site] = parent_seq[site]
            else
                #otherwise, we take any element from the set
                reconstructed[name][site] = first(fs)
            end
        end
    end

    return reconstructed

end


function reconstruction_error_hamming(reconstructed, seq_dict)
    total_distance = 0
    n = length(reconstructed)

    for (name, recon_seq) in reconstructed
        true_seq = seq_dict[name]
        total_distance += hamming_distance(recon_seq, true_seq)
    end

    return total_distance / n
end


function reconstruction_error_prob_folding(reconstructed, seq_dict, all_ctc_maps, target_number)
    total_distance = 0
    n = length(reconstructed)

    for (name, recon_seq) in reconstructed
        true_seq = seq_dict[name]

        old_proba_folding = proba_of_seq(true_seq, all_ctc_maps, target_number)
        new_proba_folding = proba_of_seq(recon_seq, all_ctc_maps, target_number)
        # println(name)
        # println(old_proba_folding)
        # println(new_proba_folding)
        # println()

        #compute how much the new sequence is better/worse
        #if the difference is positive, it means that the reconstructed seq has a better probability of folding
        total_distance += new_proba_folding - old_proba_folding
    end

    return total_distance / n

end


function plot_reconstruction_error_hamming(root_seq, distance_btw_nodes, number_of_nodes, target_number, all_ctc_maps; beta=100)
    recon_errors = Float64[]
    for d in distance_btw_nodes
        newick_str, seq_dict = create_genealogy(root_seq, d, number_of_nodes, target_number, all_ctc_maps, beta=beta)
        tree = parse_newick_string(newick_str)
        reconstruction_result = fitch_reconstruction(tree, seq_dict)
        error = reconstruction_error_hamming(reconstruction_result, seq_dict)
        push!(recon_errors, error)
    end

    p = plot(distance_btw_nodes, recon_errors,
        xlabel="Distance between nodes",
        ylabel="Average Hamming distance",
        title="Hamming reconstruction error given the distance between nodes",
    )
    return p
end


function plot_reconstruction_error_proba(root_seq, distance_btw_nodes, number_of_nodes, target_number, all_ctc_maps; beta=100)
    recon_errors = Float64[]
    for d in distance_btw_nodes
        newick_str, seq_dict = create_genealogy(root_seq, d, number_of_nodes, target_number, all_ctc_maps, beta=beta)
        tree = parse_newick_string(newick_str)
        reconstruction_result = fitch_reconstruction(tree, seq_dict)
        error = reconstruction_error_prob_folding(reconstruction_result, seq_dict, all_ctc_maps, target_number)
        push!(recon_errors, error)
    end

    p = plot(distance_btw_nodes, recon_errors,
        xlabel="Distance between nodes",
        ylabel="Average folding probability difference",
        title="Probabilistic reconstruction error given the distance between nodes",
    )
    return p
end