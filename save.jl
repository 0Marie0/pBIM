# ── Précalcul (à faire une seule fois avant l'algo) ──────────────────────────

function precompute_pos_contacts(all_ctc_maps)
    # pos_contacts[p] = liste de (structure_idx, partenaire) pour la position p
    pos_contacts = [Vector{Tuple{Int,Int}}() for _ in 1:27]
    for (s, ctc_map) in enumerate(all_ctc_maps)
        for (a, b) in ctc_map
            push!(pos_contacts[a], (s, b))
            push!(pos_contacts[b], (s, a))
        end
    end
    return pos_contacts
end

function init_weights(all_ctc_maps, seq)
    # weights[s] = exp(−E_s(seq)) pour chaque structure s
    return [exp(-energy_one_struct(ctc_map, seq)) for ctc_map in all_ctc_maps]
end


# ── Metropolis optimisé ────────────────────────────────────────────────────────

function metropolis(mut, seq, target_number, weights, pos_contacts, ttl_energy, log_proba_target, beta=1000)
    pos, new_aa = mut
    old_aa = seq[pos]

    # Calcul incrémental : seules les structures ayant un contact en `pos` changent
    affected = Dict{Int, Float64}()
    for (s, j) in pos_contacts[pos]
        affected[s] = get(affected, s, 0.0) + MJ[new_aa, seq[j]] - MJ[old_aa, seq[j]]
    end

    # Mise à jour de ttl_energy sans tout recalculer
    new_ttl_energy = ttl_energy
    for (s, dE) in affected
        new_ttl_energy += weights[s] * (exp(-dE) - 1.0)
    end

    # Nouvelle log-proba pour la structure cible uniquement (pas tout le vecteur !)
    dE_target = get(affected, target_number, 0.0)
    new_log_proba = log(weights[target_number] * exp(-dE_target)) - log(new_ttl_energy)

    # Critère Metropolis
    accept = new_log_proba >= log_proba_target ||
             rand() <= exp(beta * (new_log_proba - log_proba_target))

    if accept
        new_seq = copy(seq)            # copy seulement si on accepte
        new_seq[pos] = new_aa
        for (s, dE) in affected        # mise à jour du vecteur weights en place
            weights[s] *= exp(-dE)
        end
        return new_seq, new_ttl_energy, new_log_proba
    else
        return seq, ttl_energy, log_proba_target
    end
end


# ── Algo optimisé ─────────────────────────────────────────────────────────────

function algo(seq, target_number, all_ctc_maps, epochs=100, beta=100)
    # Précalculs (une seule fois)
    pos_contacts = precompute_pos_contacts(all_ctc_maps)
    weights      = init_weights(all_ctc_maps, seq)
    ttl_energy   = sum(weights)

    actual_seq        = copy(seq)
    log_proba_target  = log(weights[target_number]) - log(ttl_energy)
    p_struct_time_t   = Float64[]

    for i in 1:epochs
        if i % 100 == 0
            println("epoch number $i")
        end

        push!(p_struct_time_t, exp(log_proba_target))  # on travaille en log, exp à la fin

        mut = mutation(actual_seq)
        actual_seq, ttl_energy, log_proba_target =
            metropolis(mut, actual_seq, target_number, weights, pos_contacts, ttl_energy, log_proba_target, beta)
    end

    return actual_seq, p_struct_time_t
end


final_seq, p_struct_time_t = algo(seq, target_number, all_ctc_maps, 1000, 100)


using Plots

t = 1:length(p_struct_time_t)

plot(t, p_struct_time_t, xlabel="Time", ylabel="Probability", title="Evolution of the probability of the sequence over time")