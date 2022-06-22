module rho_update

using MKL
include("./FEA.jl")
include("./bcs.jl")

import Main.FEA: generate_t_matrices, fea
import Main.bcs: neumann
using LinearAlgebra
using SparseArrays
using DataFrames

function get_Area(ρᵢⱼ, A₀, cell_connectivity, connectivities_phase_id, redundancy, N_CONNECTIONS)
    Aₑ = zeros(N_CONNECTIONS);
    for (cell_id, connection_id) ∈ eachrow(cell_connectivity)
        phase_id = connectivities_phase_id[connection_id];
        Aₑ[connection_id] += ρᵢⱼ[cell_id, phase_id]*A₀[phase_id]/redundancy[connection_id];
    end
    return Aₑ
end

# for point forces. Not used.
# function apply_neumann_bc!(F, i, neumann_bc, node_ids)
#     Fx, Fy, Fz = neumann_bc.F[:,i];
#     bc_nodes = node_ids[neumann_bc.boundaries_boolean];
#     F[((bc_nodes.-1).*6).+1] .+= Fx;
#     F[((bc_nodes.-1).*6).+2] .+= Fy;
#     F[((bc_nodes.-1).*6).+3] .+= Fz;
# end

function apply_neumann_bc!(F, i, neumann_bc, connectivities, Lₑ)
    qx, qy, qz = neumann_bc.q[:,i];
    dfn = 6;
    for connection_id in neumann_bc.connections
        connectivity = connectivities[connection_id, :];
        inds = vcat([[(1:dfn).+(i-1)*dfn;] for i in connectivity]...);
        l = Lₑ[connection_id];
        pq = [qx*l/2, qy*l/2, qz*l/2, 0, -qz*l^2/12, qy*l^2/12, qx*l/2, qy*l/2, qz*l/2, 0, qz*l^2/12, -qy*l^2/12];
        F[inds] += pq;
    end
end


function get_iₜₕ_stiffness_sensitivity(i, mp, dfn, node_ids, connectivities_phase_id, connectivities, cell_connectivity, Lₑ, Aₑ, A₀, redundancy, t_matrices, tangents, neumann_bcs, active_indices, N_NODES, N_CELLS, N_CONNECTIONS, N_PHASES)
    F = zeros(dfn*N_NODES);
    println("Generating force")
    for neumann_bc ∈ neumann_bcs
        # apply_neumann_bc!(F, i, neumann_bc, node_ids)
        apply_neumann_bc!(F, i, neumann_bc, connectivities, Lₑ)
    end
    U, Compliance, dK = fea(mp, N_NODES, connectivities_phase_id, connectivities, Lₑ, Aₑ, A₀, redundancy, t_matrices, tangents, F, active_indices);

    Uₑ_∂Kₑ_Uₑ = zeros(N_CONNECTIONS);
    U_∂K_U = zeros(N_CELLS, N_PHASES);

    Threads.@threads for connection_id ∈ 1:N_CONNECTIONS
        connectivity = connectivities[connection_id, :];
        U_inds = vcat([[(1:dfn).+(i-1)*dfn;] for i in connectivity]...);
        U_e = U[U_inds];
        ∂K_e = dK[connection_id];
        Uₑ_∂Kₑ_Uₑ[connection_id] = U_e'*∂K_e*U_e;
    end

    println("Generating sensitivity");
    for (cell_id, connection_id) ∈ eachrow(cell_connectivity)
        phase_id = connectivities_phase_id[connection_id];
        U_∂K_U[cell_id, phase_id] += Uₑ_∂Kₑ_Uₑ[connection_id];
    end

    return Compliance, U_∂K_U, U

end

function get_C_U_∂K_U(mp, dfn, node_ids, connectivities_phase_id, connectivities, cell_connectivity, Lₑ, Aₑ, A₀, redundancy, t_matrices, tangents, neumann_bcs, active_indices, weights, N_NODES, N_CELLS, N_CONNECTIONS, N_PHASES)
    N_STAGES = 8;
    U_∂K_U_vals = Vector{Array{Float64, 2}}(undef, N_STAGES);
    Compliance = Vector{Float64}(undef, N_STAGES);
    U = Vector(undef, N_STAGES);
    for i ∈ 1:N_STAGES
        println("On Stage ", i)
        Compliance[i], U_∂K_U_vals[i], U[i] = get_iₜₕ_stiffness_sensitivity(i, mp, dfn, node_ids, connectivities_phase_id, connectivities, cell_connectivity, Lₑ, Aₑ, A₀, redundancy, t_matrices, tangents, neumann_bcs, active_indices, N_NODES, N_CELLS, N_CONNECTIONS, N_PHASES);
        println("Stage $i done \n \n")
    end
    U_∂K_U_val = sum(U_∂K_U_vals.*weights);
    U_val = sum(U.*weights);
    return Compliance, U_∂K_U_val, U_val
end

function update_densities(ρᵢⱼ, κ, β, ρ₀, ρ⁺, ρ⁻, U_∂K_U, ∂Vₖ_∂ρᵢⱼ, Lₑ, Vꜛ, Vᵤ, A₀, cell_connectivity, dirichlet_bc_cells, connectivities_phase_id, redundancy, N_CONNECTIONS, N_CELLS, N_PHASES)
    # ρᵢⱼ has the shape (N_CELLS, N_PHASES);
    sensitivity_ratio = (U_∂K_U./∂Vₖ_∂ρᵢⱼ).^β;

    V =sum(get_Area(ρᵢⱼ.*sensitivity_ratio, A₀, cell_connectivity, connectivities_phase_id, redundancy, N_CONNECTIONS).*Lₑ);

    λ = (V/Vꜛ)^(1/β);

    OC = sensitivity_ratio./(λ^β);

    ρᵢⱼ = ρᵢⱼ.*OC;

    for cell_id ∈ dirichlet_bc_cells
        ρᵢⱼ[cell_id, ρᵢⱼ[cell_id, :].<ρ₀] .= ρ₀;
    end

    ρ⁻_ind = ρᵢⱼ .<= ρ⁻;
    ρ⁺_ind = ρᵢⱼ .>= ρ⁺;

    ρᵢⱼ[ρ⁻_ind] .= ρ⁻;
    ρᵢⱼ[ρ⁺_ind] .= ρ⁺;

    ρ_bounds = ρ⁻_ind .| ρ⁺_ind;

    Vᵢₙ = ∂Vₖ_∂ρᵢⱼ.*ρᵢⱼ; # volume of individual phase in each cell. (shape N_CELLS×N_PHASES)

    vrholim = sum(Vᵢₙ[ρ_bounds]);
    Vst = Vꜛ - vrholim;

    ρᵢⱼ[.~ρ_bounds] .= Vst./sum(Vᵢₙ[.~ρ_bounds]).*ρᵢⱼ[.~ρ_bounds];

    # Now ρ satisfies the upper and lower bound as well as the volume constraint.

    Vᵢₙ = ∂Vₖ_∂ρᵢⱼ.*ρᵢⱼ;
    vcellea = sum(Vᵢₙ, dims=2)./Vᵤ;
    lgac = vec(vcellea .> κ);

    n_cell_grt_κ = count(lgac);

    facs = ones(length(lgac));
    facs[lgac] = κ./vcellea[lgac];
    ρᵢⱼ = ρᵢⱼ.*facs;
    Vactcell = κ*Vᵤ*n_cell_grt_κ;
    ρ⁻_ind = ρᵢⱼ .<= ρ⁻;
    ρ⁺_ind = ρᵢⱼ .>= ρ⁺;

    ρ⁻_ind[.~lgac, :] .= false;
    ρ⁺_ind[.~lgac, :] .= false;

    ρᵢⱼ[ρ⁻_ind] .= ρ⁻;
    ρᵢⱼ[ρ⁺_ind] .= ρ⁺;

    Vᵢₙ = ∂Vₖ_∂ρᵢⱼ.*ρᵢⱼ;

    ρ_bounds = ρ⁻_ind .| ρ⁺_ind;
    V_act_bounds = sum(Vᵢₙ[ρ_bounds]);

    lginact = .~ρ_bounds.&repeat(.~lgac, outer=(1,N_PHASES));
    facs = ones(N_CELLS, N_PHASES);
    facs[lginact] .= (Vꜛ-Vactcell-V_act_bounds)/sum(Vᵢₙ[lginact]);

    ρᵢⱼ .*= facs;

    return ρᵢⱼ;
end

end
