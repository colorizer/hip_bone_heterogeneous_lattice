using MKL
include("./msh_frm_binvox.jl")
include("./gen_ltc_msh.jl")
include("./bcs.jl")
include("./FEA.jl")
include("./rho_update.jl")

using Plots
using LinearAlgebra
using SparseArrays
using JLD2
using DataFrames

import Main.gen_ltc_msh: generate_t_matrices, generate_lattice
import Main.FEA: fea
import Main.bcs: generate_neumann_bcs, generate_dirichlet_bcs, weights, neumann
import Main.rho_update: get_Area, get_C_U_∂K_U, update_densities

filename = "hip_bone_hex.msh";
nodeCoords, node_ids, cell_ids, connectivities, cell_connectivity, phases, phase_counts, lengths, redundancy, boundary_nodes_bool, connectivities_phase_id, Vᵤ = generate_lattice(filename);

t_matrices, tangents = generate_t_matrices(nodeCoords, connectivities); # tangents for each element and the corresponding transformation matrices. Used in FEA

# Material properties
mp = Dict(:E=>17e3, #MPa
            :ν=>0.3); #Poisson's ratio
mp[:G] = mp[:E]/(2*(1+mp[:ν]));

# Problem Parameters
ρ⁻ = 1e-8; # ρ min
ρ⁺ = 2; # ρ max

N_CELLS = size(cell_ids, 1);
N_CONNECTIONS = size(connectivities,1);
N_PHASES = length(phases);
N_NODES = length(node_ids);
dfn = 6;

β = 0.5; # parameter to control the speed of convergence.

rho_cortical = 1900 # kg/m^3
mass = 0.2904 # kg
mass_trabecular = 0.0892 # kg

Vₜ = Vᵤ*N_CELLS; # total volume within the boundary
case1 = mass*1e9/(Vₜ*rho_cortical);
case2 = mass_trabecular*1e9/(Vₜ*rho_cortical);

for μ in round.([case1, case2], digits=3)
    ρᵢⱼ = ones(N_CELLS, N_PHASES);
    κ = 0.4;
    Vꜛ = μ*Vₜ;

    println("""
    Total volume: $Vₜ
    V target: $Vꜛ
    volume ratio: $μ
    porosity constraint: $κ
    β: $β
    ρ max: $ρ⁺
    ρ min: $ρ⁻
    """)

    A₀ = zeros(N_PHASES);
    for phase in values(phases)
        d = lengths[phase]/3;
        A₀[phase] = (π*d^2)/4;
    end

    Lₑ = zeros(N_CONNECTIONS);
    for i in 1:N_CONNECTIONS
        Lₑ[i] = lengths[connectivities_phase_id[i]];
    end

    println("Generating ∑Lₖ/nₖ")
    ∑Lₖ_nₖ = zeros(N_CELLS, N_PHASES);
    for (cell_id, connection_id) ∈ eachrow(cell_connectivity)
        phase_id = connectivities_phase_id[connection_id];
        ∑Lₖ_nₖ[cell_id, phase_id] += Lₑ[connection_id]/redundancy[connection_id];
    end

    ∂Vₖ_∂ρᵢⱼ = ∑Lₖ_nₖ.*repeat(A₀', inner=(N_CELLS,1));

    get_area(ρᵢⱼ) = get_Area(ρᵢⱼ, A₀, cell_connectivity, connectivities_phase_id, redundancy, N_CONNECTIONS);

    Aₑ = get_area(ρᵢⱼ);
    V = sum(Aₑ.*Lₑ);
    println("Initial volume ", V)
    ρᵢⱼ .*= Vꜛ/V;
    ρ₀=maximum(ρᵢⱼ);
    println("Density val: ", ρ₀);
    Aₑ = get_area(ρᵢⱼ);
    V = sum(Aₑ.*Lₑ);
    V₀ = V;

    V_vals = Float64[];
    C_vals = Float64[];

    # GENERATE BOUNDARY CONDITIONS

    dirichlet_bc_cells, dirichlet_bc_nodes_bool = generate_dirichlet_bcs(cell_connectivity, connectivities, boundary_nodes_bool);
    println("Dirichlet bc generated. Fixed nodes count: ", count(dirichlet_bc_nodes_bool))
    non_dirichlet_node_ids = node_ids[.~dirichlet_bc_nodes_bool];
    active_indices = vcat([[(1:dfn).+(i-1)*dfn;] for i in non_dirichlet_node_ids]...);

    #**************** Generate Neumann boundaries ***************#
    # each one is a neumann_bc object
    neumann_bcs = generate_neumann_bcs(Lₑ, cell_connectivity, boundary_nodes_bool);
    #***********************************************#
    N_STAGES = 8;
    N_ITER = 30;
    mC_vals = Array{Float64,2}(undef,N_STAGES, N_ITER);


    for iter ∈ 1:N_ITER
        println("\n \n On iteration ",iter, "\n \n")
        global U;

        Aₑ = get_area(ρᵢⱼ);
        Compliances, U_∂K_U, U = get_C_U_∂K_U(mp, dfn, node_ids, connectivities_phase_id, connectivities, cell_connectivity, Lₑ, Aₑ, A₀, redundancy, t_matrices, tangents, neumann_bcs, active_indices, weights, N_NODES, N_CELLS, N_CONNECTIONS, N_PHASES);

        ρᵢⱼ = update_densities(ρᵢⱼ, κ, β, ρ₀, ρ⁺, ρ⁻, U_∂K_U, ∂Vₖ_∂ρᵢⱼ, Lₑ, Vꜛ, Vᵤ, A₀, cell_connectivity, dirichlet_bc_cells, connectivities_phase_id, redundancy, N_CONNECTIONS, N_CELLS, N_PHASES);

        mC_vals[:, iter] = Compliances;
        C = weights'Compliances;

        push!(C_vals, C);

        Vᵢₙ = ∂Vₖ_∂ρᵢⱼ.*ρᵢⱼ;
        V = sum(Vᵢₙ);
        push!(V_vals, V);

        println("Compliance: ", C)
        println("volume ", V)
        println("Max displacement", maximum(U));
        println("Minimum displacement", minimum(U));
    end

    println("Initial volume ", V₀)
    Compliances_df = DataFrame(mC_vals, :auto)
    println("Mean Compliance values for all stages")
    println(Compliances_df)
    println("Weighted mean Compliance")
    df = DataFrame(:C=>C_vals, :V=>V_vals)
    println(df)
    println("Max displacement ", maximum(U));
    println("Minimum displacement ", minimum(U));
    plot(Matrix(Compliances_df)', label=reshape(["Stage $i" for i in 1:8], (1,8)), lw=2, title="Mean Compliance", markershape = :auto)
    plot!(C_vals, lw=4, label="Weighted mC", title="Mean Compliance", markershape = :auto)
    savefig("Compliance_$(μ)_$(κ).png")
    plot(C_vals, lw=2, label="Weighted mC", title="Mean Compliance", markershape = :auto)
    savefig("Compliance_$(μ)_$(κ)_wmC.png")

    get_dia(area) = 2√(area/π);
    Aₑ = get_area(ρᵢⱼ);
    dₑ = get_dia.(Aₑ);

    save("hip_bone_lattice_outputs_$(μ)_$(κ).jld2", "nodeCoords", nodeCoords,
                                            "connectivities", connectivities,
                                            "connectivities_phase_id", connectivities_phase_id,
                                            "d_e", dₑ,
                                            "L_e", Lₑ,
                                            "A_e", Aₑ,
                                            "df", df, 
                                            "c_df", Compliances_df)

end
