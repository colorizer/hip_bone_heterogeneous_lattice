module FEA

using MKL
using LinearAlgebra
using SparseArrays
using Krylov

export fea

function local_stiffness(properties, l, A)
    E, G, ν = [properties[i] for i in [:E, :G, :ν]];
    k = 30/34;

    mI = A^2/(4π);
    J = 2mI;

    ϕ = 12E*mI/(l^2*G*k*A);
    μ = 1/(1+ϕ);
    μ2 = 3μ-1;
    μ4 = 3μ+1;

    a1 = E*A/l;
    a2 = 6E*mI*μ/l^2;
    a3 = 2a2/l;
    a4 = G*J/l;
    a5 = E*mI*μ4/l;
    a6 = E*mI*μ2/l;

    inds = [1 1; 7 7;
            1 7; 7 1;
            2 2; 3 3; 8 8; 9 9;
            2 8; 3 9; 8 2; 9 3;
            2 6; 2 12; 5 9; 9 11; 6 2; 12 2; 9 5; 11 9;
            3 5; 3 11; 5 3; 6 8; 8 6; 8 12; 11 3; 12 8;
            4 4; 10 10;
            4 10; 10 4;
            5 5; 6 6; 11 11; 12 12;
            5 11; 6 12; 11 5; 12 6;]
    vals = [a1,a1,
            -a1,-a1,
            a3,a3,a3,a3,
            -a3,-a3,-a3,-a3,
            a2,a2,a2,a2,a2,a2,a2,a2,
            -a2,-a2,-a2,-a2,-a2,-a2,-a2,-a2,
            a4,a4,
            -a4,-a4,
            a5,a5,a5,a5,
            a6,a6,a6,a6];

    K = sparse(inds[:,1], inds[:,2], vals, 12, 12);
    return K
end

function dK_dρₖ(properties, l, A, A₀, nₖ)
    E, G, ν= [properties[i] for i in [:E, :G, :ν]];
    k = 30/34;
    dAₖ_dρₖ = A₀/nₖ;

    inds=[1 1; 7 7; 1 7; 7 1;
          4 4; 10 10; 4 10; 10 4;
          2 6; 2 12; 5 9; 6 2; 9 5; 9 11; 11 9; 12 2;
          3 5; 3 11; 5 3; 6 8; 8 6; 8 12; 11 3; 12 8;
          5 11; 6 12; 12 6; 11 5;
          5 5; 6 6; 11 11; 12 12;
          2 2; 3 3; 8 8; 9 9;
          2 8; 3 9; 8 2; 9 3;];

    c = (3*(1+ν)/(π*k))/l^2;
    kGc²Al²=  π*G*k*c^2*A*l^2;
    A2cp1l = π*l*(2*A*c+1)^2;

    h0 = E/l;
    h1 = G*A/(π*l);
    h7 = 3*A*E + 2*kGc²Al²*A;
    h2 = 2*h7/((l^2)*A2cp1l);
    h4 = h2*l/2;
    h8 = 2*E*(A*c + 1. + (A*c)^2);
    h9 = 2*E*(A*c - 0.5 + (A*c)^2);
    h3 = A*(h8 + kGc²Al²)/A2cp1l;
    h5 = -A*(h9 - kGc²Al²)/A2cp1l;


    vals = vcat([h0, h0, -h0, -h0],
                [h1, h1, -h1, -h1],
                h4*ones(8), -h4*ones(8),
                h5*ones(4), h3*ones(4),
                h2*ones(4), -h2*ones(4)).*dAₖ_dρₖ;

    dK = sparse(inds[:,1], inds[:,2], vals, 12, 12);
    return dK
end

function make_global_indices!(rows, cols, nodes)
    r = rows.≤6;
    c = cols.≤6;
    rows[r] .+= (6*(nodes[1]-1));
    rows[.~r] .+= (6*(nodes[2]-2));
    cols[c] .+= (6*(nodes[1]-1));
    cols[.~c] .+= (6*(nodes[2]-2));
end


function get_K_el(mp, connection_id, l, A, connectivities, R)
    K_local = local_stiffness(mp, l, A);
    ks = sum(K_local);
    K_el = R' * K_local * R;
    connectivity_nodes = connectivities[connection_id, :];
    Krow, Kcol, Kdata = findnz(K_el);
    make_global_indices!(Krow, Kcol, connectivity_nodes);
    return (Krow, Kcol, Kdata)
end

function get_dK(mp, l, A, A₀, nₖ, R)
    dK_local = dK_dρₖ(mp, l, A, A₀, nₖ);
    dK_el = R' * dK_local * R;
    return dK_el
end

function assemble_stiffness(N_CONNECTIONS, N_NODES, connectivities, t_matrices, tangents, connectivities_phase_id, Lₑ, Aₑ, A₀, redundancy, mp)

    Kr = Vector{Vector{Int64}}(undef, N_CONNECTIONS);
    Kc = Vector{Vector{Int64}}(undef, N_CONNECTIONS);
    Kd = Vector{Vector{Float64}}(undef, N_CONNECTIONS);

    dK = Vector{SparseMatrixCSC{Float64, Int64}}(undef, N_CONNECTIONS);

    println("Generating local stiffness matrices")
    Threads.@threads for connection_id ∈ 1:N_CONNECTIONS
        R = t_matrices[tangents[connection_id]];
        l = Lₑ[connection_id];
        A = Aₑ[connection_id];
        nₖ = redundancy[connection_id];

        Kr[connection_id], Kc[connection_id], Kd[connection_id] = get_K_el(mp, connection_id, l, A, connectivities, R);
        dK[connection_id]= get_dK(mp, l, A, A₀[connectivities_phase_id[connection_id]], nₖ, R);
    end
    println("Assembling stiffness matrix")

    vcat2(i) = reduce(vcat, i);
    Kr = vcat2(Kr);
    Kc = vcat2(Kc);
    Kd = vcat2(Kd);

    dfn = 6;
    K = sparse(Kr, Kc, Kd, dfn*N_NODES, dfn*N_NODES);
    return K, dK
end

function fea(mp, N_NODES, connectivities_phase_id, connectivities, Lₑ, Aₑ, A₀, redundancy, t_matrices, tangents, F, active_indices)

    dfn = 6;
    N_CONNECTIONS = size(connectivities, 1);

    K, dK = assemble_stiffness(N_CONNECTIONS, N_NODES, connectivities, t_matrices, tangents, connectivities_phase_id, Lₑ, Aₑ, A₀, redundancy, mp);
    K = (K+K')/2;

    println("Trying to solve")
    U = zeros(6*N_NODES);
    println(size(K[active_indices, active_indices]));
    X =  cholesky(K[active_indices, active_indices]) \ F[active_indices];
    println("Solved FEA")
    U[active_indices] = X;
    Compliance = U'F;
    println("Compliance: $Compliance")
    return U, Compliance, dK

end

end
