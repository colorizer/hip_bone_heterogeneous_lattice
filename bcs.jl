module bcs
export generate_neumann_bcs, generate_dirichlet_bcs, weights, neumann

using LinearAlgebra
using DataFrames
using CSV

weights = [0.04,0.11,0.22,0.13,0.04,0.11,0.22,0.13];

ltg =  [0.21113  0.36257   -0.90773; # local to global coordinates
        0.0      0.92866    0.37093;
        -0.97746  0.078313  -0.19607];

forces = Dict(
        "iliacus" => [0 0 0 -75.42 -101.55 -89.975 0 0;
                    0 0 0 -213.69 -287.73 -254.93 0 0;
                    0 0 0 25.14 33.851 29.992 0 0],
        "psoas" => [-49.288 0 -104.53 -57.888 -29.11 -57.888 -34.733 -46.311;
                    -139.65 0 -296.17 -164.02 -82.477 -164.02 -98.41 -131.21;
                    16.429 0 34.843 19.296 9.7032 19.296 11.578 15.437],
        "pectineus" => [0 0 -59.465 -32.621 0 -50.631 0 0;
                        0 0 -137.51 -75.436 0 -117.08 0 0;
                        0 0 90.437 49.611 0 77.001 0 0],
        "sartorius" => [0 -10.087 0 0 -4.0118 -18.11 -10.087 -10.087;
                        0 -85.993 0 0 -34.202 -154.4 -85.993 -85.993;
                        0 -15.728 0 0 -6.2556 -28.24 -15.728 -15.728],
        "rectus femoris" => [0 0 0 0 0 0 0 0;
                            0 -122.83 0 0 0 -174.75 -104.85 -95.865;
                            0 -6.5218 0 0 0 -9.279 -5.5674 -5.0902],
        "abductor longus" => [0 -15.939 0 0 -15.939 -28.617 -12.679 -25.357;
                            0 -76.594 0 0 -76.594 -137.52 -60.927 -121.85;
                            0 40.29 0 0 40.29 72.338 32.049 64.097],
        "abductor magnus 1" => [0 0 0 0 -11.498 -22.909 0 0;
                                0 0 0 0 -78.396 -156.2 0 0;
                                0 0 0 0 105.57 210.35 0 0],
        "abductor magnus 2" => [0 0 0 0 25.686 51.178 0 0;
                                0 0 0 0 -119.16 -237.41 0 0;
                                0 0 0 0 50.659 100.93 0 0],
        "abductor magnus 3" => [0 0 0 0 18.561 36.981 0 0;
                                0 0 0 0 -130.68 -260.38 0 0;
                                0 0 0 0 1.1364 2.2641 0 0],
        "abductor brevis" => [0 -19.047 0 0 0 -33.75 0 -19.047;
                                0 -71.219 0 0 0 -126.2 0 -71.219;
                                0 86.954 0 0 0 154.08 0 86.954],
        "gracillis" => [0 0 0 0 -5.4617 -9.8063 -4.3446 -8.6891;
                        0 0 0 0 -87.615 -157.31 -69.694 -139.39;
                        0 0 0 0 6.1444 11.032 4.8876 9.7752],
        "glutus maximus" => [314.51 347.38 62.379 140.82 170.33 183.4 42.582 180.04;
                            -689.5 -761.56 -136.75 -308.72 -373.41 -402.07 -93.353 -394.7;
                            366.93 405.28 72.775 164.29 198.72 213.97 49.679 210.05],
        "glutus medius 1" => [-336.53 -348.1 -487.27 -498.84 -466.78 -324.63 -34.711 -139.17;
                                -957.24 -990.15 -1386 -1418.9 -1327.7 -923.39 -98.733 -395.87;
                                82.263 85.091 119.11 121.94 114.1 79.354 8.4849 34.02],
        "glutus medius 2" => [-96.916 -100.25 -140.33 -143.66 -134.43 -93.489 -9.9963 -40.08;
                                -957.05 -989.95 -1385.7 -1418.7 -1327.5 -923.2 -98.713 -395.79;
                                333.15 344.6 482.38 493.83 462.09 321.37 34.362 137.78],
        "glutus medius 3" => [198.07 204.88 286.8 293.61 274.73 191.07 20.43 81.915;
                            -812.1 -840.02 -1175.9 -1203.8 -1126.4 -783.38 -83.763 -335.85;
                            581.02 600.99 841.28 861.25 805.89 560.47 59.928 240.28],
        "glutus minimus 1" => [-53.364 -32.768 -61.556 -53.364 -40.959 -28.789 -26.682 -51.258;
                                -213.46 -131.07 -246.22 -213.46 -163.84 -115.15 -106.73 -205.03;
                                59.768 36.7 68.943 59.768 45.875 32.243 29.884 57.409],
        "glutus minimus 2" => [14.562 8.9414 16.797 14.562 11.177 7.8557 7.2809 13.987;
                                -209.33 -128.53 -241.46 -209.33 -160.67 -112.93 -104.66 -201.06;
                                89.191 54.766 102.88 89.191 68.458 48.116 44.595 85.67],
        "glutus minimus 3" => [55.361 33.994 63.86 55.361 42.492 29.866 27.681 53.176;
                            -180.85 -111.05 -208.61 -180.85 -138.81 -97.562 -90.424 -173.71;
                            127.33 78.186 146.88 127.33 97.732 68.692 63.666 122.31],
        "tensor fascia lata" => [0 -5.8948 -3.9299 -7.0559 -6.654 -3.9299 -3.126 -4.2871;
                                0 -131.74 -87.824 -157.68 -148.7 -87.824 -69.86 -95.808;
                                0 -5.8948 -3.9299 -7.0559 -6.654 -3.9299 -3.126 -4.2871],
        "piriformis" => [111.47 151.76 0 0 0 0 67.877 125.82;
                        -81.071 -110.37 0 0 0 0 -49.365 -91.506;
                        147.67 201.03 0 0 0 0 89.915 166.67],
        "obturator internus" => [97.58 71.871 0 35.643 35.643 87.063 71.871 0;
                                12.457 9.175 0 4.5502 4.5502 11.114 9.175 0;
                                134.95 99.395 0 49.294 49.294 120.41 99.395 0],
        "glemellus superior" => [76.1 47.834 66.859 42.942 0 0 85.884 109.8;
                                -15.531 -9.7621 -13.645 -8.7637 0 0 -17.527 -22.409;
                                116.48 73.216 102.34 65.728 0 0 131.46 168.06],
        "glemellus inferior" => [0 0 0 0 0 84.848 47.878 90.302;
                                0 0 0 0 0 13.812 7.7941 14.7;
                                0 0 0 0 0 110.5 62.353 117.6],
        "quadratus femoris" => [6.8122 10.721 0 0 9.8274 20.548 0 0;
                                5.839 9.1893 0 0 8.4235 17.613 0 0;
                                60.337 94.956 0 0 87.043 182 0 0],
        "obturator externus" => [0 0 0 0 -30.719 -41.708 -32.966 -30.719;
                                0 0 0 0 22.448 30.479 24.091 22.448;
                                0 0 0 0 116.97 158.81 125.53 116.97],
        "semi tendinosus" => [0 10.874 8.1554 19.107 24.544 28.583 8.1554 0;
                                0 -139.26 -104.44 -244.69 -314.32 -366.04 -104.44 0;
                                0 -9.4708 -7.1031 -16.641 -21.377 -24.895 -7.1031 0],
        "semi membranosus" => [21.814 13.865 12.546 13.865 15.862 11.227 2.2982 15.862;
                                -577.35 -366.95 -332.05 -366.95 -419.8 -297.15 -60.826 -419.8;
                                -37.811 -24.032 -21.746 -24.032 -27.493 -19.461 -3.9836 -27.493],
        "hip joint force" => [0 0 0 0 0 0 0 0;
                            0 0 0 0 0 0 0 0;
                            -426 -2158 -1876 -1651 -1180 -187 -87 -379.0]
)

cell_ids_bc = DataFrame(CSV.File("cell_ids_bc.csv"));

get_cell_ids(muscle) = cell_ids_bc[!, muscle] |> skipmissing |> collect |> unique |> sort;

function get_connections(name, connection_ids, ret_cells=false)
    bc_cells = get_cell_ids(name);
    connections = hcat(connection_ids[bc_cells]...) |> unique |> sort;
    if ret_cells
        return bc_cells, connections
    else
        return connections
    end
end

struct neumann
    q::Array{Float64,2};
    connections;
    
    function neumann(name, Lₑ, cell_connection_ids, on_boundary)
        Fₗ = forces[name];
        connections = get_connections(name, cell_connection_ids);
        Lₜ = sum(Lₑ[connections]);
        return new(ltg*Fₗ./Lₜ, connections)
    end 
end

function generate_neumann_bcs(Lₑ, cell_connectivity, on_boundary)
    println("Generating Neumann boundary conditions")
    cc_df = DataFrame(cell_connectivity, ["cell_id", "connection_id"]);
    cc_gdf = combine(groupby(cc_df, :cell_id), :connection_id => Ref =>:connection_ids);
    connection_ids = map(x->hcat(x...), cc_gdf[:, :connection_ids]);
    bcs = [];
    for name in keys(forces)
        push!(bcs, neumann(name, Lₑ, connection_ids, on_boundary))
    end
    return bcs
end

function generate_dirichlet_bcs(cell_connectivity, connectivities, on_boundary)
    println("Generating Dirichlet boundary conditions")
    cc_df = DataFrame(cell_connectivity, ["cell_id", "connection_id"]);
    cc_gdf = combine(groupby(cc_df, :cell_id), :connection_id => Ref =>:connection_ids);
    connection_ids = map(x->hcat(x...), cc_gdf[:, :connection_ids]);
    names = ["fixed support"];
    dirichlet_bc = falses(size(on_boundary));
	bc_cells, connections = get_connections("fixed support", connection_ids, true);
	node_ids = reshape(connectivities[connections,:], :) |> unique |> sort;
	dirichlet_bc[node_ids] .= true;
    return bc_cells, dirichlet_bc .&& on_boundary
end

end
