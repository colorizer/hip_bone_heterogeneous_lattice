using MKL
include("./msh_frm_binvox.jl")
include("./gen_ltc_msh.jl")
import Main.gen_ltc_msh: generate_lattice
using JLD2

filename = "hip_bone_hex"

nodeCoords, node_ids, cell_ids, connectivities, cell_connectivity, phases, phase_counts, lengths, redundancy, boundary_nodes_bool, connectivities_phase_id, Vᵤ = generate_lattice(filename*".msh");

save(filename*".jld2", "nodeCoords", nodeCoords,
                    "node_ids", node_ids,
                    "cell_ids", cell_ids,
                    "connectivities", connectivities,
                    "cell_connectivity", cell_connectivity,
                    "phases", phases,
                    "phase_counts", phase_counts,
                    "lengths", lengths,
                    "redundancy", redundancy,
                    "boundary_nodes_bool", boundary_nodes_bool,
                    "connectivities_phase_id", connectivities_phase_id,
                    "Vᵤ", Vᵤ)
