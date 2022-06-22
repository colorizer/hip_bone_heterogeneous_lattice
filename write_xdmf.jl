include("./msh_frm_binvox.jl")
include("./gen_ltc_msh.jl")
import Main.gen_ltc_msh: generate_lattice
using PyCall
meshio = pyimport("meshio")

nodeCoords, node_ids, cell_id, connectivities, phase_cell_connectivity, phases, phase_counts, lengths, redundancy, boundary, connectivities_phase_id, Vᵤ = generate_lattice("hip_bone_hex.msh");

points = nodeCoords';
cells = [("line", connectivities.-1)];
mesh = meshio.Mesh(points, cells,
                    cell_data=PyDict(Dict("phase"=>[connectivities_phase_id], "redundancy"=>[redundancy])));
mesh.write("hip_bone_lattice.xdmf");
println("Generated xdmf file.")
