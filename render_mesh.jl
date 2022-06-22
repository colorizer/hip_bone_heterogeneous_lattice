using JLD2
using GeometryBasics
using GLMakie
using Colors
using DataFrames

mu1 = 0.355
mu2 = 0.109
# based on https://juliaplots.org/MakieReferenceImages/gallery//merged_color_mesh/index.html
filename = "hip_bone_lattice_outputs_$(mu1)_0.4"

nodeCoords, connectivities, connectivities_phase_id, d_e, df = load(filename*".jld2", "nodeCoords", "connectivities", "connectivities_phase_id", "d_e", "df");
println("loaded values")

clr_list = distinguishable_colors(8, [RGB(1,1,1), RGB(0,0,0)], dropseed=true);
# opc = (d_e .- minimum(d_e))./(maximum(d_e)-minimum(d_e)); # Opacity from 0 to 1 based on diameter.
# opc = 0.5*((d_e .- minimum(d_e))./(maximum(d_e)-minimum(d_e))) .+ 0.5; # Opacity from 0.5 to 1 based on diameter.
opc = fill(0.8, length(d_e)); # constant opacity of 0.8

each_clr = clr_list[connectivities_phase_id];
clr_fnl = collect(zip(each_clr, opc));

r_e = d_e./2;
points = [];
for p in eachcol(nodeCoords)
    push!(points, Point(p...));
end
cylinders = Array{GeometryBasics.Mesh}(undef, size(d_e));

Threads.@threads for i in 1:length(d_e)
    start, stop = connectivities[i,:];
    obj = Cylinder(points[start], points[stop], r_e[i]);
    cylinders[i] = GeometryBasics.mesh(obj);
end
println("Cylinders ready")
fig = GLMakie.poly(cylinders, color=clr_fnl)
display(fig)
