using JLD2
using GLMakie
using Colors
using DataFrames
μ1 = 0.355
μ2 = 0.109
combo(μ) = "$(μ)_0.4"

# https://levelup.gitconnected.com/makie-barplots-in-julia-2821dbe35366
constraints = combo(μ1);

nodeCoords, connectivities, connectivities_phase_id, d_e, l_e, A_e, df, c_df = load("hip_bone_lattice_outputs_"*constraints*".jld2", "nodeCoords", "connectivities", "connectivities_phase_id", "d_e", "L_e", "A_e", "df", "c_df");
println("loaded values")

phases = Dict(:VE=>1, :FE=>2, :VF=>3, :EE=>4, :BV=>5, :BF=>6, :BE=>7, :FF=>8);
clr_list = distinguishable_colors(8, [RGB(1,1,1), RGB(0,0,0)], dropseed=true);

clr_dict = Dict();

for phase in keys(phases)
    clr_dict[phase] = clr_list[phases[phase]];
end

V_e = A_e.*l_e;
V_phases = Dict();
V_t = sum(V_e);

for phase in keys(phases)
    V_phases[phase] = sum(V_e[connectivities_phase_id.==phases[phase]])/V_t*100;
end

function plot_composition()
    fig = Figure(resolution = (1024, 1024));
    xs = [values(phases)...];
    ys = [values(V_phases)...];
    clrs=[values(clr_dict)...];
    names = String.([keys(V_phases)...]);

    axis = Axis(fig[1,1],
                xticks = (xs, names),
                ylabel = "Volume %",
                title = "Phase composition μ=$μ1, κ=0.4")
    ylims!(axis, (0,50))
    barplot!(axis, xs, ys, color=clrs, bar_labels=:y, width=0.7)

    return fig
end

fig = plot_composition()
save("phase_comp_"*constraints*".png", fig)
display(fig)
