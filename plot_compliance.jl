using JLD2
using GLMakie
using Colors
using DataFrames

μ1 = 0.355
μ2 = 0.109
combo(μ) = "$(μ)_0.4"

# https://levelup.gitconnected.com/makie-barplots-in-julia-2821dbe35366
constraints = combo(μ1);
nodeCoords, connectivities, connectivities_phase_id, d_e, l_e, A_e, df, c_df = load("final_results_beta_05/hip_bone_lattice_outputs_"*constraints*".jld2", "nodeCoords", "connectivities", "connectivities_phase_id", "d_e", "L_e", "A_e", "df", "c_df");
println("loaded values")

fig = Figure(resolution = (1024, 900));
xs = 1:30;
ys = df[!, :C];
axis = Axis(fig[1,1],
            xlabel="Iteration",
            ylabel="Compliance in Nmm",
            title = "Weighted mean Compliance μ=$μ, κ=0.4")
lines!(axis, 1:30, df[!, :C])
plot!(axis, 1:30, df[!,:C])
text!(axis, "$(round(df[1, :C], digits=2))", position=(1.5, df[1, :C]+1))
text!(axis, "$(round(df[end, :C], digits=2))", position=(28.5, df[end, :C]+15))
save("Compliance_"*constraints*".png", fig)
display(fig)
