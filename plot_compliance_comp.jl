using JLD2
using GLMakie
using Colors
using DataFrames

μ1 = 0.355
μ2 = 0.109
combo(μ) = "$(μ)_0.4"
# https://levelup.gitconnected.com/makie-barplots-in-julia-2821dbe35366

c_df1, df1 = load("hip_bone_lattice_outputs_"*combo(μ1)*".jld2", "c_df", "df");
c_df2, df2 = load("hip_bone_lattice_outputs_"*combo(μ2)*".jld2", "c_df", "df");
println("loaded values")

sf = df2[1, :C]/df1[1,:C]
fig = Figure(resolution = (1024, 900));
xs = 1:30;
ys = [df1[!, :C] df2[!,:C] sf.*df1[!,:C]]';
axis = Axis(fig[1,1],
            xlabel="Iteration",
            ylabel="Compliance in Nmm",
            title = "Weighted mean Compliance")
series!(axis, xs, ys, markersize=8, labels=["Case 1", "Case 2", "What if"])
axislegend(axis)
save("Compliance_comp.png", fig)
display(fig)
