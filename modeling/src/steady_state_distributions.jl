using miR_iFFL
using Catalyst
import PolynomialRoots
using Makie
using CairoMakie
using Printf
using DataFrames
using Parquet2
using Random, Distributions

# Make the output directory if it does not exist
outdir = "$(@__DIR__)/../../output/modeling/julia_distribution_sweeps"
mkpath(outdir)


param_dict = Dict(Symbol(k) => v for (k,v) ∈ ModelingToolkit.get_defaults(single_transcript_iFFL))

# RI = R / (κ1 + κ2 R)
# RIM = R / (κ3 + κ4 R)
pre(p) = p[:r_splicing] / p[:r_dicer] * (p[:regulated_copy] * p[:α_im]) / (p[:δ_im] + p[:r_splicing])
r_regulated(p) = p[:regulated_copy] * p[:α_im]
κ1(p) = (p[:k_miRNA_deg] * p[:δ_mi]) / (p[:k_miRNA_bind] * p[:r_dicer] * pre(p))
κ2(p) = (p[:k_miRNA_deg]) / (p[:r_dicer] * pre(p))
κ3(p) = (p[:δ_m] * p[:δ_mi] * p[:k_miRNA_deg]) * (p[:k_deg] + p[:k_mRNA_unbind]) / (p[:k_mRNA_bind] * p[:k_miRNA_bind] * p[:r_dicer] * pre(p) * r_regulated(p))
κ4(p) = p[:δ_m] * p[:k_miRNA_deg] * (p[:k_deg] + p[:k_mRNA_unbind]) / (p[:k_mRNA_bind] * p[:r_dicer] * pre(p) * r_regulated(p)) + (p[:k_deg] / r_regulated(p))
mRNA(R, p) = r_regulated(p) / (
        p[:δ_m] + 
            (p[:k_deg] * p[:k_mRNA_bind] * p[:k_miRNA_bind] * p[:r_dicer] * pre(p) * R)
            / (p[:k_miRNA_deg] * (p[:k_deg] + p[:k_mRNA_unbind]) * (p[:δ_mi] + p[:k_miRNA_bind] * R))
    )
protein(R,p) = p[:α_p] / p[:δ_p] * (mRNA(R,p) + p[:ζ] * R / (κ3(p) + κ4(p) * R))
# Then
# So Rtot = R + RI + RIM
# Rtot (κ1 + κ2 R) (κ3 + κ4 R) = R (κ1 + κ2 R) (κ3 + κ4 R) + (κ3 + κ4 R) R + (κ1 + κ2 R) R
# 0 = κ2 κ4 R^3 + κ1 κ4 R^2 + κ2 κ3 R^2 + κ1 κ3 R + κ4 R^2 + κ3 R + κ2 R^2 + κ1 R - Rtot κ2 κ4 R^2 - Rtot κ1 κ4 R - Rtot κ2 κ3 R - Rtot κ1 κ3
# (κ2 κ4) R^3 + (κ1 κ4 + κ2 κ3 + κ4 + κ2 - Rtot κ2 κ4) R^2 + (κ1 κ3 + κ3  + κ1 - Rtot κ1 κ4 - Rtot κ2 κ3) R - Rtot κ1 κ3
function find_R_vals(params)
    κ = [κ1(params), κ2(params), κ3(params), κ4(params)]
    Rtot = params[:Rtot]
    # Coeffs ordered from R^0 up to R^3
    coeffs = [
        -Rtot * κ[1] * κ[3],
        κ[1] * κ[3] + κ[3] + κ[1] - Rtot * κ[1] * κ[4] - Rtot * κ[2] * κ[3],
        κ[1] * κ[4] + κ[2] * κ[3] + κ[4] + κ[2] - Rtot * κ[2] * κ[4],
        κ[2] * κ[4]
    ]
    d, c, b, a = (coeffs[1], coeffs[2], coeffs[3], coeffs[4])
    discriminant = b * b * c * c - 4 * a * c * c * c - 4 * b * b * b * d - 27 * a * a * b * b + 18 * a * b * c * d
    @assert discriminant > 0.0
    roots = PolynomialRoots.roots(coeffs)
    # Return all positive real roots
    return map(real, filter(x->imag(x)==0.0 && real(x) > 0, roots))
end

function calculate_protein(copy_number, α_p, is_closed_loop)
    # Is this deepcopy inefficent? You betcha, but we'll use a more optimized
    # version of this when we go to fit
    params = deepcopy(param_dict)
    params[:regulated_copy] = copy_number
    params[:Rtot] = 1e5
    params[:α_p] = α_p
    if !is_closed_loop
        params[:k_deg] = 0
    end
    R_vals = find_R_vals(params)
    @assert length(R_vals) == 1
    protein(R_vals[1], params)
end

function calculate_slope(Rtot, is_closed_loop)
    params = deepcopy(param_dict)
    params[:regulated_copy] = 10000
    params[:Rtot] = Rtot
    if !is_closed_loop
        params[:k_deg] = 0
    end
    low_R_vals = find_R_vals(params)
    @assert length(low_R_vals) == 1
    low_protein = protein(low_R_vals[1], params)
    params[:regulated_copy] = 10100
    high_R_vals = find_R_vals(params)
    @assert length(high_R_vals) == 1
    high_protein = protein(high_R_vals[1], params)

    (high_protein - low_protein) / 100
end

logrange_around_param(p, decades, n) = 10.0 .^ range(log10(p) - (decades / 2.0), log10(p) + (decades / 2.0), n)


begin
f = Figure()
ax = Axis(f[1,1], xlabel="R_tot", xscale=log10, ylabel="Slope", title="Post-transition (RISC limited) slope")

r_range = logrange_around_param(1e5, 2, 50)
lines!(ax, r_range, calculate_slope.(r_range, false), label="Open loop")
lines!(ax, r_range, calculate_slope.(r_range, true), label="Closed loop")
hidespines!(ax, :r)
hidespines!(ax, :t)
axislegend(ax, position=:rc)
f
end

begin
f = Figure()
ax = Axis(f[1,1], xlabel="R_tot", xscale=log10, ylabel="Δ slope", title="Post-transition (RISC limited) Δ slope")

r_range = logrange_around_param(1e5, 2, 50)
lines!(ax, r_range, calculate_slope.(r_range, false) .- calculate_slope.(r_range, true))
hidespines!(ax, :r)
hidespines!(ax, :t)
f
end

# Sanity check that our closed loop/open loop predictions look correct
begin
f = Figure()
ax = Axis(f[1,1], xlabel="Copy number", ylabel="Steady state protein", title="Compare open/closed loop")
#ylims!(ax, 400, 3000)
xlims!(ax, 0, 100)
lines!(ax, 1:100, calculate_protein.(1:100, 0.000333, false), label="Open loop")
lines!(ax, 1:100, calculate_protein.(1:100, 0.000333, true), label="Closed loop")
hidespines!(ax, :r)
hidespines!(ax, :t)
axislegend(ax, position=:lt)
#save("$outdir/rtot_sweep.pdf", f)
f
end

# Plot over different distribution ranges
begin
f = Figure()
ax = Axis(f[1,1], xlabel="α_p", xscale=log10, ylabel="Steady state protein", title="Sweep α_p, within pre-linear")
xrange = logrange_around_param(0.000333, 1, 50)
lines!(ax, xrange, calculate_protein.(30, xrange, false), label="Open loop")
lines!(ax, xrange, calculate_protein.(30, xrange, true), label="Closed loop")
hidespines!(ax, :r)
hidespines!(ax, :t)
axislegend(ax, position=:lt)
f
end

begin
f = Figure()
ax = Axis(f[1,1], xlabel="α_p", xscale=log10, yscale=log10, ylabel="Steady state protein", title="Sweep α_p, within pre-linear")
xrange = logrange_around_param(0.000333, 1, 50)
lines!(ax, xrange, calculate_protein.(30, xrange, false), label="Open loop")
lines!(ax, xrange, calculate_protein.(30, xrange, true), label="Closed loop")
hidespines!(ax, :r)
hidespines!(ax, :t)
axislegend(ax, position=:lt)
f
end

begin
f = Figure()
ax = Axis(f[1,1], xlabel="α_p", xscale=log10, ylabel="Steady state protein", title="Sweep α_p, post-linear")
xrange = logrange_around_param(0.000333, 1, 50)
lines!(ax, xrange, calculate_protein.(150, xrange, false), label="Open loop")
lines!(ax, xrange, calculate_protein.(150, xrange, true), label="Closed loop")
hidespines!(ax, :r)
hidespines!(ax, :t)
axislegend(ax, position=:lt)
f
end

begin
f = Figure()
ax = Axis(f[1,1], xlabel="α_p", xscale=log10, yscale=log10, ylabel="Steady state protein", title="Sweep α_p, post-linear")
xrange = logrange_around_param(0.000333, 1, 50)
lines!(ax, xrange, calculate_protein.(150, xrange, false), label="Open loop")
lines!(ax, xrange, calculate_protein.(150, xrange, true), label="Closed loop")
hidespines!(ax, :r)
hidespines!(ax, :t)
axislegend(ax, position=:lt)
f
end

# Plot explicit output distributions, pulling from (zero-truncated) normal distributions
begin
f = Figure()
ax = Axis(f[1,1], xlabel="Steady state protein", xscale=log10, ylabel="Density", title="OL and CL distributions")
for (idx, d_width) ∈ enumerate([2e-5, 3e-5, 1e-5])
    distribution = truncated(Normal(0.000333, d_width), lower=1e-7)
    density!(ax, calculate_protein.(30, rand(distribution, 10000), true),
        color=(:black,0.0), strokecolor=("#008080", 0.8), strokewidth=3,
        label="Closed loop"
    )
    density!(ax, calculate_protein.(30, rand(distribution, 10000), false),
        color=(:black,0.0), strokecolor=("#808080", 0.8), strokewidth=3,
        label="Open loop"
    )
end
hidespines!(ax, :r)
hidespines!(ax, :t)
#axislegend(ax, position=:rb)
f
end

# Plot explicit output distributions, pulling from (zero-truncated) normal distributions
begin
f = Figure()
ax = Axis(f[1,1], xlabel="Steady state protein", xscale=log10, ylabel="Density", title="?OL and CL distributions")
d_width = 1e-5
for (idx, center) ∈ enumerate(logrange_around_param(0.000333, 1, 50))
    distribution = truncated(Normal(center, d_width), lower=1e-7)
    density!(ax, calculate_protein.(30, rand(distribution, 10000), true),
        color=(:black,0.0), strokecolor=("#008080", 0.8), strokewidth=3,
        label="Closed loop"
    )
    density!(ax, calculate_protein.(30, rand(distribution, 10000), false),
        color=(:black,0.0), strokecolor=("#808080", 0.8), strokewidth=3,
        label="Open loop"
    )
end
hidespines!(ax, :r)
hidespines!(ax, :t)
#axislegend(ax, position=:rb)
f
end

# Generate a whole dataframe of parameter values
N=3*10000*2
copy_number_vals = zeros(N)
is_closed_loop_vals = zeros(N)
distribution_width_values = zeros(N)
protein_vals = zeros(N)
idx = 1
for d_width        = [2e-5, 3e-5, 1e-5],
    is_closed_loop = [true, false]

    distribution = truncated(Normal(0.000333, d_width), lower=1e-7)

    for n = 1:10000
        copy_number = 30
        copy_number_vals[idx] = copy_number
        is_closed_loop_vals[idx] = is_closed_loop
        distribution_width_values[idx] = d_width
        protein_vals[idx] = calculate_protein(copy_number, rand(distribution), is_closed_loop)
    end
end

distribution_df = DataFrame(
    copy_number=copy_number_vals,
    is_closed_loop=is_closed_loop_vals,
    distribution_width=distribution_width_values,
    protein=protein_vals
)


Parquet2.writefile("$outdir/distribution_df.gzip", distribution_df; compression_codec=:gzip)

# Generate a dataframe for the slope dataframe
N=2*500
is_closed_loop_vals = zeros(N)
rtot_values = zeros(N)
slope_values = zeros(N)
idx = 1
for rtot           = logrange_around_param(1e5, 2, 500),
    is_closed_loop = [true, false]

    is_closed_loop_vals[idx] = is_closed_loop
    rtot_values[idx] = rtot
    slope_values[idx] = calculate_slope(rtot, is_closed_loop)
end

rtot_slope_df = DataFrame(
    is_closed_loop=is_closed_loop_vals,
    rtot=rtot_values,
    slope=slope_values
)


Parquet2.writefile("$outdir/slope_df.gzip", distribution_df; compression_codec=:gzip)