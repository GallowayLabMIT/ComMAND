using miR_iFFL
using Catalyst
import PolynomialRoots
using Makie
using CairoMakie
using Printf
using DataFrames
using Parquet2

# Make the output directory if it does not exist
outdir = "$(@__DIR__)/../../output/modeling/julia_param_sweeps"
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

function calculate_protein(Rtot, copy_number, ζ, k_mRNA_bind, k_mRNA_unbind)
    # Is this deepcopy inefficent? You betcha, but we'll use a more optimized
    # version of this when we go to fit
    params = deepcopy(param_dict)
    params[:Rtot] = Rtot
    params[:regulated_copy] = copy_number
    params[:ζ] = ζ
    params[:k_mRNA_bind] = k_mRNA_bind
    params[:k_mRNA_unbind] = k_mRNA_unbind
    R_vals = find_R_vals(params)
    @assert length(R_vals) == 1
    protein(R_vals[1], params)
end


logrange_around_param(p, decades, n) = 10.0 .^ range(log10(p) - (decades / 2.0), log10(p) + (decades / 2.0), n)

# Calculate initial sweeps to understand behavior
begin
f = Figure()
ax = Axis(f[1,1], xlabel="Copy number", ylabel="Steady state protein", title="Sweep of Rtot")
ylims!(ax, 400, 3000)
xlims!(ax, 0, 200)
for (idx, Rtot) ∈ enumerate([1e5, 1.7e5, 3e5, 1e6])
    lines!(ax, 1:200, calculate_protein.(
        Rtot,
        1:200,
        param_dict[:ζ],
        param_dict[:k_mRNA_bind],
        param_dict[:k_mRNA_unbind]
    ), label="Rtot=$(@sprintf("%.2e", Rtot))",
    color = idx, colormap = :viridis, colorrange = (0, 4.5))
end
hidespines!(ax, :r)
hidespines!(ax, :t)
axislegend(ax, position=:rb)
save("$outdir/rtot_sweep.pdf", f)
f
end

begin
f = Figure()
ax = Axis(f[1,1], xlabel="Copy number", ylabel="Steady state protein", title="Sweep of ζ")
ylims!(ax, 400, 7000)
xlims!(ax, 1, 250)
for (idx, ζ) ∈ enumerate(0:0.1:1.0)
    lines!(ax, 1:250, calculate_protein.(
        170000,
        1:250,
        ζ,
        param_dict[:k_mRNA_bind],
        param_dict[:k_mRNA_unbind]
    ), label="ζ=$(@sprintf("%.1f", ζ))",
    color = idx, colormap = :viridis, colorrange = (1, 12))
end
hidespines!(ax, :r)
hidespines!(ax, :t)
axislegend(ax, position=:rb)
save("$outdir/zeta_sweep.pdf", f)
save("$outdir/zeta_sweep.svg", f)
f
end

begin
f = Figure()
ax = Axis(f[1,1], xlabel="Copy number", ylabel="Steady state protein", title="Sweep of ζ", xscale=log10, yscale=log10)
ylims!(ax, 400, 10000)
xlims!(ax, 1, 500)
for (idx, ζ) ∈ enumerate(0:0.1:1.0)
    lines!(ax, 1:500, calculate_protein.(
        170000,
        1:500,
        ζ,
        param_dict[:k_mRNA_bind],
        param_dict[:k_mRNA_unbind]
    ), label="ζ=$(@sprintf("%.1f", ζ))",
    color = idx, colormap = :viridis, colorrange = (1, 12))
end
hidespines!(ax, :r)
hidespines!(ax, :t)
axislegend(ax, position=:rb)
save("$outdir/zeta_sweep_log.pdf", f)
save("$outdir/zeta_sweep_log.svg", f)
f
end

begin
f = Figure()
ax = Axis(f[1,1], xlabel="Copy number", ylabel="Steady state protein", title="Sweep of RISC-mRNA binding")
#ylims!(ax, 400, 3000)
xlims!(ax, 0, 200)
for (idx, k_mRNA_bind) ∈ enumerate(logrange_around_param(param_dict[:k_mRNA_bind], 1, 10))
    lines!(ax, 1:200, calculate_protein.(
        170000,
        1:200,
        param_dict[:ζ],
        k_mRNA_bind,
        param_dict[:k_mRNA_unbind]
    ), label="k_b=$(@sprintf("%.2e", k_mRNA_bind))",
    color = idx, colormap = :vikO, colorrange = (0, 12))
end
hidespines!(ax, :r)
hidespines!(ax, :t)
axislegend(ax, position=:rt)
save("$outdir/k_mRNA_bind_sweep.pdf", f)
f
end

begin
f = Figure()
ax = Axis(f[1,1], xlabel="Copy number", ylabel="Steady state protein", title="Sweep of RISC-mRNA unbinding")
#ylims!(ax, 400, 3000)
xlims!(ax, 0, 200)
for (idx, k_mRNA_unbind) ∈ enumerate(logrange_around_param(param_dict[:k_mRNA_unbind], 1, 10))
    lines!(ax, 1:200, calculate_protein.(
        170000,
        1:200,
        param_dict[:ζ],
        param_dict[:k_mRNA_bind],
        k_mRNA_unbind
    ), label="k_ub=$(@sprintf("%.2e", k_mRNA_unbind))",
    color = idx, colormap = :vikO, colorrange = (0, 12))
end
hidespines!(ax, :r)
hidespines!(ax, :t)
axislegend(ax, position=:rt)
save("$outdir/k_mRNA_unbind_sweep.pdf", f)
f
end

begin
f = Figure()
ax = Axis(f[1,1], xlabel="Copy number", ylabel="Steady state protein", title="Co-sweep of RISC-mRNA parameters")
ylims!(ax, 200, 10000)
xlims!(ax, 0, 200)
for (idx, (k_mRNA_bind,k_mRNA_unbind)) ∈ enumerate(zip(
    logrange_around_param(param_dict[:k_mRNA_bind], 1, 10),
    reverse(logrange_around_param(param_dict[:k_mRNA_unbind], 1, 10))))
    lines!(ax, 1:200, calculate_protein.(
        170000,
        1:200,
        param_dict[:ζ],
        k_mRNA_bind,
        k_mRNA_unbind
    ), label="k_b / k_ub=$(@sprintf("%.2e", k_mRNA_bind / k_mRNA_unbind))",
    color = idx, colormap = :vikO, colorrange = (0, 12))
end
hidespines!(ax, :r)
hidespines!(ax, :t)
axislegend(ax, position=:rt)
save("$outdir/k_mRNA_cosweep.pdf", f)
f
end

# Generate a whole dataframe of parameter values
N=4*500*11*10*10
Rtot_vals = zeros(N)
copy_n_vals = zeros(N)
ζ_vals = zeros(N)
bind_vals = zeros(N)
unbind_vals = zeros(N)
protein_vals = zeros(N)
idx = 1
for Rtot          = [1e5, 1.7e5, 3e5, 1e6],
    copy_n        = 1:500,
    ζ             = 0:0.1:1.0,
    k_mRNA_bind   = logrange_around_param(param_dict[:k_mRNA_bind], 1, 10),
    k_mRNA_unbind = logrange_around_param(param_dict[:k_mRNA_unbind], 1, 10)

    Rtot_vals[idx] = Rtot
    copy_n_vals[idx] = copy_n
    ζ_vals[idx] = ζ
    bind_vals[idx] = k_mRNA_bind
    unbind_vals[idx] = k_mRNA_unbind
    protein_vals[idx] = calculate_protein(Rtot, copy_n, ζ, k_mRNA_bind, k_mRNA_unbind)
    idx += 1
end

sweep_df = DataFrame(
    Rtot=Rtot_vals,
    copy_number=copy_n_vals,
    zeta=ζ_vals,
    k_mRNA_bind=bind_vals,
    k_mRNA_unbind=unbind_vals,
    protein=protein_vals
)

Parquet2.writefile("$outdir/sweep_df.gzip", sweep_df; compression_codec=:gzip)