using miR_iFFL
using Catalyst
import PolynomialRoots
using Makie
using CairoMakie
using Printf
using DataFrames
using Parquet2

# Make the output directory if it does not exist
outdir = "$(@__DIR__)/../../output/modeling/julia_param_sweeps/per_param"
mkpath(outdir)


param_dict = Dict(Symbol(k) => v for (k,v) ∈ ModelingToolkit.get_defaults(single_transcript_iFFL))
param_dict[:Rtot] = 1e5

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

function calculate_protein(param, value, copy_number)
    params = deepcopy(param_dict)
    params[param] = value
    params[:regulated_copy] = copy_number
    R_vals = find_R_vals(params)
    @assert length(R_vals) == 1
    protein(R_vals[1], params)
end


logrange_around_param(p, decades, n) = 10.0 .^ range(log10(p) - (decades / 2.0), log10(p) + (decades / 2.0), n)

sweep_df = DataFrame(
    param=String[],
    param_val=Float64[],
    copy_num=Int[],
    protein=Float64[]
)

# Calculate protein values & plot
for (key, val) ∈ param_dict
    if key ∈ [:regulated_copy]
        continue
    end

    if key == :ζ
        param_range = 0.0:0.1:1.0
    else
        param_range = logrange_around_param(val, 1, 11)
    end
    f = Figure()
    ax = Axis(f[1,1], xlabel="Copy number", ylabel="Steady state protein", title="Sweep of $(String(key))")
    for (idx, iterval) ∈ enumerate(param_range)
        copy_n = 1:100
        prot = calculate_protein.(key, iterval, copy_n)
        lines!(
            ax, copy_n, prot,
            label="$(@sprintf("%.1f", iterval / val))x",
            color = idx, colormap = :viridis, colorrange = (0, 12))

        df = DataFrame(
            param=String(key),
            param_val=iterval,
            copy_num=copy_n,
            protein=prot
        )
        append!(sweep_df,df)
    end

    hidespines!(ax, :r)
    hidespines!(ax, :t)
    axislegend(ax, position=:rb)
    save("$outdir/$(String(key)).pdf", f)
    save("$outdir/$(String(key)).svg", f)
end

Parquet2.writefile("$outdir/sweep_df.gzip", sweep_df; compression_codec=:gzip)

function calculate_slope(param, value, is_closed_loop)
    params = deepcopy(param_dict)
    params[param] = value
    params[:regulated_copy] = 100000
    if !is_closed_loop
        params[:k_deg] = 0
    end
    low_R_vals = find_R_vals(params)
    @assert length(low_R_vals) == 1
    low_protein = protein(low_R_vals[1], params)
    params[:regulated_copy] = 100100
    high_R_vals = find_R_vals(params)
    @assert length(high_R_vals) == 1
    high_protein = protein(high_R_vals[1], params)

    (high_protein - low_protein) / 100
end

slope_sweep_df = DataFrame(
    param=String[],
    param_val=Float64[],
    is_closed_loop=Bool[],
    slope=Float64[],
)

# Calculate protein values & plot
for (key, val) ∈ param_dict
    if key ∈ [:regulated_copy]
        continue
    end

    if key == :ζ
        param_range = 0.0:0.1:1.0
    else
        param_range = logrange_around_param(val, 1, 11)
    end
    f = Figure()
    ax = Axis(f[1,1], xlabel="$(String(key))", xscale=log10, ylabel="Slope", title="Sweep of $(String(key))")
    lines!(ax, param_range, calculate_slope.(key, param_range, false), label="Open loop")
    lines!(ax, param_range, calculate_slope.(key, param_range, true), label="Closed loop")

    df = DataFrame(
        param=String(key),
        param_val=param_range,
        is_closed_loop=false,
        slope=calculate_slope.(key, param_range, false)
    )
    append!(slope_sweep_df,df)
    df = DataFrame(
        param=String(key),
        param_val=param_range,
        is_closed_loop=true,
        slope=calculate_slope.(key, param_range, true)
    )
    append!(slope_sweep_df,df)

    hidespines!(ax, :r)
    hidespines!(ax, :t)
    axislegend(ax, position=:rb)
    save("$outdir/slope_$(String(key)).pdf", f)
    save("$outdir/slope_$(String(key)).svg", f)
end
outdir

Parquet2.writefile("$outdir/slope_sweep_df.gzip", sweep_df; compression_codec=:gzip)
