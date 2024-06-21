using miR_iFFL
using Catalyst
import PolynomialRoots
using Makie
using CairoMakie
using Printf
using DataFrames
using Parquet2

# Make the output directory if it does not exist
outdir = "$(@__DIR__)/../../output/modeling/steady_state_low_copy_number"
mkpath(outdir)
datadir = open(io->read(io, String),"$(@__DIR__)/../../datadir.txt")



struct SingleTranscript end
struct DualTranscript end
struct DualTranscriptU6 end

param_dict = Dict{Symbol,Any}(Symbol(k) => v for (k,v) ∈ ModelingToolkit.get_defaults(single_transcript_iFFL))
param_dict[:Rtot] = 1e5
param_dict[:transcript_type] = SingleTranscript()
param_dict[:ζ] = 0.7

# RI = R / (κ1 + κ2 R)
# RIM = R / (κ3 + κ4 R)
pre(p, _::SingleTranscript) = p[:r_splicing] / p[:r_dicer] * (p[:regulated_copy] * p[:α_im]) / (p[:δ_im] + p[:r_splicing])
pre(p, _::DualTranscript) = p[:r_splicing] / p[:r_dicer] * (p[:unregulated_copy] * p[:α_im]) / (p[:δ_im] + p[:r_splicing])
pre(p, _::DualTranscriptU6) = p[:u6_copy] * p[:α_u6] / p[:r_dicer]
pre(p) = pre(p, p[:transcript_type])
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

function calculate_protein(r_tot, transcript_type, copy_number)
    params = deepcopy(param_dict)
    params[:Rtot] = r_tot
    params[:transcript_type] = transcript_type
    params[:regulated_copy] = copy_number
    params[:unregulated_copy] = copy_number
    R_vals = find_R_vals(params)
    @assert length(R_vals) == 1
    protein(R_vals[1], params)
end

# Load stocahstic sim dataframes
stochastic_sims = DataFrame(Parquet2.readfile("$datadir/projects/miR-iFFL/modeling/julia_stochastic_simulations/stochastic_sims.gzip"))
unique(stochastic_sims.design)
# Compare "Design 1", "Dual Transcript", and "Dual Vector"
design_mapping = Dict(["Design 1" => SingleTranscript(), "Dual Transcript" => DualTranscript(), "Dual Vector" => DualTranscript()])
for (design, moi, R_tot) ∈ Iterators.product(
    ["Design 1", "Dual Transcript", "Dual Vector"],
    unique(stochastic_sims.moi),
    unique(stochastic_sims.risc),
)
    f = Figure()
    ax = Axis(f[1,1], title="$design: moi=$moi, RISC=$R_tot")
    plot_df = stochastic_sims[
        (stochastic_sims.risc .== R_tot) .&
        (stochastic_sims.design .== design) .&
        (stochastic_sims.moi .== moi), :]
    if size(plot_df, 1) == 0
        continue
    end
    scatter!(ax, plot_df.copynum, plot_df.reg_gene, color=(:blue, 0.05))
    scatter!(ax, plot_df.copynum, plot_df.unreg_gene, color=(:gray, 0.05))

    max_copy_number = maximum(plot_df.copynum)
    x_range = 1:max_copy_number
    ss_result = calculate_protein.(R_tot, [design_mapping[design]], x_range)
    lines!(ax, x_range, ss_result)
    display(f)
end
