import PolynomialRoots

param_dict = Dict(Symbol(k) => v for (k,v) ∈ ModelingToolkit.get_defaults(single_transcript_iFFL))

# RI = R / (κ1 + κ2 R)
# RIM = R / (κ3 + κ4 R)
pre(p) = p[:r_splicing] / p[:r_dicer] * (p[:regulated_copy] * p[:α_im]) / (p[:δ_im] + p[:r_splicing])
κ1(p) = (p[:k_miRNA_deg] * p[:δ_mi]) / (p[:k_miRNA_bind] * p[:r_dicer] * pre(p))
κ2(p) = (p[:k_miRNA_deg]) / (p[:r_dicer] * pre(p))
κ3(p) = (p[:δ_m] * p[:δ_mi] * p[:k_miRNA_deg]) * (p[:k_deg] + p[:k_mRNA_unbind]) / (p[:k_mRNA_bind] * p[:k_miRNA_bind] * p[:r_dicer] * pre(p))
κ4(p) = p[:δ_m] * p[:k_miRNA_deg] * (p[:k_deg] + p[:k_mRNA_unbind]) / (p[:k_mRNA_bind] * p[:r_dicer] * pre(p)) + p[:k_deg]
# Then
# So Rtot = R + RI + RIM
# Rtot (κ1 + κ2 R) (κ3 + κ4 R) = R (κ1 + κ2 R) (κ3 + κ4 R) + (κ3 + κ4 R) R + (κ1 + κ2 R) R
# 0 = κ2 κ4 R^3 + κ1 κ4 R^2 + κ2 κ3 R^2 + κ1 κ3 R + κ4 R^2 + κ3 R + κ2 R^2 + κ1 R - Rtot κ2 κ4 R^2 - Rtot κ1 κ4 R - Rtot κ2 κ3 R - Rtot κ1 κ3
# (κ2 κ4) R^3 + (κ1 κ4 + κ2 κ3 + κ4 + κ2 - Rtot κ2 κ4) R^2 + (κ1 κ3 + κ3  + κ1 - Rtot κ1 κ4 - Rtot κ2 κ3) R - Rtot κ1 κ3
function find_R_val(params)
    κ = [κ1(params), κ2(params), κ3(params), κ4(params)]
    Rtot = 170000
    roots = PolynomialRoots.roots([
        κ[2] * κ[4],
        κ[1] * κ[4] + κ[2] * κ[3] + κ[4] + κ[2] - Rtot * κ[2] * κ[4],
        κ[1] * κ[3] + κ[3] + κ[1] - Rtot * κ[1] * κ[4] - Rtot * κ[2] * κ[3],
        -Rtot * κ[1] * κ[3]
    ])
end

[κ1(param_dict), κ2(param_dict), κ3(param_dict), κ4(param_dict)]
find_R_val(param_dict)

"""
\text{RIM} = \frac{k_{mRNA\_b} k_{miR\_b} r_\text{dicer} \cdot \text{pre} \cdot \text{R}}{\delta_m \delta_{mi} k_{miR\_ub} (k_\text{deg} +

k_{mRNA\_ub}) + k_{miR\_b} (\delta_m k_\text{deg} k_{miR\_ub} + \delta_m k_{miR\_ub} k_{mRNA\_ub} + k_\text{deg} k_{mRNA\_m} r_\text{dicer} \cdot \text{pre}) R}
# NOTE: update the latex from k_{mRNA\_m} to k_{mRNA\_bind}
"""