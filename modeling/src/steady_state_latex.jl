using miR_iFFL
using Catalyst
using Latexify

# Make the output directory if it does not exist
outdir = "$(@__DIR__)/../../output/modeling/latex"
mkpath(outdir)

variable_replacements = (
    raw"\mathrm{pre}_{miRNA}\left( t \right)" => raw"\text{pre}(t)",
       raw"\mathrm{d} \mathrm{miRNA}\left( t \right)" => raw"\text{miR}(t)",
       raw"\mathrm{risc}\left( t \right)" => raw"\text{R}(t)",
       raw"\mathrm{risc}_{miRNA}\left( t \right)" => raw"\text{RI}(t)",
       raw"\mathrm{regulated}_{mRNA}\left( t \right)" => raw"\text{mRNA}(t)",
       raw"\mathrm{risc}_{miRNA\_mRNA}\left( t \right)" => raw"\text{RIM}(t)",
       raw"\mathrm{protein}\left( t \right)" => raw"P(t)",
       raw"\mathrm{immature}_{mRNA}\left( t \right)" => raw"\text{mRNA}_i(t)",
       raw"\mathrm{pri}_{miRNA}\left( t \right)" => raw"\text{pri}(t)",
       raw"{dicer}" => raw"\text{dicer}",
       raw"{drosha}" => raw"\text{drosha}",
       raw"{splicing}" => raw"\text{splicing}",
       raw"regulated_{copy}" => raw"c_\text{regulated}",
       "miRNA" => "miR",
       "unbind" => "ub",
       "bind" => "b",
)

open("$(outdir)/single_transcript_odes.tex", "w") do file
    ode_system_latex = latexify(convert(ODESystem, single_transcript_iFFL))
    mapped = replace(ode_system_latex, variable_replacements...)
    reordered = getindex(split(mapped, "\n"), [1, 9, 10, 2, 3, 4, 5, 7, 6, 8, 11])
    write(file, join(reordered, "\n"))
end

open("$(outdir)/shared_odes.tex", "w") do file
    ode_system_latex = latexify(convert(ODESystem, single_transcript_iFFL))
    mapped = replace(ode_system_latex, variable_replacements...)
    reordered = getindex(split(mapped, "\n"), [1, 2, 3, 4, 5, 7, 6, 8, 11])
    write(file, join(reordered, "\n"))
end