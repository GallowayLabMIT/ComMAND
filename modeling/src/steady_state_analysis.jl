using miR_iFFL
using Catalyst
using Latexify

# Make the output directory if it does not exist
outdir = "$(@__DIR__)/../../output/modeling/latex"
mkpath(outdir)

single_transcript_iFFL
open("$(outdir)/single_transcript_odes.tex", "w") do file
    ode_system_latex = latexify(convert(ODESystem, single_transcript_iFFL))
    write(file, replace(ode_system_latex,
       raw"\mathrm{pre}_{miRNA}\left( t \right)" => raw"\text{pre}(t)",
       raw"\mathrm{d} \mathrm{miRNA}\left( t \right)" => raw"\text{miR}(t)",
       raw"\mathrm{risc}\left( t \right)" => raw"\text{R}(t)",
       raw"\mathrm{risc}_{miRNA}\left( t \right)" => raw"\text{RI}(t)",
       raw"\mathrm{regulated}_{mRNA}\left( t \right)" => raw"\text{mRNA}(t)",
       raw"\mathrm{risc}_{miRNA\_mRNA}\left( t \right)" => raw"\text{RIM}(t)",
       raw"\mathrm{protein}\left( t \right)" => raw"P(t)",
       raw"\mathrm{immature}_{mRNA}\left( t \right)" => raw"\text{mRNA}_i(t)",
       raw"\mathrm{pri}_{miRNA}\left( t \right)" => raw"\text{pri}(t)",
       "unbind" => "ub",
       "bind" => "b",
    ))
end