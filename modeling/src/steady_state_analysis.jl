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
        "immature" => "imm",
        "unbind" => "ub",
        "bind" => "b",
        "RNA" => "R",
        "splicing" => "splice",
        "regulated" => "reg"
    ))
end