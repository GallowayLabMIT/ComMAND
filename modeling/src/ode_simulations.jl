using miR_iFFL
using DifferentialEquations
using Makie
using CairoMakie

system = convert(ODESystem, single_transcript_iFFL)

problem = ODEProblem(single_transcript_iFFL, 
    [:pre_miRNA => 0, :miRNA => 0, :risc => 100, :risc_miRNA => 0, :regulated_mRNA => 0, :risc_miRNA_mRNA => 0, :protein => 0, :immature_mRNA => 0, :pri_miRNA => 0],
    (0., 50000.)
)

soln = solve(problem, Tsit5())
plot(soln)