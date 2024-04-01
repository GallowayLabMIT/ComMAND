module miR_iFFL

# See SI for yangSyntheticCircuitBuffering2021 10.1038/s41467-021-23889-0 (https://www.nature.com/articles/s41467-021-23889-0)
# and carignanoExtrinsicNoiseSuppression2018   10.1109/CDC.2018.8619371   (https://ieeexplore.ieee.org/document/8619371)

using Catalyst

# α => production
# k => binding
# δ => mRNA degradation
# γ => protein degradation
nuclear_rn = @reaction_network NuclearReactions begin
    @parameters begin
        α_im=4.67e-2 #Transcription rate constant [1/s]
        δ_im=2.88e-4 #Immature mRNA degradation rate constant [1/s]
        δ_m=2.88e-4 #Mature mRNA degradation rate constant [1/s]
        δ_mi=2.88e-4 #miRNA degradation rate constant [1/s]
        δ_p=9.67e-5 #Protein degradation rate constant [1/s] 
        α_p=3.33e-4 #Translation rate constant [1/s]
        ζ=0.0 #Translation efficiency for RISC complex? [unitless]
        r_splicing=2.0e-3 #Splicing rate constant [1/s]
        r_drosha=1e-2 #pri-miRNA processing rate constant [1/s] 
        r_dicer=1e-3 #pre-miRNA processing rate constant [1/s]
        k_miRNA_bind=1e-5 #miRNA-RISC binding rate constant [1/mol*s]
        k_miRNA_unbind=2.16e-5 #miRNA-RISC unbinding rate constant [1/s]
        k_mRNA_bind=1.84e-6 #RISC complex formation rate constant[1/mol*s]
        k_mRNA_unbind=0.303 #RISC complex deformation rate constant [1/s]
        k_deg=7e-3 #RISC-mediated mRNA degradation rate constant [1/s]
    end
    # Transcription
    α_im, ∅ --> immature_mRNA # Should we include a copy number parameter?
    # Degradation of products
    δ_im, immature_mRNA --> ∅
    δ_m, mRNA --> ∅
    δ_mi, miRNA --> ∅
    δ_p, protein --> ∅
    # Production of protein
    α_p, mRNA --> protein
    α_p * ζ, risc_miRNA_mRNA --> protein # Is there an estimate for this parameter?
    # Splicing
    r_splicing, immature_mRNA --> mRNA + pri_miRNA
    # miRNA processing
    r_drosha, pri_miRNA --> pre_miRNA # drosha
    r_dicer, pre_miRNA --> miRNA # dicer
    # Knockdown
    (k_miRNA_bind, k_miRNA_unbind), risc + miRNA <--> risc_miRNA #The equalizer model has the unbinding going to just RISC, meaning the miRNA decays when it unbinds?
    (k_mRNA_bind, k_mRNA_unbind), risc_miRNA + mRNA <--> risc_miRNA_mRNA
    k_deg, risc_miRNA_mRNA --> risc_miRNA
end

end