module miR_iFFL

# See SI for yangSyntheticCircuitBuffering2021 10.1038/s41467-021-23889-0 (https://www.nature.com/articles/s41467-021-23889-0)
# and carignanoExtrinsicNoiseSuppression2018   10.1109/CDC.2018.8619371   (https://ieeexplore.ieee.org/document/8619371)

using Catalyst

# α => production
# k => binding
# δ => mRNA degradation
# γ => protein degradation
nuclear_rn = @reaction_network NuclearReactions begin
    # Transcription
    α_im, ∅ --> immature_mRNA
    # Degradation of products
    δ_im, immature_mRNA --> ∅
    δ_m, mRNA --> ∅
    δ_mi, miRNA --> ∅
    δ_p, protein --> ∅
    # Production of protein
    α_p, mRNA --> protein
    α_p * ζ, risc_miRNA_mRNA --> protein
    # Splicing
    r_splicing, immature_mRNA --> mRNA + pri_miRNA
    # miRNA processing
    r_drosha, pri_miRNA --> pre_miRNA # drosha
    r_dicer, pre_miRNA --> miRNA # dicer
    # Knockdown
    (k_miRNA_bind, k_miRNA_unbind), risc + miRNA <--> risc_miRNA
    (k_mRNA_bind, k_mRNA_unbind), risc_miRNA + mRNA <--> risc_miRNA_mRNA
    k_deg, risc_miRNA_mRNA --> risc_miRNA
end

end
