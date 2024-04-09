module miR_iFFL

# See SI for yangSyntheticCircuitBuffering2021 10.1038/s41467-021-23889-0 (https://www.nature.com/articles/s41467-021-23889-0)
# and carignanoExtrinsicNoiseSuppression2018   10.1109/CDC.2018.8619371   (https://ieeexplore.ieee.org/document/8619371)

export single_transcript_iFFL, dual_transcript_iFFL, dual_transcript_U6_iFFL
using Catalyst
#using Latexify
#latexify(single_iFFL_rn)
#latexify(convert(ODESystem, single_iFFL_rn))

# RISC-mediated knockdown of the regulated gene, `regulated_mRNA`
dicer_risc_rn = @reaction_network dicer_RISC_complex begin
    @parameters begin
        r_dicer=1e-3            # pre-miRNA processing rate constant [1/s]
        δ_mi=2.88e-4            # miRNA degradation rate constant [1/s]
        k_miRNA_bind=1e-5       # miRNA-RISC binding rate constant [1/mol*s]
        k_miRNA_deg=2.16e-5  # miRNA-RISC unbinding rate constant [1/s]
        k_mRNA_bind=1.84e-6     # RISC complex formation rate constant[1/mol*s]
        k_mRNA_unbind=0.303     # RISC complex deformation rate constant [1/s]
        k_deg=7e-3              # RISC-mediated mRNA degradation rate constant [1/s]
    end
    # dicer
    r_dicer, pre_miRNA --> miRNA
    # degradation
    δ_mi, miRNA --> ∅
    # Knockdown
    k_miRNA_bind, risc + miRNA --> risc_miRNA
    k_miRNA_deg, risc_miRNA --> risc
    (k_mRNA_bind, k_mRNA_unbind), risc_miRNA + regulated_mRNA <--> risc_miRNA_mRNA
    k_deg, risc_miRNA_mRNA --> risc_miRNA
end

# Regulated gene dynamics
regulated_gene_rn = @reaction_network regulated_gene_dynamics begin
    @parameters begin
        δ_m=2.88e-4 # Mature mRNA degradation rate constant [1/s]
        δ_p=9.67e-5 # Protein degradation rate constant [1/s] 
        α_p=3.33e-4 # Translation rate constant [1/s]
        ζ=0.0       # Translation efficiency for RISC complex [unitless], Range: 0 to 1
    end
    # Production of protein
    α_p, regulated_mRNA --> protein
    α_p * ζ, risc_miRNA_mRNA --> protein + risc_miRNA_mRNA
    # Degradation of products
    δ_m, regulated_mRNA --> ∅
    δ_p, protein --> ∅
end

# Direct production of the miRNA from a U6 promoter
u6_rn = @reaction_network u6_miRNA_expression begin
    @parameters begin
        α_u6=4.67e-2  # Transcription rate constant [1/s]
        u6_copy=1     # Number of copies of U6
    end
    α_u6 * u6_copy, ∅ --> pre_miRNA
end

# Production of a regulated gene transcript, without an internal intron
nonintronic_regulated_gene_rn = @reaction_network gene_ts_expression begin
    @parameters begin
        α_m=4.67e-2     # Transcription rate constant [1/s]
        regulated_copy=1 # Number of copies of the regulated gene
    end
    α_m * regulated_copy, ∅ --> regulated_mRNA
end

# Production of an regulated transcript with its own miRNA intron
intronic_regulated_gene_rn = @reaction_network intronic_regulated_expression begin
    @parameters begin
        α_im=4.67e-2        # Transcription rate constant [1/s]
        δ_im=2.88e-4        # Immature mRNA degradation rate constant [1/s]
        r_splicing=2.0e-3   # Splicing rate constant [1/s]
        r_drosha=1e-2       # pri-miRNA processing rate constant [1/s] 
        regulated_copy=1    # Number of copies of the regulated gene
    end
    # Transcription
    α_im * regulated_copy, ∅ --> immature_mRNA
    # Degradation of products
    δ_im, immature_mRNA --> ∅
    # Splicing
    r_splicing, immature_mRNA --> regulated_mRNA + pri_miRNA
    # miRNA processing: drosha
    r_drosha, pri_miRNA --> pre_miRNA
end

# Production of an unregulated transcript with the miRNA intron
intronic_unregulated_gene_rn = @reaction_network intronic_unregulated_expression begin
    @parameters begin
        α_im=4.67e-2       # Transcription rate constant [1/s]
        α_p=3.33e-4        # Translation rate constant [1/s]
        δ_m=2.88e-4        # Mature mRNA degradation rate constant [1/s]
        δ_p=9.67e-5        # Protein degradation rate constant [1/s] 
        δ_im=2.88e-4       # Immature mRNA degradation rate constant gene 2 [1/s]
        r_splicing=2.0e-3   # Splicing rate constant [1/s]
        r_drosha=1e-2       # pri-miRNA processing rate constant [1/s] 
        unregulated_copy=1  # Number of copies of the unregulated gene
    end
    # Transcription
    α_im * unregulated_copy, ∅ --> immature_mRNA
    # Protein production
    α_p, unregulated_mRNA --> protein2
    # Degradation of products
    δ_im, immature_mRNA --> ∅
    δ_m, unregulated_mRNA --> ∅
    δ_p, unregulated_protein --> ∅
    # Splicing
    r_splicing, immature_mRNA --> unregulated_mRNA + pri_miRNA
    # miRNA processing: drosha
    r_drosha, pri_miRNA --> pre_miRNA
end

common_rn = extend(regulated_gene_rn, dicer_risc_rn)

@named single_transcript_iFFL = extend(intronic_regulated_gene_rn, common_rn)

@named dual_transcript_iFFL = extend(
    extend(nonintronic_regulated_gene_rn, intronic_unregulated_gene_rn),
    common_rn
)

@named dual_transcript_U6_iFFL = extend(
    extend(nonintronic_regulated_gene_rn, u6_rn),
    common_rn
)

single_iFFL_rn = @reaction_network single_transcript_iFFL begin
    @parameters begin
        α_im=4.67e-2 #Transcription rate constant [1/s]
        δ_im=2.88e-4 #Immature mRNA degradation rate constant [1/s]
        δ_m=2.88e-4 #Mature mRNA degradation rate constant [1/s]
        δ_mi=2.88e-4 #miRNA degradation rate constant [1/s]
        δ_p=9.67e-5 #Protein degradation rate constant [1/s] 
        α_p=3.33e-4 #Translation rate constant [1/s]
        ζ=0.0 #Translation efficiency for RISC complex [unitless], Range: 0 to 1
        r_splicing=2.0e-3 #Splicing rate constant [1/s]
        r_drosha=1e-2 #pri-miRNA processing rate constant [1/s] 
        r_dicer=1e-3 #pre-miRNA processing rate constant [1/s]
        k_miRNA_bind=1e-5 #miRNA-RISC binding rate constant [1/mol*s]
        k_miRNA_deg=2.16e-5 #miRNA-RISC unbinding rate constant [1/s]
        k_mRNA_bind=1.84e-6 #RISC complex formation rate constant[1/mol*s]
        k_mRNA_unbind=0.303 #RISC complex deformation rate constant [1/s]
        k_deg=7e-3 #RISC-mediated mRNA degradation rate constant [1/s]
    end
    # Transcription
    α_im, gene --> immature_mRNA + gene #Set copy number as the initial number of "gene"
    # Degradation of products
    δ_im, immature_mRNA --> ∅
    δ_m, mRNA --> ∅
    δ_mi, miRNA --> ∅
    δ_p, protein --> ∅
    # Production of protein
    α_p, mRNA --> protein + mRNA
    α_p * ζ, risc_miRNA_mRNA --> protein + risc_miRNA_mRNA
    # Splicing
    r_splicing, immature_mRNA --> mRNA + pri_miRNA
    # miRNA processing
    r_drosha, pri_miRNA --> pre_miRNA # drosha
    r_dicer, pre_miRNA --> miRNA # dicer
    # Knockdown
    k_miRNA_bind, risc + miRNA --> risc_miRNA
    k_miRNA_deg, risc_miRNA --> risc
    (k_mRNA_bind, k_mRNA_unbind), risc_miRNA + mRNA <--> risc_miRNA_mRNA
    k_deg, risc_miRNA_mRNA --> risc_miRNA
end

double_iFFL_rn = @reaction_network double_transcript_iFFL begin
    # gene 1 is regulated gene w miRNA binding single_transcript_iFFL
    # gene 2 is gene containing miRNA intron
    @parameters begin
        α_im1=4.67e-2 #Transcription rate constant gene 1 [1/s]
        α_im2=4.67e-2 #Transcription rate constant gene 2 [1/s]
        δ_im2=2.88e-4 #Immature mRNA degradation rate constant gene 2 [1/s]
        δ_m1=2.88e-4 #Mature mRNA degradation rate constant gene 1 [1/s]
        δ_m2=2.88e-4 #Mature mRNA degradation rate constant gene 2 [1/s]
        δ_mi=2.88e-4 #miRNA degradation rate constant [1/s]
        δ_p1=9.67e-5 #Protein degradation rate constant gene 1 [1/s] 
        δ_p2=9.67e-5 #Protein degradation rate constant gene 2 [1/s] 
        α_p1=3.33e-4 #Translation rate constant gene 1 [1/s]
        α_p2=3.33e-4 #Translation rate constant gene 2 [1/s]
        ζ=0.0 #Translation efficiency for RISC complex [unitless], Range: 0 to 1
        r_splicing=2.0e-3 #Splicing rate constant [1/s]
        r_drosha=1e-2 #pri-miRNA processing rate constant [1/s] 
        r_dicer=1e-3 #pre-miRNA processing rate constant [1/s]
        k_miRNA_bind=1e-5 #miRNA-RISC binding rate constant [1/mol*s]
        k_miRNA_deg=2.16e-5 #miRNA-RISC unbinding rate constant [1/s]
        k_mRNA_bind=1.84e-6 #RISC complex formation rate constant[1/mol*s]
        k_mRNA_unbind=0.303 #RISC complex deformation rate constant [1/s]
        k_deg=7e-3 #RISC-mediated mRNA degradation rate constant [1/s]
    end
    # Transcription
    α_im1, gene1 --> mRNA1 + gene1 #Set copy number as the initial number of "gene1"
    α_im2, gene2 --> immature_mRNA2 + gene2 #Set copy number as the initial number of "gene2"
    # Degradation of products
    δ_im2, immature_mRNA2 --> ∅
    δ_m1, mRNA1 --> ∅
    δ_m2, mRNA2 --> ∅
    δ_mi, miRNA --> ∅
    δ_p1, protein1 --> ∅
    δ_p2, protein2 --> ∅
    # Production of protein
    α_p1, mRNA1 --> protein1
    α_p2, mRNA2 --> protein2
    α_p1 * ζ, risc_miRNA_mRNA --> protein1
    # Splicing
    r_splicing, immature_mRNA2 --> mRNA2 + pri_miRNA
    # miRNA processing
    r_drosha, pri_miRNA --> pre_miRNA # drosha
    r_dicer, pre_miRNA --> miRNA # dicer
    # Knockdown
    k_miRNA_bind, risc + miRNA --> risc_miRNA
    k_miRNA_deg, risc_miRNA --> risc
    (k_mRNA_bind, k_mRNA_unbind), risc_miRNA + mRNA1 <--> risc_miRNA_mRNA
    k_deg, risc_miRNA_mRNA --> risc_miRNA
end

U6_iFFL_rn = @reaction_network double_transcript_U6_iFFL begin
    # gene 1 is regulated gene w miRNA binding single_transcript_iFFL
    # gene 2 is U6-driven miRNA expression
    @parameters begin
        α_im1=4.67e-2 #Transcription rate constant gene 1 [1/s]
        α_im2=4.67e-2 #Transcription rate constant gene 2 [1/s]
        δ_m=2.88e-4 #Mature mRNA degradation rate constant gene 1 [1/s]
        δ_mi=2.88e-4 #miRNA degradation rate constant [1/s]
        δ_p=9.67e-5 #Protein degradation rate constant gene 1 [1/s] 
        α_p=3.33e-4 #Translation rate constant gene 1 [1/s]
        ζ=0.0 #Translation efficiency for RISC complex [unitless], Range: 0 to 1
        r_dicer=1e-3 #pre-miRNA processing rate constant [1/s]
        k_miRNA_bind=1e-5 #miRNA-RISC binding rate constant [1/mol*s]
        k_miRNA_deg=2.16e-5 #miRNA-RISC unbinding rate constant [1/s]
        k_mRNA_bind=1.84e-6 #RISC complex formation rate constant[1/mol*s]
        k_mRNA_unbind=0.303 #RISC complex deformation rate constant [1/s]
        k_deg=7e-3 #RISC-mediated mRNA degradation rate constant [1/s]
    end
    # Transcription
    α_im1, gene1 --> mRNA + gene1 #Set copy number as the initial number of "gene1"
    α_im2, gene2 --> pre_miRNA + gene2 #Set copy number as the initial number of "gene2"
    # Degradation of products
    δ_m, mRNA --> ∅
    δ_mi, miRNA --> ∅
    δ_p, protein --> ∅
    # Production of protein
    α_p, mRNA --> protein
    α_p * ζ, risc_miRNA_mRNA --> protein
    # miRNA processing
    r_dicer, pre_miRNA --> miRNA # dicer
    # Knockdown
    k_miRNA_bind, risc + miRNA --> risc_miRNA
    k_miRNA_deg, risc_miRNA --> risc
    (k_mRNA_bind, k_mRNA_unbind), risc_miRNA + mRNA1 <--> risc_miRNA_mRNA
    k_deg, risc_miRNA_mRNA --> risc_miRNA
end

end