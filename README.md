# gcat_post_analysis_script


## Results

### Somatic Mutation

* **Chr**: Mutation candidate positions
* **Pos**: Mutation candidate positions
* **Ref**: Reference base against the mutation candidate.
* **Alt**: Base sequences of the mutation candidate.
* **Read_Depth_Tumor**: Tumor read depth
* **Read_Depth_Forward_Strand_Tumor**: Tumor read depth in the forward strand
* **Read_Depth_Reverse_Strand_Tumor**: Tumor read depth in the reverse strand
* **Allelic_Depth_Tumor**: Tumor allelic depth
* **Allelic_Depth_Forward_Strand_Tumor**: Tumor allelic depth in the forward strand
* **Allelic_Depth_Reverse_Strand_Tumor**: Tumor allelic depth in the reverse strand
* **Read_Depth_Normal**: Normal read depth
* **Read_Depth_Forward_Strand_Normal**: Normal read depth in the forward strand
* **Read_Depth_Reverse_Strand_Normal**: Normal read depth in the reverse strand
* **Allelic_Depth_Normal**: Normal allelic depth
* **Allelic_Depth_Forward_Strand_Normal**: Normal allelic depth in the forward strand
* **Allelic_Depth_Reverse_Strand_Normal**: Normal allelic depth in the reverse strand
* **Mismatch_Rate_Tumor**: Tumor mismatch rate
* **Strand_Bias_Tumor**: Percentage of mismatched reads that are mapped along the + strand of Tumor
* **Mismatch_Rate_Normal**: Normal mismatch rate
* **Strand_Bias_Normal**: Percentage of mismatched reads that are mapped along the + strand of Normal
* **Realignment_Non-Allelic_Reads_Tumor**: Tumor number of non-allelic reads
* **Realignment_Allelic_Reads_Tumor**: Tumor number of allelic reads
* **Other_Reads_Tumor**: Tumor number of other reads
* **Realignment_Non-Allelic_Reads_Normal**: Normal number of non-allelic reads
* **Realignment_Allelic_Reads_Normal**: Normal number of allelic reads
* **Other_Reads_Normal**: Normal number of other reads
* **Number_of_Indel_Around_Mutation**: The number of indel-containing reads around ALT by the matched normal sample
* **Ratio_of_Indel_Around_Mutation**: The ratio of indel-containing reads around ALT by matched normal sample
* **Number_of_Breakpoint_Around_Mutation**: The number of breakpoint-containing reads around ALT by the matched normal sample
* **Length_From_Breakpoint_to_Mutation_Position**: The genome length from breakpoint-position to ALT position by matched normal sample
* **Fisher_Pval**: Minus logarithm of the p-value by Fishers exact test
* **Realignment_Fisher_Pval**: Minus logarithm of the p-value by Fishers exact test processed with Realignment Filter
* **Simple_Repeat_Region**: Simple repeat region
* **LOD_Score_of_Hotspot_Call**: LOD Score of Hotspot Call
* **EBCall_Pval**: EBCall Score

From here on down, the output of Emsemble VEP continues. For more information, please refer to the following website: https://m.ensembl.org/info/docs/tools/vep/vep_formats.html

* **Allele**: Emsemble VEP output field.
* **Consequence**: Emsemble VEP output field.
* **IMPACT**: Emsemble VEP output field.
* **SYMBOL**: Emsemble VEP output field.
* **Gene**: Emsemble VEP output field.
* **Feature_type**: Emsemble VEP output field.
* **Feature**: Emsemble VEP output field.
* **BIOTYPE**: Emsemble VEP output field.
* **EXON**: Emsemble VEP output field.
* **INTRON**: Emsemble VEP output field.
* **HGVSc**: Emsemble VEP output field.
* **HGVSp**: Emsemble VEP output field.
* **cDNA_position**: Emsemble VEP output field.
* **CDS_position**: Emsemble VEP output field.
* **Protein_position**: Emsemble VEP output field.
* **Amino_acids**: Emsemble VEP output field.
* **Codons**: Emsemble VEP output field.
* **Existing_variation**: Emsemble VEP output field.
* **ALLELE_NUM**: Emsemble VEP output field.
* **DISTANCE**: Emsemble VEP output field.
* **STRAND**: Emsemble VEP output field.
* **FLAGS**: Emsemble VEP output field.
* **MINIMISED**: Emsemble VEP output field.
* **SYMBOL_SOURCE**: Emsemble VEP output field.
* **HGNC_ID**: Emsemble VEP output field.
* **REFSEQ_MATCH**: Emsemble VEP output field.
* **REFSEQ_OFFSET**: Emsemble VEP output field.
* **GIVEN_REF**: Emsemble VEP output field.
* **USED_REF**: Emsemble VEP output field.
* **BAM_EDIT**: Emsemble VEP output field.
* **SOURCE**: Emsemble VEP output field.
* **CLIN_SIG**: Emsemble VEP output field.
* **SOMATIC**: Emsemble VEP output field.
* **PHENO**: Emsemble VEP output field.
* **VAR_SYNONYMS**: Emsemble VEP output field.
* **SpliceAI_pred_DP_AG**: SpliceAI predicted effect on splicing. Delta position for acceptor gain.
* **SpliceAI_pred_DP_AL**: SpliceAI predicted effect on splicing. Delta position for acceptor loss.
* **SpliceAI_pred_DP_DG**: SpliceAI predicted effect on splicing. Delta position for donor gain.
* **SpliceAI_pred_DP_DL**: SpliceAI predicted effect on splicing. Delta position for donor loss.
* **SpliceAI_pred_DS_AG**: SpliceAI predicted effect on splicing. Delta score for acceptor gain.
* **SpliceAI_pred_DS_AL**: SpliceAI predicted effect on splicing. Delta score for acceptor loss.
* **SpliceAI_pred_DS_DG**: SpliceAI predicted effect on splicing. Delta score for donor gain.
* **SpliceAI_pred_DS_DL**: SpliceAI predicted effect on splicing. Delta score for donor loss.
* **SpliceAI_pred_SYMBOL**: SpliceAI gene symbol.
* **LoF**: Loss-of-function annotation (HC = High Confidence; LC = Low Confidence).
* **LoF_filter**: Reason for LoF not being HC.
* **LoF_flags**: Possible warning flags for LoF.
* **LoF_info**: Info used for LoF annotation.
* **CADD_PHRED**: PHRED-like scaled CADD score.
* **CADD_RAW**: Raw CADD score.
* **gnomADg**: gnomad.genomes.r3.0.sites.vcf.bgz (exact).
* **gnomADg_AF**: AF field from gnomad.genomes.r3.0.sites.vcf.bgz.
* **gnomADg_AF_eas**: AF_eas field from gnomad.genomes.r3.0.sites.vcf.bgz".
* **ClinVar**: clinvar_20210102.vcf.gz (exact).
* **ClinVar_CLNSIG**: CLNSIG field from clinvar_20210102.vcf.gz.
* **ClinVar_CLNREVSTAT**: CLNREVSTAT field from clinvar_20210102.vcf.gz.
* **ClinVar_CLNDN**: CLNDN field from clinvar_20210102.vcf.gz.
* **ToMMo**: tommo-8.3kjpn-20200831-af_snvall_merged_liftedoverGRCh38.vcf.gz or tommo-8.3kjpn-20200831-af_indelall_merged_liftedoverGRCh38.vcf.gz (exact).
* **ToMMo_AN**: AN field from tommo-8.3kjpn-20200831-af_${snv/indel}all_merged_liftedoverGRCh38.vcf.gz.
* **ToMMo_AC**: AC field from tommo-8.3kjpn-20200831-af_${snv/indel}all_merged_liftedoverGRCh38.vcf.gz.
* **ToMMo_AF**: AF field from tommo-8.3kjpn-20200831-af_${snv/indel}all_merged_liftedoverGRCh38.vcf.gz.
* **Cancer_Gene_Census**: Role in Cancer and  Mutation Types field from cancer_gene_census_20210409.txt, separated by ":"
* **Target_Mutation_Flag**: 

### Somatic SV

* **Chr_1**: chromosome for the 1st breakpoint
* **Pos_1**: coordinate for the 1st breakpoint
* **Dir_1**: direction of the 1st breakpoint
* **Chr_2**: chromosome for the 2nd breakpoint
* **Pos_2**: coordinate for the 2nd breakpoint
* **Dir_2**: direction of the 2nd breakpoint
* **Inserted_Seq**: inserted nucleotides within the breakpoints
* **Variant_Type**: type of the structural variation
* **Gene_1**: gene overlapping the 1st breakpoint
* **Gene_2**: gene overlapping the 2nd breakpoint
* **Exon_1**: exon overlapping the 1st breakpoint
* **Exon_2**: exon overlapping the 2nd breakpoint
* **Num_Tumor_Ref_Read_Pair**: #read_pairs not supporting the variant (reference read pairs) for the tumor sample
* **Num_Tumor_Var_Read_Pair**: #read_pairs supporting the variant (variant read paris) for the tumor sample
* **Tumor_VAF**: frequency of variant read pairs for the tumor sample 
* **Num_Control_Ref_Read_Pair**: #read_pairs not supporting the variant for the matched control sample
* **Num_Control_Var_Read_Pair**: #read_pairs supporting the variant for the matched control sample
* **Control_VAF**: frequency of variant read pairs for the matched control sample 
* **Minus_Log_Fisher_P_value**: minus common logarithm of p-value for the Fisher's exact text (on contingency table of (tumor v.s. matched control) and (reference v.s. variant read pairs)
* **Non-Matched_Control_Sample_With_Max_Junction**: sample name with the maximum number of junction read pairs
* **Num_Max_Non-Matched_Control_Junction**: the maximum number of junction read pairs among non-matched control samples
* **Max_Over_Hang_1**: maximum overlang size of supporting read pairs from the 1st breakpoint
* **Max_Over_Hang_2**: maximum overlang size of supporting read pairs from the 2nd breakpoint
* **Role_in_Cancer_1**: Role in Cancer field from cancer_gene_census_20210409.txt from the 1st breakpoint
* **Mutation_Types_1**: Mutation Types field from cancer_gene_census_20210409.txt from the 1st breakpoint
* **Translocation_Partner_1**: Translocation Partnerfield from cancer_gene_census_20210409.txt from the 1st breakpoint
* **Role_in_Cancer_2**: Role in Cancer field from cancer_gene_census_20210409.txt from the 2nd breakpoint
* **Mutation_Types_2**: Mutation Types field from cancer_gene_census_20210409.txt from the 2nd breakpoint
* **Translocation_Partner_2**: Translocation Partner field from cancer_gene_census_20210409.txt from the 2nd breakpoint
* **Target_SV_Flag**: 

### Germline Mutation

* **CHROM**: Mutation candidate positions
* **POS**: Mutation candidate positions
* **REF**: Reference base against the mutation candidate.
* **ALT**: Base sequences of the mutation candidate. 
* **QUAL**: Quality field of VCF
* **FILTER**: Filter field of VCF
* **DP_INFO**: Approximate read depth; some reads may have been filtered
* **AC_INFO**: Allele count in genotypes, for each ALT allele, in the same order as listed
* **AN_INFO**: Total number of alleles in called genotypes
* **AF_INFO**: Allele Frequency, for each ALT allele, in the same order as listed
* **DS_INFO**: Were any of the samples downsampled?
* **END_INFO**: Stop position of the interval
* **ExcessHet_INFO**: Phred-scaled p-value for exact test of excess heterozygosity
* **FS_INFO**: Phred-scaled p-value using Fisher's exact test to detect strand bias
* **InbreedingCoeff_INFO**: Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation
* **MLEAC_INFO**: Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed
* **MLEAF_INFO**: Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed
* **MQ_INFO**: RMS Mapping Quality
* **BaseQRankSum_INFO**: Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities
* **MQRankSum_INFO**: Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities
* **ReadPosRankSum_INFO**: Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias
* **QD_INFO**: Variant Confidence/Quality by Depth
* **RAW_MQandDP_INFO**: Raw data (sum of squared MQ and total depth) for improved RMS Mapping Quality calculation. Incompatible with deprecated RAW_MQ formulation.
* **SOR_INFO**: Symmetric Odds Ratio of 2x2 contingency table to detect strand bias
* **GT_FORMAT**: Genotype
* **DP_FORMAT**: Approximate read depth (reads with MQ=255 or with bad mates are filtered)
* **AD_FORMAT**: Allelic depths for the ref and alt alleles in the order listed
* **GQ_FORMAT**: Genotype Quality
* **PGT_FORMAT**: Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another
* **PID_FORMAT**: Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group
* **PL_FORMAT**: Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification
* **PS_FORMAT**: Phasing set (typically the position of the first variant in the set)
* **RGQ_FORMAT**: Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)
* **SB_FORMAT**: Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.

From here on down, the output of Emsemble VEP continues. For more information, please refer to the following website: https://m.ensembl.org/info/docs/tools/vep/vep_formats.html

* **Allele**: Emsemble VEP output field.
* **Consequence**: Emsemble VEP output field.
* **IMPACT**: Emsemble VEP output field.
* **SYMBOL**: Emsemble VEP output field.
* **Gene**: Emsemble VEP output field.
* **Feature_type**: Emsemble VEP output field.
* **Feature**: Emsemble VEP output field.
* **BIOTYPE**: Emsemble VEP output field.
* **EXON**: Emsemble VEP output field.
* **INTRON**: Emsemble VEP output field.
* **HGVSc**: Emsemble VEP output field.
* **HGVSp**: Emsemble VEP output field.
* **cDNA_position**: Emsemble VEP output field.
* **CDS_position**: Emsemble VEP output field.
* **Protein_position**: Emsemble VEP output field.
* **Amino_acids**: Emsemble VEP output field.
* **Codons**: Emsemble VEP output field.
* **Existing_variation**: Emsemble VEP output field.
* **ALLELE_NUM**: Emsemble VEP output field.
* **DISTANCE**: Emsemble VEP output field.
* **STRAND**: Emsemble VEP output field.
* **FLAGS**: Emsemble VEP output field.
* **MINIMISED**: Emsemble VEP output field.
* **SYMBOL_SOURCE**: Emsemble VEP output field.
* **HGNC_ID**: Emsemble VEP output field.
* **REFSEQ_MATCH**: Emsemble VEP output field.
* **REFSEQ_OFFSET**: Emsemble VEP output field.
* **GIVEN_REF**: Emsemble VEP output field.
* **USED_REF**: Emsemble VEP output field.
* **BAM_EDIT**: Emsemble VEP output field.
* **SOURCE**: Emsemble VEP output field.
* **CLIN_SIG**: Emsemble VEP output field.
* **SOMATIC**: Emsemble VEP output field.
* **PHENO**: Emsemble VEP output field.
* **VAR_SYNONYMS**: Emsemble VEP output field.
* **SpliceAI_pred_DP_AG**: SpliceAI predicted effect on splicing. Delta position for acceptor gain
* **SpliceAI_pred_DP_AL**: SpliceAI predicted effect on splicing. Delta position for acceptor loss
* **SpliceAI_pred_DP_DG**: SpliceAI predicted effect on splicing. Delta position for donor gain
* **SpliceAI_pred_DP_DL**: SpliceAI predicted effect on splicing. Delta position for donor loss
* **SpliceAI_pred_DS_AG**: SpliceAI predicted effect on splicing. Delta score for acceptor gain
* **SpliceAI_pred_DS_AL**: SpliceAI predicted effect on splicing. Delta score for acceptor loss
* **SpliceAI_pred_DS_DG**: SpliceAI predicted effect on splicing. Delta score for donor gain
* **SpliceAI_pred_DS_DL**: SpliceAI predicted effect on splicing. Delta score for donor loss
* **SpliceAI_pred_SYMBOL**: SpliceAI gene symbol
* **LoF**: Loss-of-function annotation (HC = High Confidence; LC = Low Confidence)
* **LoF_filter**: Reason for LoF not being HC
* **LoF_flags**: Possible warning flags for LoF
* **LoF_info**: Info used for LoF annotation
* **CADD_PHRED**: PHRED-like scaled CADD score
* **CADD_RAW**: Raw CADD score
* **gnomADg**: gnomad.genomes.r3.0.sites.vcf.bgz (exact)
* **gnomADg_AF**: AF field from gnomad.genomes.r3.0.sites.vcf.bgz
* **gnomADg_AF_eas**: AF_eas field from gnomad.genomes.r3.0.sites.vcf.bgz"
* **ClinVar**: clinvar_20210102.vcf.gz (exact)
* **ClinVar_CLNSIG**: CLNSIG field from clinvar_20210102.vcf.gz
* **ClinVar_CLNREVSTAT**: CLNREVSTAT field from clinvar_20210102.vcf.gz
* **ClinVar_CLNDN**: CLNDN field from clinvar_20210102.vcf.gz
* **ToMMo**: tommo-8.3kjpn-20200831-af_snvall_merged_liftedoverGRCh38.vcf.gz or tommo-8.3kjpn-20200831-af_indelall_merged_liftedoverGRCh38.vcf.gz (exact)
* **ToMMo_AN**: AN field from tommo-8.3kjpn-20200831-af_${snv/indel}all_merged_liftedoverGRCh38.vcf.gz
* **ToMMo_AC**: AC field from tommo-8.3kjpn-20200831-af_${snv/indel}all_merged_liftedoverGRCh38.vcf.gz
* **ToMMo_AF**: AF field from tommo-8.3kjpn-20200831-af_${snv/indel}all_merged_liftedoverGRCh38.vcf.gz
* **Cancer_Gene_Census**: Role in Cancer and  Mutation Types field from cancer_gene_census_20210409.txt, separated by ":"
* **Target_Mutation_Flag**: 
