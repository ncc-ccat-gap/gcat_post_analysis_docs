# gcat_post_analysis_script


## Somatic Mutation Results

* **Chr**: Mutation candidate positions
* **Pos**: Mutation candidate positions
* **Ref**: Reference base against the mutation candidate. "-" for insertions.
* **Alt**: Base sequences of the mutation candidate. "-" for deletions.
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
For more information about the Emsemble VEP output fields, please refer to the following website: https://m.ensembl.org/info/docs/tools/vep/vep_formats.html
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
* **SpliceAI_pred_DP_AG**: SpliceAI predicted effect on splicing. Delta position for acceptor gain. Customized emsemble VEP output field.
* **SpliceAI_pred_DP_AL**: SpliceAI predicted effect on splicing. Delta position for acceptor loss. Customized emsemble VEP output field.
* **SpliceAI_pred_DP_DG**: SpliceAI predicted effect on splicing. Delta position for donor gain. Customized emsemble VEP output field.
* **SpliceAI_pred_DP_DL**: SpliceAI predicted effect on splicing. Delta position for donor loss. Customized emsemble VEP output field.
* **SpliceAI_pred_DS_AG**: SpliceAI predicted effect on splicing. Delta score for acceptor gain. Customized emsemble VEP output field.
* **SpliceAI_pred_DS_AL**: SpliceAI predicted effect on splicing. Delta score for acceptor loss. Customized emsemble VEP output field.
* **SpliceAI_pred_DS_DG**: SpliceAI predicted effect on splicing. Delta score for donor gain. Customized emsemble VEP output field.
* **SpliceAI_pred_DS_DL**: SpliceAI predicted effect on splicing. Delta score for donor loss. Customized emsemble VEP output field.
* **SpliceAI_pred_SYMBOL**: SpliceAI gene symbol. Customized emsemble VEP output field.
* **LoF**: Loss-of-function annotation (HC = High Confidence; LC = Low Confidence). Customized emsemble VEP output field.
* **LoF_filter**: Reason for LoF not being HC. Customized emsemble VEP output field.
* **LoF_flags**: Possible warning flags for LoF. Customized emsemble VEP output field.
* **LoF_info**: Info used for LoF annotation. Customized emsemble VEP output field.
* **CADD_PHRED**: PHRED-like scaled CADD score. Customized emsemble VEP output field.
* **CADD_RAW**: Raw CADD score. Customized emsemble VEP output field.
* **gnomADg**: gnomad.genomes.r3.0.sites.vcf.bgz (exact). Customized emsemble VEP output field.
* **gnomADg_AF**: AF field from gnomad.genomes.r3.0.sites.vcf.bgz. Customized emsemble VEP output field.
* **gnomADg_AF_eas**: AF_eas field from gnomad.genomes.r3.0.sites.vcf.bgz". Customized emsemble VEP output field.
* **ClinVar**: clinvar_20210102.vcf.gz (exact). Customized emsemble VEP output field.
* **ClinVar_CLNSIG**: CLNSIG field from clinvar_20210102.vcf.gz. Customized emsemble VEP output field.
* **ClinVar_CLNREVSTAT**: CLNREVSTAT field from clinvar_20210102.vcf.gz. Customized emsemble VEP output field.
* **ClinVar_CLNDN**: CLNDN field from clinvar_20210102.vcf.gz. Customized emsemble VEP output field.
* **ToMMo**: tommo-8.3kjpn-20200831-af_snvall_merged_liftedoverGRCh38.vcf.gz or tommo-8.3kjpn-20200831-af_indelall_merged_liftedoverGRCh38.vcf.gz (exact). Customized emsemble VEP output field.
* **ToMMo_AN**: AN field from tommo-8.3kjpn-20200831-af_${snv/indel}all_merged_liftedoverGRCh38.vcf.gz. Customized emsemble VEP output field.
* **ToMMo_AC**: AC field from tommo-8.3kjpn-20200831-af_${snv/indel}all_merged_liftedoverGRCh38.vcf.gz. Customized emsemble VEP output field.
* **ToMMo_AF**: AF field from tommo-8.3kjpn-20200831-af_${snv/indel}all_merged_liftedoverGRCh38.vcf.gz. Customized emsemble VEP output field.
* **Cancer_Gene_Census**: 
* **Target_Mutation_Flag**: 
* 

## Somatic SV Results

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
