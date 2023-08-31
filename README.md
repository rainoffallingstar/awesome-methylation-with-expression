# awesome methylation with expression
Methylation :tw-1f493: Expression

### Paper list

1. (Integrated single-cell RNA-seq and DNA methylation reveal the effects of air pollution in patients with recurrent spontaneous abortion)[https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9400245/] (code)[] (data)[]
> We identified 98 DEGs with aberrant methylation by overlapping the RRBS and RNA-seq data. 
> A total of 4133 differentially methylated regions (DMRs) targeting 3526 differentially methylated genes (DMGs) were identified (P < 0.05 and |t|> 2.5), of which 1944 DMRs targeting 1645 DMGs showed hypermethylation, whereas 2187 DMRs targeting 1881 DMGs were hypomethylated in RSA compared to the controls
> A total of 928 differentially expressed genes (DEGs) were identified, including 380 downregulated DEGs and 548 upregulated DEGs.
> To explore how methylation affects expression, we used an integrated analysis strategy. After integrating methylomes and transcriptomes, ninety-eight DEGs with abnormal methylation were identified, comprising fifty-five hypomethyl-upregulated DEGs and forty-three hypermethyl-downregulated DEGs (Fig. 3A). According to the literature and enrichment analysis, six genes (ADAM12, FLT1, DLX3, IGF2BP1, F13A1, and FSTL3) were screened and verified by qRT-PCR. The qRT-PCR validation was performed in 16 controls and 16 patients with RSA. As shown in Fig. 3B, the expression levels of FLT1 and IGF2BP1 in the decidua were significantly higher in patients with RSA than in controls. Subsequently, a pyrosequencing assay was performed on two genes (FLT1 and IGF2BP1). We identified that DNA methylation levels of seven CpG sites in the IGF2BP1 promoter region in the decidua were significantly different between the patients with RSA and controls (Fig. 3C, Additional file 1: Table 3), while no CpG sites in the FLT1 promoter showed a statistical difference (Fig. 3D). Next, we performed receiver operating characteristic (ROC) curve analysis to examine the potential diagnostic value of methylation of IGF2BP1 in distinguishing patients with RSA and controls, indicating that the methylation level of IGF2BP1 has considerable potential as a diagnostic biomarker (Additional file 6: Fig. 4A). Interestingly, IGF2BP1 was also observed in 369 upregulated DEGs of scRNA-seq and was highly expressed in dM1 and dNK1 (Additional file 6: Fig. 4B). Collectively, these findings illustrated the abnormal methylation and expression level of IGF2BP1 may relate with RSA.
2. (Identification of Subtypes of Barrett’s Esophagus and Esophageal Adenocarcinoma Based on DNA Methylation Profiles and Integration of Transcriptome and Genome Data)[https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7305027/] (code)[] (data:Methylation data is accessible from European Genome-phenome Archive under accession numbers EGAD00010001822, EGAD00010001838 and EGAD00010001834.)[]
> We analyzed methylation profiles of all BE and EAC tissues and assigned them to subgroups using non-negative matrix factorization with k-means clustering. Data from whole-genome sequencing and transcriptome studies were then incorporated; we performed integrative methylation and RNA-seq analyses to identify genes that were suppressed with increased methylation in promoter regions.
> For assessing which genes undergo transcriptional repression under the influence of gaining methylation in promoter regions, we performed integrative methylation and RNA-seq analysis. For this analysis, we considered samples for which both RNA-seq and methylation were available. For each gene, we identified all probes located 1500 bp both up and downstream from the transcription start site (TSS). We selectively removed all CpG sites that were methylated in normal tissues (mean β-value >0.2). Methylation data was then dichotomised using β-value of ≥0.3 as a threshold (as used in TCGA studies13, 31) for positive DNA methylation, and discarded CpG sites methylated in fewer than 10% of samples. For each probe/gene pair, we then applied the following conditions: 1) categorized samples as either methylated (β ≥0.3) or unmethylated (β <0.3); 2) Compare expression in the methylated and unmethylated groups using the Mann-Whitney test; 3) Compute the correlation between methylation beta and expression TPM. We labelled each individual tumour sample as epigenetically silenced for a specific probe/gene pair selected above if for the probes there is a difference in beta (>0.2) between two groups, difference in distribution of expression of (adjusted p-value < 0.05) and negative correlation between methylation and expression (r < -0.1, adjusted p-value < 0.05). Only genes with multiple probes were considered for this analysis and a sample considered as epigenetically silenced if more than thirty percent of probes for the corresponding gene was also labelled as epigenetically silenced.
> We used ELMER32 for understanding which transcription factors are regulated upon perturbations from regulatory regions. Briefly, this method is based on initially identifying differentially methylated distal probes and predicting enriched motifs across them. Methylation levels from motif associated probes are then correlated with expression levels of transcription factor and ranked for any significant associations. We performed supervised analysis where each subtype was compared with others. On doing so we did not find significant results for most of the comparisons except for one, that between Subtype 2 and Subtype 3.
> We also observe that a few immune regulators (BLNK, CD40, VAV3, IRS2) are also affected by methylation.

3. (Genome-wide DNA methylation and RNA expression differences correlate with invasiveness in melanoma cell lines)[https://www.futuremedicine.com/doi/full/10.2217/epi-2020-0440?rfr_dat=cr_pub++0pubmed&url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org] (data:GSE153595)[]
> Integration of methylome & transcriptome data: mapping of ERV1 DMFs with DEGs & lncRNAs.To identify the effect of ERV1 hypomethylation on expression, the authors mapped the RRBS fragments overlapping with ERV1 to differentially expressed protein-coding genes, lncRNAs and pseudogenes. 
> To determine the association between gene expression and DNA methylation, DEGs (439 genes) were compared with DMFs (49 RRBS fragments). Overlaps based on the closest genes to the DMFs revealed six common genes with differential methylation and differential expression (Table 1). The overlapping DMFs were identified in the promoter region (two DMFs) as well as in the intronic region (six DMFs). 
> Validation of expression of selected genes using qRT-PCR

4. ***(Identification of epigenetic modulators in human breast cancer by integrated analysis of DNA methylation and RNA-Seq data)[https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6291302/] 
> flow chart of the computational pipeline for identifying modulator genes that are involved in dysregulated DNA methylation in breast cancer. In the first data preprocessing step, Illumina 450K methylation microarray intensities of both breast tumors and normal tissues are normalized; the batch effect is removed, and the effect of unknown and unmeasured variables is also removed with surrogate variable analysis (SVA). In the second step, DMRs in tumors relative to normal tissues are detected. In the third step, co-regulated DMRs are determined with hierarchical clustering, and the motifs of DNA binding proteins that are significantly enriched in each cluster of DMRs are identified. In the fourth step, gene expression levels of those proteins whose binding motifs are enriched in DMRs are correlated with the DNA methylation, and based on such correlation, methylation modulator genes for DMRs are determined. Finally, network analysis is performed to find network modules connected to modulator genes. More detailed description of each step is given in the following.
> (flowgraph)[https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6291302/figure/F0006/]

5. ***(An integrative analysis of DNA methylation and RNA-Seq data for human heart, kidney and liver)[https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3287572/]
> Next, we wanted to determine to what degree the gene expression differences among tissues are affected by epigenetic changes. We used cuffdiff to find expression variation (in FPKM) between tissues and then associated the selected genes with their methylation levels in HumanMethylation27 array. In total, we selected 1296 genes that were differentially expressed between any pair of the three tissues and that exist in HumanMethylatin27 array. We further looked at whether these 1296 genes showed DNA methylation variation between tissues as well. Unpaired t test was used to detect the methylation differences among tissues and the results are summarized in Figure ​Figure3.3
>  The distance of each CpG markers to the transcriptional start site (TSS) was used a covariate to fit into a linear model, but it was not identified as a confounding factor to influence gene expression. In short, the majority of differentially expressed genes showed significant changes in DNA methylation, implying DNA methylation plays an important role in mediating tissue differentiation.

6. (Identification of DNA methylation-regulated differentially expressed genes in RA by integrated analysis of DNA methylation and RNA-Seq data )[https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9588210/]

7. *("iNETgrate": integrating DNA methylation and gene expression data in a single gene network )[https://pubmed.ncbi.nlm.nih.gov/37645739/]

8. (Integrative analysis identifies two molecular and clinical subsets in Luminal B breast cancer)[https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10448479/]

9. (Integrative analysis of DNA methylome and transcriptome reveals epigenetic regulation of bisphenols-induced cardiomyocyte hypertrophy)[https://www.sciencedirect.com/science/article/pii/S0147651323008953?via%3Dihub]

### Methods, Materials & Resources
#### Data preprocess in expr

1. Sequencing data were aligned using STAR aligner

2. Using ENSEMBL gene annotation, counts of individual genes for all samples were computed using GenomicAlignments package from Bioconductor. Based on the counts, sequencing depth of individual samples and gene annotation, Transcripts Per Kilobase Million (TPM) for individual genes was computed across all samples. 

3. TPM were further corrected for batch effects using Combat

4. Tophat2 was used for alignment of RNA-seq data. A prebuilt human genome GRCh37 Bowtie2 index from Illumina's iGenomes collection was used as a reference genome.

5. Gencode 19 annotations were used to annotate all RNAs.

6. The read counts were generated using the htseq-count script, which is a part of the HTSeq Python library. 

7.  lncRNAs were annotated using the Gencode 19 annotation file.

8.  Fragments per kb million were generated using the Cufflinks pipeline.
9.  We used Tophat to map short reads onto the reference sequence of human genome 19 (GRCh37/hg19)
10. Software Cufflinks [23] was then used to quantify the expression levels of transcripts
11. Cutadapt software (https://cutadapt.readthedocs.io/en/stable/,version:cutadapt-1.9) was used to remove the reads containing adaptor contamination, (command line: ∼cutadapt -a ADAPT1 -A ADAPT2 -o out1.fastq -p out2.fastq in1.fastq in2.fastq -O 5 -m 100). And after removing the low quality bases and undetermined bases,
12. Cutadapt software (https://cutadapt.readthedocs.io/en/stable/,version:cutadapt-1.9) was used to remove the reads containing adaptor contamination, (command line: ∼cutadapt -a ADAPT1 -A ADAPT2 -o out1.fastq -p out2.fastq in1.fastq in2.fastq -O 5 -m 100). And after removing the low quality bases and undetermined bases



#### Data preprocess in methy

1. The sequenced reads were mapped against the complete human reference genome GRCh37/h19 using Bismark software (https://www.bioinformatics.babraham.ac.uk/projects/bismark/) 
2.  All raw data were processed using minfi
3.  Processed methylation data were further normalized using BETA mixture model BMIQ15 implemented in ChAMP package
4. Non-negative matrix factorization (NMF)18 on 5,000 most variable probes together with k-means clustering. Through NMF we first estimated optimal ranks/metagenes by executing it in combinations of 2–10 metagenes over 200 runs.
5. Differential analysis on individual probes was performed using linear models implemented in limma
6. (whole geneset)We used Strelka20 for calling somatic mutations, ASCAT21 for calling copy number and Manta22 for calling structural variants under similar settings as previously described
7. Based on the quality reports, adapter and hard trimming was carried out using the cleanadaptors program from the differential methylation analysis package
8. the effect of unknown and unmeasured variables is also removed with surrogate variable analysis (SVA),The data matrix output from ComBat and known variables including age and issue status were used by the SVA algorithm to generate values of surrogate variables, which were then used in the downstream analysis.
9. The raw data were preprocessed with software minifi  and normalized with SWAN to obtain the methylation percentages (-values) of CpG sites. All -values were transformed to M-values, as M-values yield more statistically robust results in finding differentially methylated CpG sites
10. The batch effect in the DNA methylation data was removed with computational method ComBat. More specifically, methylation data samples, the batch information (obtained from the ‘plate’ field in the TCGA bar code of each sample), and other covariates including age and tissue statues were input to the ComBat function. Additive and multiplicative batch parameters in the L/S model were estimated, and then the M-values were adjusted with these estimated parameters to correct the batch effect

11. First, we performed hierarchical clustering using all 27,578 CpG markers. We tested different dissimilarity measures, including Euclidean distance, Pearson's dissimilarity and Spearman's dissimilarity. Pearson's and Spearman's dissimilarity are both based on data correlation, thus, they have similar performance that is better than Euclidean distance. we intended to filter out the markers with low variance prior to clustering. We selected only the CpG markers that have a standard deviation greater than 0.2 and we obtained 488 markers. 

12. to further investigate methylation clustering, we used an alternative method, principal components analysis (PCA). One of the advantages of PCA is that it disassociates the correlation between markers when methylation data are transformed into principal components. We plotted the 18 samples in the first three principal components

13. Cutadapt (Valenzuela et al., 2011) and Perl scripts in house were utilized to cut out the reads that contained inaccurate bases. 

#### QC methods in expr

1. Quality of the raw reads was assessed using the FastQC tool. 

2. Trimmomatic was used to discard low-quality reads

3.  we drew quantile-quantile plots to compare the expression distributions between tissues

#### QC methods in methy
1. FastQC application (http://www.bioinformatics.babraham.ac.uk/projects/fastqc) to evaluate the quality of sequencing data


#### DEGs analysis in expr
1. RNA-seq data were normalized and analyzed using the limma package to identify DEGs.
2. converted the obtained matrix into a Seurat object using the Seurat package (Version 4.0.2) and integrated the single-cell data using the Harmony package. Finally, the RunPCA and RunTSNE functions were used to linearly scale the expression data and nonlinearly scale dimensionality reduction.
3. Differential analysis of each individual subtype over all other subtypes was performed on counts using the edgeRpackage.
4. the DESeq2 R package was used for the differential expression analysis according to the standard protocol. 
5. DEGs were performed using the DESeq R package. Adjust P-value (False Discovery Rate) < 0.05 and |log2FoldChange| ≥1.5 were set as the threshold for significant differential expression. 



#### MDRs,MDGs analysis in methy

1. Data processing and analyses were performed using methylKit package, a package is used to descript the methylation statistics of all samples and find differentially methylated regions from sorted Bismark-aligned BAM files in software R (Version 4.1.0).
>  In methylKit, tileMethylCounts function was used to discover de novo differentially methylated genes (DMRs). DMRs were defined as a sliding window with a default window size of 1000 bp. Specifically, statistical tests for differential methylation at each region were performed by function “calculateDiffMeth” with default parameters, which output was then processed using function “getMethylDiff” to call DMRs by comparing the RSA group with the control group.
2.  ChIPseeker was used to annotate the differentially methylated regions.
3. Stockwell PA, Chatterjee A, Rodger EJ, Morison IM. DMAP: differential methylation analysis package for RRBS and WGBS data. Bioinformatics 30(13), 1814–1822 
4. Based on the co-methylated regions determined from the a-clustering algorithm and the ‘s estimated from model (1) with Limma, the bump hunting algorithm  was employed to find DMRs. Specifically, each co-methylated region was assigned a new test statistic , where stands for the th co-methylated region. The null distribution of this statistic was determined by random permutation of the disease status of samples, and DMRs were found at a false discovery rate (FDR) of , which yielded a total of DMRs. The test statistic was also used to identify a set of non-differentially methylated regions

5. Differentially methylated regions (DMRs) were estimated by the R package-MethylKit (Akalin et al., 2012) with default parameters (1000 bp slide windows, 500 bp overlap, p value< 0.05).

#### The intergration of MDRs,MDGs, and DEGs

1. co-regulated DMRs are determined with hierarchical clustering, and the motifs of DNA binding proteins that are significantly enriched in each cluster of DMRs are identified. 
2. gene expression levels of those proteins whose binding motifs are enriched in DMRs are correlated with the DNA methylation, and based on such correlation, methylation modulator genes for DMRs are determined.
3. For the DNA binding proteins that we compiled from 7 databases, a network was constructed with FIMO. Specifically, the DNA sequence of 900 bp, starting at 700 bp before the transcription start site (TSS) and ending at 200 bp after the TSS, of each gene encoding a DNA binding protein was extracted from the human genome. For each of 964 genes, FIMO was employed to search over the 900 bp long DNA sequences of all other 963 genes to determine if the DNA motif of the gene is present. 
4. used cuffdiff to find expression variation (in FPKM) between tissues and then associated the selected genes with their methylation levels in HumanMethylation27 array.(what degree the gene expression differences among tissues are affected by epigenetic changes)
5. The distance of each CpG markers to the transcriptional start site (TSS) was used a covariate to fit into a linear model, but it was not identified as a confounding factor to influence gene expression.  In short, the majority of differentially expressed genes showed significant changes in DNA methylation, implying DNA methylation plays an important role in mediating tissue differentiation.
6. Venn diagrams of genes and pathways were drawn using OmicStudio tools (details shown in www.omicstudio.cn).

#### The commom machine learning workflows

1. we performed univariate Cox regression analysis and multivariate Cox regression analysis using the “Survival” R package (Therneau 2023). The nodes of the network were included in the initial Cox model. Subsequently, the “glmnet” package in R was used to further construct the model through least absolute shrinkage and selection operator (LASSO) regression (Friedman et al. 2010). The genes with the lowest cross-validation points were selected to prevent overfitting. Then forest plot and nomogram were plotted. To evaluate the discrimination and accuracy, the receiver operating characteristic (ROC) curve was plotted, and the areas under the curve (AUC) of patients after 1, 3 and 5 years were calculated (Kamarudin et al. 2017). The calibration curve was also plotted. Finally, we calculated the survival risk score of each patient and divided them into high-risk and low-risk groups according to the median for Kaplan–Meier (K–M) survival curve analysis. Gene prognostic factors were identified by Cox risk model.

#### The deeper layer of enrichment analysis, network anslysis,and model construction

1. The Database for Annotation, Visualization, and Integrated Discovery (DAVID, https://david.ncifcrf.gov/home.jsp) was used to perform Gene Ontology (GO) annotations and Kyoto Encyclopedia of Genes and Genomes (KEGG) pathways analysis.
2.  Pathway analysis was performed on ranked data from differential analysis using Gene Set Enrichment Analysis (GSEA)
3. [ppi network] The STRING database was used for protein–protein interactions. An extended string network was created with 557 differentially expressed genes (DEGs).
4. [co-express network] To determine the potential functions of the lncRNAs, a protein-coding gene–lncRNA co-expression network was constructed. 
5. network analysis is performed to find network modules connected to modulator genes
6.  A network of these 964 genes was inferred from their expression levels using the ACRANE algorithm with a p value threshold of for the mutual information and a tolerance equal to . The network created by FIMO (named FIMO network) was modified with the network created by the ACRANE (named ACRANE network). Specifically, an edge between two genes in the FIMO network was removed if no edge existed between the same two genes in the ACRANE network.
7. The analysis was conducted using the clusterProfiler package of R software, and the data set was from the Molecular Signatures Database v7.2 (MSigDB) downloaded from the GSEA-MSigDB website. The MSigDB is a database of gene sets for performing gene set enrichment analysis
8. Weighted gene co-expression network analysis (WGCNA) is used to construct the gene coexpression network and identify the functional modules [10]. The WGCNA was down by the WGCNA R software package
9. The protein–protein interaction (PPI) network of MeDEGs was constructed by the STRING database.  The PPI network was displayed by Reactome functional interactions (FI) Cytoscape Plugin. Reactome FI Cytoscape Plugin can be used to find network patterns by the Reactome database 
10. Correlation analysis between methylation level and mRNA expression,
R software was used to explore the association between methylation level and expression by spearman correlation analysis. R ≥ 0.4, P < 0.05 was used as the criterion for screening correlation. 
11.  The “estimate” package of R language (Yoshihara et al. 2013) was loaded and used to estimate the proportion of immune–stromal component in TME of each sample, expressed in the form of several types of immunological indicators: immune score, StromalScore, ESTIMATEScore and TumorPurity, which positively correlated with the proportion of immune, stromal, the sum of both and tumor, respectively.

#### others
1. We used ELMER32 for understanding which transcription factors are regulated upon perturbations from regulatory regions. Briefly, this method is based on initially identifying differentially methylated distal probes and predicting enriched motifs across them. 

### General Workflow





