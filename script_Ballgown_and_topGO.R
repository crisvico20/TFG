## 001-  Load the modules
library(ballgown) 
library(RColorBrewer) 
library(genefilter) 
library(dplyr) 
library(devtools) 
library(knitr)
library(topGO)

## 002- Load the Phenotypic data, expression data, annotations and sequences
pheno_data = read.csv("PhenoData.csv")
bg = ballgown(dataDir = "Ballgown", samplePattern = "PeameRNASeq_", pData = pheno_data)

## 003- Filter the expression for rows were the FPKM is 1 or more
bg_filt = subset(bg, "rowVars(texpr(bg)) > 1", genomesubset=TRUE)

## 004- Get dataframes with the expression for genes and transcripts
gene_expr = as.data.frame(gexpr(bg_filt))
gene_expr$id = row.names(gene_expr)
transcript_expr = as.data.frame(texpr(bg_filt))
transcript_expr$id = row.names(transcript_expr)

## 005- Get the equivalences for gene_id = gene_name and t_id = t_name
t_id2t_names = bg_filt@expr$trans[,c("t_id", "t_name")]
colnames(t_id2t_names) = c("id", "transcript_name")
t_id2gene_id = bg_filt@expr$trans[,c("t_id", "gene_id")]
gene_id2gene_name = bg_filt@expr$trans[,c("gene_id", "gene_name")]

###########################
### EXPRESSION ANALYSIS ###
###########################

## The expression analysis will be performed in 6 comparisons

## 01- DE for genes and transcripts by treatment (CTR vs INF)
results_genes_by_treatment = stattest(bg_filt, feature="gene", covariate="treatment", adjustvars = c("accession"), getFC=TRUE, meas="FPKM")
results_transcripts_by_treatment = stattest(bg_filt, feature="transcript", covariate="treatment", adjustvars = c("accession"), getFC=TRUE, meas="FPKM")
de_ge_by_treatment = subset(results_genes_by_treatment, results_genes_by_treatment$pval < 0.01)
length(de_ge_by_treatment$id)
## 1202 genes were expressed differentially between treatments
de_te_by_treatment = subset(results_transcripts_by_treatment, results_transcripts_by_treatment$pval < 0.01)
length(de_te_by_treatment$id)
## 1007 transcripts were expressed differentially between treatments
new_results_genes_by_treatment = merge(gene_id2gene_name, results_genes_by_treatment, by.x="gene_id", by.y="id")

## 02- DE for genes and transcripts by accession (BG vs INF)
results_genes_by_accession = stattest(bg_filt, feature="gene", covariate="accession", adjustvars = c("treatment"), getFC=TRUE, meas="FPKM")
results_transcripts_by_accession = stattest(bg_filt, feature="transcript", covariate="accession", adjustvars = c("treatment"), getFC=TRUE, meas="FPKM")
de_ge_by_accession = subset(results_genes_by_accession, results_genes_by_accession$pval < 0.01)
length(de_ge_by_accession$id)
## 3798 genes were expressed differentially between treatments
de_te_by_accession = subset(results_transcripts_by_accession, results_transcripts_by_accession$pval < 0.01)
length(de_te_by_accession$id)
## 3598 transcripts were expressed differentially between treatments
new_results_genes_by_accession = merge(gene_id2gene_name, results_genes_by_accession, by.x="gene_id", by.y="id")

## 03- Now it will get a subset of the samples BG and it will calculate the same (BG: CTR vs INF)
bg_subset_BG = subset(bg_filt, "accession == 'BG'", genomesubset=FALSE)
resultsBG_genes_by_treatment = stattest(bg_subset_BG, feature="gene", covariate="treatment", getFC=TRUE, meas="FPKM")
resultsBG_transcripts_by_treatment = stattest(bg_subset_BG, feature="transcript", covariate="treatment", getFC=TRUE, meas="FPKM")
deBG_ge_by_treatment = subset(resultsBG_genes_by_treatment, resultsBG_genes_by_treatment$pval < 0.01)
length(deBG_ge_by_treatment$id)
## 404 genes were expressed differentially between treatments in BG
deBG_te_by_treatment = subset(resultsBG_transcripts_by_treatment, resultsBG_transcripts_by_treatment$pval < 0.01)
length(deBG_te_by_treatment$id)
## 505 transcripts were expressed differentially between treatments
new_resultsBG_genes_by_treatment = merge(gene_id2gene_name, resultsBG_genes_by_treatment, by.x="gene_id", by.y="id")

## 04- Same thing with DUSA (Dusa CTR vs INF)
bg_subset_DUSA = subset(bg_filt, "accession == 'DUSA'", genomesubset=FALSE)
resultsDUSA_genes_by_treatment = stattest(bg_subset_DUSA, feature="gene", covariate="treatment", getFC=TRUE, meas="FPKM")
resultsDUSA_transcripts_by_treatment = stattest(bg_subset_DUSA, feature="transcript", covariate="treatment", getFC=TRUE, meas="FPKM")
deDUSA_ge_by_treatment = subset(resultsDUSA_genes_by_treatment, resultsDUSA_genes_by_treatment$pval < 0.01)
length(deDUSA_ge_by_treatment$id)
## 800 genes
deDUSA_te_by_treatment = subset(resultsDUSA_transcripts_by_treatment, resultsDUSA_transcripts_by_treatment$pval < 0.01)
length(deDUSA_te_by_treatment$id)
## 877 transcripts
new_resultsDUSA_genes_by_treatment = merge(gene_id2gene_name, resultsDUSA_genes_by_treatment, by.x="gene_id", by.y="id")

## Subset CTR
## 05- CTR (BG vs DUSA)
bg_subset_CTR = subset(bg_filt, "treatment == 'control'", genomesubset=FALSE)
resultsCTR_genes_by_accession = stattest(bg_subset_CTR, feature="gene", covariate="accession", getFC=TRUE, meas="FPKM")
resultsCTR_transcripts_by_accession = stattest(bg_subset_CTR, feature="transcript", covariate="accession", getFC=TRUE, meas="FPKM")
deCTR_ge_by_accession = subset(resultsCTR_genes_by_accession, resultsCTR_genes_by_accession$pval < 0.01)
length(deCTR_ge_by_accession$id)
## 2490 genes 
deCTR_te_by_accession = subset(resultsCTR_transcripts_by_accession, resultsCTR_transcripts_by_accession$pval < 0.01)
length(deCTR_te_by_accession$id)
## 3021 transcripts
new_resultsCTR_genes_by_accession = merge(gene_id2gene_name, resultsCTR_genes_by_accession, by.x="gene_id", by.y="id")

## Subset CTR
## 06- INF (BG vs DUSA)
bg_subset_INF = subset(bg_filt, "treatment == 'infected'", genomesubset=FALSE)
resultsINF_genes_by_accession = stattest(bg_subset_INF, feature="gene", covariate="accession", getFC=TRUE, meas="FPKM")
resultsINF_transcripts_by_accession = stattest(bg_subset_INF, feature="transcript", covariate="accession", getFC=TRUE, meas="FPKM")
deINF_ge_by_accession = subset(resultsINF_genes_by_accession, resultsINF_genes_by_accession$pval < 0.01)
length(deINF_ge_by_accession$id)
## 1396 genes 
deINF_te_by_accession = subset(resultsINF_transcripts_by_accession, resultsINF_transcripts_by_accession$pval < 0.01)
length(deINF_te_by_accession$id)
## 1804 transcripts
new_resultsINF_genes_by_accession = merge(gene_id2gene_name, resultsINF_genes_by_accession, by.x="gene_id", by.y="id")

  
table.genes.de = c(length(de_ge_by_treatment$id),length(deBG_ge_by_treatment$id),length(deDUSA_ge_by_treatment$id),length(deCTR_ge_by_accession$id),length(deINF_ge_by_accession$id),length(de_ge_by_accession$id))
table.transcripts.de = c(length(de_te_by_treatment$id),length(deBG_te_by_treatment$id),length(deDUSA_te_by_treatment$id),length(deCTR_te_by_accession$id),length(deINF_te_by_accession$id),length(de_te_by_accession$id))
table.de = (cbind(table.genes.de,table.transcripts.de))
colnames(table.de) = c("Genes_DE","Transcripts_DE")
rownames(table.de) = c("CTR/INF","BG_CTR/INF","DUSA_CTR/INF","CTR_BG/DUSA","INF_BG/DUSA","BG/DUSA")


results_transcripts_by_treatment = data.frame(geneNames=ballgown::geneNames(bg_filt), geneIDs=ballgown::geneIDs(bg_filt), results_transcripts_by_treatment)
results_transcripts_by_treatment = arrange(results_transcripts_by_treatment,pval)
results_genes_by_treatment = arrange(results_genes_by_treatment,pval)
results_genes_by_accession = arrange(results_genes_by_accession,pval)


## Get the list of intersect IDs and the unique for each category
intersect_ge_ids_treatment = intersect(deBG_ge_by_treatment$id, deDUSA_ge_by_treatment$id);
length(intersect_ge_ids_treatment)
## 25 genes DE in both BD and DUSA
intersect_te_ids_treatment = intersect(deBG_te_by_treatment$id, deDUSA_te_by_treatment$id);
length(intersect_te_ids_treatment)
## 19 transcripts DE in both BD and DUSA
unique_de_ge_BG_ids = intersect(deBG_ge_by_treatment$id, setdiff(deBG_ge_by_treatment$id, deDUSA_ge_by_treatment$id));
length(unique_de_ge_BG_ids)
## 379 genes DE only in BG
unique_de_te_BG_ids = intersect(deBG_te_by_treatment$id, setdiff(deBG_te_by_treatment$id, deDUSA_te_by_treatment$id));
length(unique_de_te_BG_ids)
## 486 transcripts DE only in BG
unique_de_ge_DUSA_ids = intersect(deDUSA_ge_by_treatment$id, setdiff(deDUSA_ge_by_treatment$id, deBG_ge_by_treatment$id));
length(unique_de_ge_DUSA_ids)
## 775 genes DE only in DUSA
unique_de_te_DUSA_ids = intersect(deDUSA_te_by_treatment$id, setdiff(deDUSA_te_by_treatment$id, deBG_te_by_treatment$id));
length(unique_de_te_DUSA_ids)
## 858 transcripts DE only in DUSA

## Get the list of unique genes and transcripts
unique_dege_BG = subset(deBG_ge_by_treatment, deBG_ge_by_treatment$id %in% unique_de_ge_BG_ids)
unique_dete_BG = subset(deBG_te_by_treatment, deBG_te_by_treatment$id %in% unique_de_te_BG_ids)
unique_dege_DUSA = subset(deDUSA_ge_by_treatment, deDUSA_ge_by_treatment$id %in% unique_de_ge_DUSA_ids)
unique_dete_DUSA = subset(deDUSA_te_by_treatment, deDUSA_te_by_treatment$id %in% unique_de_te_DUSA_ids)

## Join the information about expression, annotation and mRNA-seq
unique_dete_BG_exp = merge(merge(unique_dete_BG, transcript_expr, by.y="id"), t_id2t_names, by.y="id")
unique_dete_DUSA_exp = merge(merge(unique_dete_DUSA, transcript_expr, by.y="id"), t_id2t_names, by.y="id")

fpkm = texpr(bg,meas="FPKM")
fpkm = log2(fpkm+1)

par(mar=c(9,9,4,2)) 
boxplot(fpkm, col=c("yellow", "yellow", "yellow", "sienna 1", "sienna 1", "sienna 1", "violet red 3", "violet red 3", "violet red 3", "dark orchid 4", "dark orchid 4", "dark orchid 4"), las=2, ylab='log2(FPKM+1)')

## REPORT ##
## It will create a "supertable" with the expression of all the genes 
sgene_expr_table = gene_expr

sgene_expr_table$MEAN.BGCTR = apply(sgene_expr_table[,c(1:3)], 1, mean)
sgene_expr_table$SD.BGCTR = apply(sgene_expr_table[,c(1:3)], 1, sd)
sgene_expr_table$MEAN.BGINF = apply(sgene_expr_table[,c(4:6)], 1, mean)
sgene_expr_table$SD.BGINF = apply(sgene_expr_table[,c(4:6)], 1, sd)
sgene_expr_table$BG.FC_INFvsCTR = merge(sgene_expr_table, resultsBG_genes_by_treatment, by.y="id")$fc 
sgene_expr_table$BG.pval_INFvsCTR = merge(sgene_expr_table, resultsBG_genes_by_treatment, by.y="id")$pval
sgene_expr_table$BG.qval_INFvsCTR = merge(sgene_expr_table, resultsBG_genes_by_treatment, by.y="id")$qval
sgene_expr_table$BG.sigDE_INFvsCTR = sgene_expr_table$BG.pval_INFvsCTR <= 0.01

sgene_expr_table$MEAN.DUCTR = apply(sgene_expr_table[,c(7:9)], 1, mean)
sgene_expr_table$SD.DUCTR = apply(sgene_expr_table[,c(7:9)], 1, sd)
sgene_expr_table$MEAN.DUINF = apply(sgene_expr_table[,c(10:12)], 1, mean)
sgene_expr_table$SD.DUINF = apply(sgene_expr_table[,c(10:12)], 1, sd)
sgene_expr_table$DU.FC_INFvsCTR = merge(sgene_expr_table, resultsDUSA_genes_by_treatment, by.y="id")$fc 
sgene_expr_table$DU.pval_INFvsCTR = merge(sgene_expr_table, resultsDUSA_genes_by_treatment, by.y="id")$pval
sgene_expr_table$DU.qval_INFvsCTR = merge(sgene_expr_table, resultsDUSA_genes_by_treatment, by.y="id")$qval
sgene_expr_table$DU.sigDE_INFvsCTR = sgene_expr_table$DU.pval_INFvsCTR <= 0.01

sgene_expr_table$CTR.FC_BGvsDU = merge(sgene_expr_table, resultsCTR_genes_by_accession, by.y="id")$fc 
sgene_expr_table$CTR.pval_BGvsDU = merge(sgene_expr_table, resultsCTR_genes_by_accession, by.y="id")$pval
sgene_expr_table$CTR.qval_BGvsDU = merge(sgene_expr_table, resultsCTR_genes_by_accession, by.y="id")$qval
sgene_expr_table$CTR.sigDE_BGvsDU = sgene_expr_table$CTR.pval_BGvsDU <= 0.01

sgene_expr_table$INF.FC_BGvsDU = merge(sgene_expr_table, resultsINF_genes_by_accession, by.y="id")$fc 
sgene_expr_table$INF.pval_BGvsDU = merge(sgene_expr_table, resultsINF_genes_by_accession, by.y="id")$pval
sgene_expr_table$INF.qval_BGvsDU = merge(sgene_expr_table, resultsINF_genes_by_accession, by.y="id")$qval
sgene_expr_table$INF.sigDE_BGvsDU = sgene_expr_table$INF.pval_BGvsDU <= 0.01

sgene_expr_table$Log2FoldChange_BG_INFvsCTR = log2(sgene_expr_table$MEAN.BGINF/sgene_expr_table$MEAN.BGCTR)
sgene_expr_table$Log2FoldChange_DU_INFvsCTR = log2(sgene_expr_table$MEAN.DUINF/sgene_expr_table$MEAN.DUCTR)
sgene_expr_table$Log2FoldChange_CTR_DUvsBG = log2(sgene_expr_table$MEAN.DUCTR/sgene_expr_table$MEAN.BGCTR)
sgene_expr_table$Log2FoldChange_INF_DUvsBG = log2(sgene_expr_table$MEAN.DUINF/sgene_expr_table$MEAN.BGINF)

sgene_expr_table$TREATMENT.FC_CTRvsINF = merge(sgene_expr_table, results_genes_by_treatment, by.y="id")$fc 
sgene_expr_table$TREATMENT.pval_CTRvsINF = merge(sgene_expr_table, results_genes_by_treatment, by.y="id")$pval
sgene_expr_table$TREATMENT.qval_CTRvsINF = merge(sgene_expr_table, results_genes_by_treatment, by.y="id")$qval
sgene_expr_table$TREATMENT.sigDE_CTRvsINF = sgene_expr_table$TREATMENT.pval_CTRvsINF <= 0.01

sgene_expr_table$ACCESSION.FC_BGvsDU = merge(sgene_expr_table, results_genes_by_accession, by.y="id")$fc 
sgene_expr_table$ACCESSION.pval_BGvsDU = merge(sgene_expr_table, results_genes_by_accession, by.y="id")$pval
sgene_expr_table$ACCESSION.qval_BGvsDU = merge(sgene_expr_table, results_genes_by_accession, by.y="id")$qval
sgene_expr_table$ACCESSION.sigDE_BGvsDU = sgene_expr_table$ACCESSION.pval_BGvsDU <= 0.01



## Figures
## 1- Venn Diagram for all the pair comparisons
library("VennDiagram")

n1 = length(sgene_expr_table[sgene_expr_table$BG.sigDE_INFvsCTR == TRUE, 13])
n2 = length(sgene_expr_table[sgene_expr_table$DU.sigDE_INFvsCTR == TRUE, 13])
n3 = length(sgene_expr_table[sgene_expr_table$CTR.sigDE_BGvsDU == TRUE, 13])
n4 = length(sgene_expr_table[sgene_expr_table$INF.sigDE_BGvsDU == TRUE, 13])

n12 = length(sgene_expr_table[sgene_expr_table$BG.sigDE_INFvsCTR == TRUE & sgene_expr_table$DU.sigDE_INFvsCTR == TRUE,13])
n13 = length(sgene_expr_table[sgene_expr_table$BG.sigDE_INFvsCTR == TRUE & sgene_expr_table$CTR.sigDE_BGvsDU == TRUE,13])
n14 = length(sgene_expr_table[sgene_expr_table$BG.sigDE_INFvsCTR == TRUE & sgene_expr_table$INF.sigDE_BGvsDU == TRUE,13])
n23 = length(sgene_expr_table[sgene_expr_table$DU.sigDE_INFvsCTR == TRUE & sgene_expr_table$CTR.sigDE_BGvsDU == TRUE,13])
n24 = length(sgene_expr_table[sgene_expr_table$DU.sigDE_INFvsCTR == TRUE & sgene_expr_table$INF.sigDE_BGvsDU == TRUE,13])
n34 = length(sgene_expr_table[sgene_expr_table$CTR.sigDE_BGvsDU == TRUE & sgene_expr_table$INF.sigDE_BGvsDU == TRUE,13])

n123 = length(sgene_expr_table[sgene_expr_table$BG.sigDE_INFvsCTR == TRUE & sgene_expr_table$DU.sigDE_INFvsCTR == TRUE & sgene_expr_table$CTR.sigDE_BGvsDU == TRUE,13])
n124 = length(sgene_expr_table[sgene_expr_table$BG.sigDE_INFvsCTR == TRUE & sgene_expr_table$DU.sigDE_INFvsCTR == TRUE & sgene_expr_table$INF.sigDE_BGvsDU == TRUE,13])
n134 = length(sgene_expr_table[sgene_expr_table$BG.sigDE_INFvsCTR == TRUE & sgene_expr_table$CTR.sigDE_BGvsDU == TRUE & sgene_expr_table$INF.sigDE_BGvsDU == TRUE,13])
n234 = length(sgene_expr_table[sgene_expr_table$DU.sigDE_INFvsCTR == TRUE & sgene_expr_table$CTR.sigDE_BGvsDU == TRUE & sgene_expr_table$INF.sigDE_BGvsDU == TRUE,13])

n1234 = length(sgene_expr_table[sgene_expr_table$BG.sigDE_INFvsCTR == TRUE & sgene_expr_table$DU.sigDE_INFvsCTR == TRUE & sgene_expr_table$CTR.sigDE_BGvsDU == TRUE & sgene_expr_table$INF.sigDE_BGvsDU == TRUE,13])

draw.quad.venn(n1, n2, n3, n4, n12, n13, n14, n23, n24, n34, n123, n124, n134, n234, n1234, category = c("BG181\nInf. vs Ctr.", "DUSA\nInf. vs Ctr.", "Control\nDUSA vs BG181", "Infectada\nDUSA vs BG181"), fill = c("yellow", "dark orchid 4", "orange red", "deep pink 3"))


## 2- Venn diagram for the intersection between treaments and accessions

s1 = length(sgene_expr_table[sgene_expr_table$TREATMENT.sigDE_CTRvsINF == TRUE, "id"])
s2 = length(sgene_expr_table[sgene_expr_table$ACCESSION.sigDE_BGvsDU == TRUE, "id"])
s12 = length(sgene_expr_table[sgene_expr_table$TREATMENT.sigDE_CTRvsINF == TRUE & sgene_expr_table$ACCESSION.sigDE_BGvsDU == TRUE, "id"])

draw.pairwise.venn(s1, s2, s12, category = c("Tratamiento\n(Ctr. vs Inf.)", "Variedad\n(BG181 vs DUSA)"), fill = c("deep pink 3", "yellow"), cat.pos=1, cat.dist = rep(0.04, 2))


## 3- Heatmaps

DEG_Infected_BGvsDU = sgene_expr_table[sgene_expr_table$INF.sigDE_BGvsDU == TRUE & sgene_expr_table$CTR.sigDE_BGvsDU == FALSE & sgene_expr_table$BG.sigDE_INFvsCTR == TRUE & sgene_expr_table$DU.sigDE_INFvsCTR == FALSE, c(14, 16, 22, 24)]
DEG_Infected_BGvsDU2 = sgene_expr_table[sgene_expr_table$TREATMENT.sigDE_CTRvsINF == TRUE & log2(sgene_expr_table$DU.FC_INFvsCTR) <= 1 & log2(sgene_expr_table$DU.FC_INFvsCTR) >= -1, c(14, 16, 22, 24)]
DEG_Infected_BGvsDU3 = sgene_expr_table[sgene_expr_table$ACCESSION.sigDE_BGvsDU == TRUE & sgene_expr_table$TREATMENT.sigDE_CTRvsINF == TRUE & sgene_expr_table$DU.FC_INFvsCTR < 2, c(1:45)]

DEG_Infected_BGvsDU = sgene_expr_table[sgene_expr_table$ACCESSION.pval_BGvsDU <= 0.01 & sgene_expr_table$TREATMENT.pval_CTRvsINF <= 0.01 & sgene_expr_table$Log2FoldChange_BG_INFvsCTR < 0, c(14, 16, 22, 24)]
DEG_Infected_BGvsDU_ids = sgene_expr_table[sgene_expr_table$ACCESSION.pval_BGvsDU <= 0.01 & sgene_expr_table$TREATMENT.pval_CTRvsINF <= 0.01 & sgene_expr_table$Log2FoldChange_BG_INFvsCTR < 0, 13]
DEG_Infected_BGvsDU_T = log10(DEG_Infected_BGvsDU)
DEG_Infected_BGvsDU_T$id = DEG_Infected_BGvsDU_ids

## It will plot with pheatmap

library(pheatmap)

## 3.1- Select the data to be plotted.
##      It will select the intersection between Accession and Treatment.
##      These are 227 genes, so it will divide them into two sets
##      A- Where log2 DUSA_INF/BG_INF > 1
##      B- Where log2 DUSA_INF/BG_INF < -1

## Scale the row data
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

DEGs_ISect_AccessionTreatment_matrix_Up = as.matrix(sgene_expr_table[sgene_expr_table$TREATMENT.sigDE_CTRvsINF == TRUE & sgene_expr_table$ACCESSION.sigDE_BGvsDU == TRUE & sgene_expr_table$Log2FoldChange_INF_DUvsBG > 0, c("MEAN.BGCTR", "MEAN.BGINF", "MEAN.DUCTR", "MEAN.DUINF")])
rownames(DEGs_ISect_AccessionTreatment_matrix_Up) = sgene_expr_table[sgene_expr_table$TREATMENT.sigDE_CTRvsINF == TRUE & sgene_expr_table$ACCESSION.sigDE_BGvsDU == TRUE & sgene_expr_table$Log2FoldChange_INF_DUvsBG > 0, c("id")]
colnames(DEGs_ISect_AccessionTreatment_matrix_Up) = c("BG181_Ctr", "BG181_Inf", "DUSA_Ctr", "DUSA_Inf")

DEGs_ISect_AccessionTreatment_norm_up <- t(apply(DEGs_ISect_AccessionTreatment_matrix_Up, 1, cal_z_score))
pheatmap(DEGs_ISect_AccessionTreatment_norm_up, cluster_cols = FALSE)


DEGs_ISect_AccessionTreatment_matrix_Down = as.matrix(sgene_expr_table[sgene_expr_table$TREATMENT.sigDE_CTRvsINF == TRUE & sgene_expr_table$ACCESSION.sigDE_BGvsDU == TRUE & sgene_expr_table$Log2FoldChange_INF_DUvsBG < 0, c("MEAN.BGCTR", "MEAN.BGINF", "MEAN.DUCTR", "MEAN.DUINF")])
rownames(DEGs_ISect_AccessionTreatment_matrix_Down) = sgene_expr_table[sgene_expr_table$TREATMENT.sigDE_CTRvsINF == TRUE & sgene_expr_table$ACCESSION.sigDE_BGvsDU == TRUE & sgene_expr_table$Log2FoldChange_INF_DUvsBG < 0, c("id")]
colnames(DEGs_ISect_AccessionTreatment_matrix_Down) = c("BG181_Ctr", "BG181_Inf", "DUSA_Ctr", "DUSA_Inf")

DEGs_ISect_AccessionTreatment_norm_down <- t(apply(DEGs_ISect_AccessionTreatment_matrix_Down, 1, cal_z_score))
pheatmap(DEGs_ISect_AccessionTreatment_norm_down, cluster_cols = FALSE)



#########################
### GO Terms Analysis ###
#########################

geneID2GO <- readMappings(file = "Universe/Persea_americana_annos1_hass.proteins.go_terms.txt")
GO2geneID <- inverseList(geneID2GO)
geneNames = names(geneID2GO)
clean_geneNames <- gsub('.{0,6}$','',geneNames)

DEAccession_GeneIDList = new_results_genes_by_accession[new_results_genes_by_accession$pval < 0.01,2]
DETreatment_GeneIDList = new_results_genes_by_treatment[new_results_genes_by_treatment$pval < 0.01,2]
DEIntersectAccTre_GeneIDList = intersect(DEAccession_GeneIDList, DETreatment_GeneIDList)

geneList4DEAccession = factor(as.integer(clean_geneNames %in% DEAccession_GeneIDList))
names(geneList4DEAccession) = geneNames
GOData4Accession_BP = new("topGOdata", ontology = "BP", allGenes = geneList4DEAccession, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOData4Accession_CC = new("topGOdata", ontology = "CC", allGenes = geneList4DEAccession, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOData4Accession_MF = new("topGOdata", ontology = "MF", allGenes = geneList4DEAccession, annot = annFUN.gene2GO, gene2GO = geneID2GO)

geneList4DETreatment = factor(as.integer(clean_geneNames %in% DETreatment_GeneIDList))
names(geneList4DETreatment) = geneNames
GOData4Treatment_BP = new("topGOdata", ontology = "BP", allGenes = geneList4DETreatment, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOData4Treatment_CC = new("topGOdata", ontology = "CC", allGenes = geneList4DETreatment, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOData4Treatment_MF = new("topGOdata", ontology = "MF", allGenes = geneList4DETreatment, annot = annFUN.gene2GO, gene2GO = geneID2GO)

geneList4DEIntersectAccTre = factor(as.integer(clean_geneNames %in% DEIntersectAccTre_GeneIDList))
names(geneList4DEIntersectAccTre) = geneNames
GOData4IntersectAccTre_BP = new("topGOdata", ontology = "BP", allGenes = geneList4DEIntersectAccTre, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOData4IntersectAccTre_CC = new("topGOdata", ontology = "CC", allGenes = geneList4DEIntersectAccTre, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOData4IntersectAccTre_MF = new("topGOdata", ontology = "MF", allGenes = geneList4DEIntersectAccTre, annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultFis_IAC_BP <- runTest(GOData4IntersectAccTre_BP, algorithm = "classic", statistic = "fisher")
resultWeight_IAC_BP <- runTest(GOData4IntersectAccTre_BP, algorithm = "weight01", statistic = "fisher")
resultKS_IAC_BP <- runTest(GOData4IntersectAccTre_BP, algorithm = "elim", statistic = "ks")

allRes_IAC_BP <- GenTable(GOData4IntersectAccTre_BP, classic = resultFis_IAC_BP, KS = resultKS_IAC_BP, weight = resultWeight_IAC_BP, orderBy = "weight", ranksOf = "classic", topNodes = 30)
write.csv(allRes_IAC_BP, "GenTable_Intersect_Acc_Tre_BP.csv")

## Lets perform the GO terms enrichment analysis with each of the categories

library(ggplot2)

gene_ids_DE_BG_INFvsCTR = new_resultsBG_genes_by_treatment[new_resultsBG_genes_by_treatment$pval < 0.01,2]
gene_ids_DE_DU_INFvsCTR = new_resultsDUSA_genes_by_treatment[new_resultsDUSA_genes_by_treatment$pval < 0.01,2]
gene_ids_DE_CTR_BGvsDU = new_resultsCTR_genes_by_accession[new_resultsCTR_genes_by_accession$pval < 0.01,2]
gene_ids_DE_INF_BGvsDU = new_resultsINF_genes_by_accession[new_resultsINF_genes_by_accession$pval < 0.01,2]


## 1- BG181 Infected vs Control

gene_list4DE_BG_INFvsCTR = factor(as.integer(clean_geneNames %in% gene_ids_DE_BG_INFvsCTR))
names(gene_list4DE_BG_INFvsCTR) = geneNames
GOData4DE_BG_INFvsCTR_BP = new("topGOdata", ontology = "BP", allGenes = gene_list4DE_BG_INFvsCTR, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOData4DE_BG_INFvsCTR_CC = new("topGOdata", ontology = "CC", allGenes = gene_list4DE_BG_INFvsCTR, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOData4DE_BG_INFvsCTR_MF = new("topGOdata", ontology = "MF", allGenes = gene_list4DE_BG_INFvsCTR, annot = annFUN.gene2GO, gene2GO = geneID2GO)

GOData4DE_BG_INFvsCTR_BP_resultCLFI = runTest(GOData4DE_BG_INFvsCTR_BP, algorithm = "classic", statistic = "fisher")
GOData4DE_BG_INFvsCTR_BP_resultCLWE = runTest(GOData4DE_BG_INFvsCTR_BP, algorithm = "weight01", statistic = "fisher")
GOData4DE_BG_INFvsCTR_BP_resultELKS = runTest(GOData4DE_BG_INFvsCTR_BP, algorithm = "elim", statistic = "ks")

pVals_GSEA_BG_INFvsCTR_BP_resultCLFI = score(GOData4DE_BG_INFvsCTR_BP_resultCLFI)[score(GOData4DE_BG_INFvsCTR_BP_resultCLFI) <= 0.05]
pVals_GSEA_BG_INFvsCTR_BP_resultCLWE = score(GOData4DE_BG_INFvsCTR_BP_resultCLWE)[score(GOData4DE_BG_INFvsCTR_BP_resultCLWE) <= 0.05]
pVals_GSEA_BG_INFvsCTR_BP_resultELKS = score(GOData4DE_BG_INFvsCTR_BP_resultELKS)[score(GOData4DE_BG_INFvsCTR_BP_resultELKS) <= 0.05]

GOData4DE_BG_INFvsCTR_BP_StatsCLWE = termStat(object = GOData4DE_BG_INFvsCTR_BP, whichGO = names(pVals_GSEA_BG_INFvsCTR_BP_resultCLWE))
GOData4DE_BG_INFvsCTR_BP_StatsCLWE$DEG = GOData4DE_BG_INFvsCTR_BP_StatsCLWE$Significant
GOData4DE_BG_INFvsCTR_BP_StatsCLWE$pValue = pVals_GSEA_BG_INFvsCTR_BP_resultCLWE
GOData4DE_BG_INFvsCTR_BP_StatsCLWE$Term = Term(rownames(GOData4DE_BG_INFvsCTR_BP_StatsCLWE))
ggplot(GOData4DE_BG_INFvsCTR_BP_StatsCLWE, aes(x = DEG/Expected, y = Term)) + geom_point(aes(color=pValue, size=DEG)) + scale_color_gradient(low="blue", high="yellow")

### 2- DUSA Infected vs Control

gene_list4DE_DU_INFvsCTR = factor(as.integer(clean_geneNames %in% gene_ids_DE_DU_INFvsCTR))
names(gene_list4DE_DU_INFvsCTR) = geneNames
GOData4DE_DU_INFvsCTR_BP = new("topGOdata", ontology = "BP", allGenes = gene_list4DE_DU_INFvsCTR, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOData4DE_DU_INFvsCTR_CC = new("topGOdata", ontology = "CC", allGenes = gene_list4DE_DU_INFvsCTR, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOData4DE_DU_INFvsCTR_MF = new("topGOdata", ontology = "MF", allGenes = gene_list4DE_DU_INFvsCTR, annot = annFUN.gene2GO, gene2GO = geneID2GO)

GOData4DE_DU_INFvsCTR_BP_resultCLFI = runTest(GOData4DE_DU_INFvsCTR_BP, algorithm = "classic", statistic = "fisher")
GOData4DE_DU_INFvsCTR_BP_resultCLWE = runTest(GOData4DE_DU_INFvsCTR_BP, algorithm = "weight01", statistic = "fisher")
GOData4DE_DU_INFvsCTR_BP_resultELKS = runTest(GOData4DE_DU_INFvsCTR_BP, algorithm = "elim", statistic = "ks")

pVals_GSEA_DU_INFvsCTR_BP_resultCLFI = score(GOData4DE_DU_INFvsCTR_BP_resultCLFI)[score(GOData4DE_DU_INFvsCTR_BP_resultCLFI) <= 0.05]
pVals_GSEA_DU_INFvsCTR_BP_resultCLWE = score(GOData4DE_DU_INFvsCTR_BP_resultCLWE)[score(GOData4DE_DU_INFvsCTR_BP_resultCLWE) <= 0.05]
pVals_GSEA_DU_INFvsCTR_BP_resultELKS = score(GOData4DE_DU_INFvsCTR_BP_resultELKS)[score(GOData4DE_DU_INFvsCTR_BP_resultELKS) <= 0.05]

GOData4DE_DU_INFvsCTR_BP_StatsCLWE = termStat(object = GOData4DE_DU_INFvsCTR_BP, whichGO = names(pVals_GSEA_DU_INFvsCTR_BP_resultCLWE))
GOData4DE_DU_INFvsCTR_BP_StatsCLWE$DEG = GOData4DE_DU_INFvsCTR_BP_StatsCLWE$Significant
GOData4DE_DU_INFvsCTR_BP_StatsCLWE$pValue = pVals_GSEA_DU_INFvsCTR_BP_resultCLWE
GOData4DE_DU_INFvsCTR_BP_StatsCLWE$Term = Term(rownames(GOData4DE_DU_INFvsCTR_BP_StatsCLWE))
ggplot(GOData4DE_DU_INFvsCTR_BP_StatsCLWE, aes(x = DEG/Expected, y = Term)) + geom_point(aes(color=pValue, size=DEG)) + scale_color_gradient(low="blue", high="yellow")

## 3- Control for BG181 vs DUSA

gene_list4DE_CTR_BGvsDU = factor(as.integer(clean_geneNames %in% gene_ids_DE_CTR_BGvsDU))
names(gene_list4DE_CTR_BGvsDU) = geneNames
GOData4DE_CTR_BGvsDU_BP = new("topGOdata", ontology = "BP", allGenes = gene_list4DE_CTR_BGvsDU, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOData4DE_CTR_BGvsDU_CC = new("topGOdata", ontology = "CC", allGenes = gene_list4DE_CTR_BGvsDU, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOData4DE_CTR_BGvsDU_MF = new("topGOdata", ontology = "MF", allGenes = gene_list4DE_CTR_BGvsDU, annot = annFUN.gene2GO, gene2GO = geneID2GO)

GOData4DE_CTR_BGvsDU_BP_resultCLFI = runTest(GOData4DE_CTR_BGvsDU_BP, algorithm = "classic", statistic = "fisher")
GOData4DE_CTR_BGvsDU_BP_resultCLWE = runTest(GOData4DE_CTR_BGvsDU_BP, algorithm = "weight01", statistic = "fisher")
GOData4DE_CTR_BGvsDU_BP_resultELKS = runTest(GOData4DE_CTR_BGvsDU_BP, algorithm = "elim", statistic = "ks")

pVals_GSEA_CTR_BGvsDU_BP_resultCLFI = score(GOData4DE_CTR_BGvsDU_BP_resultCLFI)[score(GOData4DE_CTR_BGvsDU_BP_resultCLFI) <= 0.05]
pVals_GSEA_CTR_BGvsDU_BP_resultCLWE = score(GOData4DE_CTR_BGvsDU_BP_resultCLWE)[score(GOData4DE_CTR_BGvsDU_BP_resultCLWE) <= 0.05]
pVals_GSEA_CTR_BGvsDU_BP_resultELKS = score(GOData4DE_CTR_BGvsDU_BP_resultELKS)[score(GOData4DE_CTR_BGvsDU_BP_resultELKS) <= 0.05]

GOData4DE_CTR_BGvsDU_BP_StatsCLWE = termStat(object = GOData4DE_CTR_BGvsDU_BP, whichGO = names(pVals_GSEA_CTR_BGvsDU_BP_resultCLWE))
GOData4DE_CTR_BGvsDU_BP_StatsCLWE$DEG = GOData4DE_CTR_BGvsDU_BP_StatsCLWE$Significant
GOData4DE_CTR_BGvsDU_BP_StatsCLWE$pValue = pVals_GSEA_CTR_BGvsDU_BP_resultCLWE
GOData4DE_CTR_BGvsDU_BP_StatsCLWE$Term = Term(rownames(GOData4DE_CTR_BGvsDU_BP_StatsCLWE))
ggplot(GOData4DE_CTR_BGvsDU_BP_StatsCLWE, aes(x = DEG/Expected, y = Term)) + geom_point(aes(color=pValue, size=DEG)) + scale_color_gradient(low="blue", high="yellow")

## 4- Infected for BG181 vs DUSA

gene_list4DE_INF_BGvsDU = factor(as.integer(clean_geneNames %in% gene_ids_DE_INF_BGvsDU))
names(gene_list4DE_INF_BGvsDU) = geneNames
GOData4DE_INF_BGvsDU_BP = new("topGOdata", ontology = "BP", allGenes = gene_list4DE_INF_BGvsDU, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOData4DE_INF_BGvsDU_CC = new("topGOdata", ontology = "CC", allGenes = gene_list4DE_INF_BGvsDU, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOData4DE_INF_BGvsDU_MF = new("topGOdata", ontology = "MF", allGenes = gene_list4DE_INF_BGvsDU, annot = annFUN.gene2GO, gene2GO = geneID2GO)

GOData4DE_INF_BGvsDU_BP_resultCLFI = runTest(GOData4DE_INF_BGvsDU_BP, algorithm = "classic", statistic = "fisher")
GOData4DE_INF_BGvsDU_BP_resultCLWE = runTest(GOData4DE_INF_BGvsDU_BP, algorithm = "weight01", statistic = "fisher")
GOData4DE_INF_BGvsDU_BP_resultELKS = runTest(GOData4DE_INF_BGvsDU_BP, algorithm = "elim", statistic = "ks")

pVals_GSEA_INF_BGvsDU_BP_resultCLFI = score(GOData4DE_INF_BGvsDU_BP_resultCLFI)[score(GOData4DE_INF_BGvsDU_BP_resultCLFI) <= 0.05]
pVals_GSEA_INF_BGvsDU_BP_resultCLWE = score(GOData4DE_INF_BGvsDU_BP_resultCLWE)[score(GOData4DE_INF_BGvsDU_BP_resultCLWE) <= 0.05]
pVals_GSEA_INF_BGvsDU_BP_resultELKS = score(GOData4DE_INF_BGvsDU_BP_resultELKS)[score(GOData4DE_INF_BGvsDU_BP_resultELKS) <= 0.05]

GOData4DE_INF_BGvsDU_BP_StatsCLWE = termStat(object = GOData4DE_INF_BGvsDU_BP, whichGO = names(pVals_GSEA_INF_BGvsDU_BP_resultCLWE))
GOData4DE_INF_BGvsDU_BP_StatsCLWE$DEG = GOData4DE_INF_BGvsDU_BP_StatsCLWE$Significant
GOData4DE_INF_BGvsDU_BP_StatsCLWE$pValue = pVals_GSEA_INF_BGvsDU_BP_resultCLWE
GOData4DE_INF_BGvsDU_BP_StatsCLWE$Term = Term(rownames(GOData4DE_INF_BGvsDU_BP_StatsCLWE))
ggplot(GOData4DE_INF_BGvsDU_BP_StatsCLWE, aes(x = DEG/Expected, y = Term)) + geom_point(aes(color=pValue, size=DEG)) + scale_color_gradient(low="blue", high="yellow")

