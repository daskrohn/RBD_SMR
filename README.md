# RBD_SMR
Using SMR for eQTL mendelian randomization.  
From Sara Bandres-Ciga:  
*The SMR software tool implements the SMR & HEIDI methods to test for pleiotropic association 
between the expression level of a gene and a complex trait of interest using summary-level data 
from GWAS and expression quantitative trait loci (eQTL) studies (Zhu et al. 2016 Nature Genetics). 
The methodology can be interpreted as an analysis to test if the effect size of a SNP on the 
phenotype is mediated by gene expression. This tool can therefore be used to prioritize genes 
underlying GWAS hits for follow-up functional studies.* 

## Set up workspace
SMR 
I used eQTLgen, GTEx v7 substantia nigra, and Brain eMeta (expression) and mMeta (methylation). SMR files were downloaded from FILL IN. 

## Preparation of files
First, IF you are working with a pathway or specific set of genes, create a gene list for use with the QTL summary stats.  
* MR_list.txt is a single column list of genes you are interested in based on GWAS summary stats.  
* genIdsFromHugo.txt is a reference file included in this repository. 

*You do not need this step if you are running with full GWAS summary statistics.*
```R
require(data.table)
geneList <- fread("MR_list.txt", header = F)
names(geneList) <- "GENE"
ensg <- fread("genIdsFromHugo.txt", header = T, sep = "\t")
data <- merge(geneList, ensg, by.x = "GENE", by.y = "Approved Symbol")
names(data)[12] <- "ensgOut"
write.table(data$ensgOut, "ensembl_MR_list.txt", quote = F, sep = "\t", row.names = F, col.names = F)
```

Next, create working data plink files, and convert plink .bim from chr:pos format to rsID format (to match the QTL summary stats):
```R
# R

require(data.table)
bim <- fread("~/runs/krohn/krohn/Main_0_Data/NeuroX/Imputed/imptd_RBD_Nx.bim")
ref <- fread("/home/daskrohn/runs/krohn/krohn/Main_0_Data/HRC.r1.1.GRCh37.for-anno.tab")

colnames(ref) = c("MarkerID", "rsID")
colnames(bim) = c("CHR", "MarkerID", "dis", "BP", "A2", "A1")

merged <- merge(bim, ref, by="MarkerID", all.x=TRUE)

info = merged[,c("CHR", "rsID", "dis", "BP", "A2", "A1", "MarkerID")]

write.table(info, file="rsID_info.for-bim.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, sep = "\t")

q()
````
Since some positions don't have an rsID, they must be removed from the bfiles. 

````
#!/bin/bash

awk '{if ($2==".") print $7} rsID_info.for-bim.txt > exclude_noRS.txt

plink --bfile ~/runs/krohn/krohn/Main_0_Data/NeuroX/Imputed/imptd_RBD_Nx \
--exclude exclude_noRS.txt \
--make-bed --out workingData

awk '{if ($2!=".") print $7,$2}' rsID_info.for-bim.txt > rename_bim.txt

plink --bfile workingData --update-name rename_bim.txt --make-bed --out workingData 
````
## Run SMR
````
./smr_Linux --bfile workingData \
--smr-multi \
--gwas-summary toSMR_summarystats.txt \
--beqtl-summary eQTLGen/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense \
--out eQTLGen_multi_RBD-meta \
--thread-num 12

./smr_Linux --bfile workingData \
--smr-multi \
--gwas-summary toSMR_summarystats.txt \
--beqtl-summary Brain-eMeta/Brain-eMeta \
--out Brain-eMeta_multi_RBD-meta \
--thread-num 12

./smr_Linux --bfile workingData \
--smr-multi \
--gwas-summary toSMR_summarystats.txt \
--beqtl-summary Brain-mMeta/Brain-mMeta \
--out Brain-mMeta_multi_RBD-meta \
--thread-num 12

./smr_Linux --bfile workingData \
--smr-multi \
--gwas-summary toSMR_summarystats.txt \
--beqtl-summary GTEx_V7_cis_eqtl_summary_lite/Brain_Substantia_nigra_1e-05 \
--out GTEx-Brain_Substantia_nigra_1e-05_RBD-meta \
--thread-num 12

mkdir results
mkdir logs 

mv *.msmr results
mv *.list logs 
````
