# RBD_SMR
Using SMR for eQTL mendelian randomization. 

## Downloading data. 
I used eQTLgen, GTEx v7 substantia nigra, and Brain eMeta (expression) and mMeta (methylation). SMR files were downloaded from FILL IN. 

## Preparation of files
Creating working data plink files, and converting plink bfiles from chr:pos format to rsID format:
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
--out eQTLGen_multi_RBD_PATHWAYS \
--genes ensembl_MR_list.txt \
--thread-num 12

./smr_Linux --bfile workingData \
--smr-multi \
--gwas-summary toSMR_summarystats.txt \
--beqtl-summary Brain-eMeta/Brain-eMeta \
--out Brain-eMeta_multi_RBD_PATHWAYS \
--genes ensembl_MR_list.txt \
--thread-num 12

./smr_Linux --bfile workingData \
--smr-multi \
--gwas-summary toSMR_summarystats.txt \
--beqtl-summary Brain-mMeta/Brain-mMeta \
--out Brain-mMeta_multi_RBD_PATHWAYS \
--genes ensembl_MR_list.txt \
--thread-num 12

./smr_Linux --bfile workingData \
--smr-multi \
--gwas-summary toSMR_summarystats.txt \
--beqtl-summary GTEx_V7_cis_eqtl_summary_lite/Brain_Substantia_nigra_1e-05 \
--out GTEx-Brain_Substantia_nigra_1e-05_PD_PATHWAYS \
--genes ensembl_MR_list.txt \
--thread-num 12
````
