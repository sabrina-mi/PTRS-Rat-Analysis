library(tidyverse)
"%&%" = function(a,b) paste(a,b,sep="")


# Path to directory with genotype, gene expression, and phyMap files
data.dir <- "/home/s1mi/Github/PTRS-Rat-Analysis/data/"
# Create new directory for GEMMA inputs and outputs
base.dir <- "/home/s1mi/Github/PTRS-Rat-Analysis/GEMMA/"

tis="Br"


bimfile <- base.dir %&% tis %&% "_bimbam" ###get SNP position information###

pheno.dir <- base.dir %&% "phenotype_files/"
ge.dir <- base.dir %&% "genotype_files/"
grm.dir <- base.dir %&% "output/"

bim <- read.table(bimfile, header=TRUE)
gex <- read.csv(data.dir %&% "gexBr.csv", header=TRUE)
gtf <- readRDS(data.dir %&% "gene_annotation.RDS")

ensidlist <- colnames(gex)[-1]
setwd(base.dir)
fn_run_gemma_grm = function(tis){
  #Read in bimbam file
  bimfile <- base.dir %&% tis %&% "_bimbam" ###get SNP position information###

  pheno.dir <- base.dir %&% "phenotype_files/"
  ge.dir <- base.dir %&% "genotype_files/"
  grm.dir <- base.dir %&% "output/"

  bim <- read.table(bimfile, header=TRUE)
  gex <- read.csv(data.dir %&% "phenotypes/gexBr.csv", header=TRUE)
  gtf <- readRDS(data.dir %&% "gene_annotation.RDS")


  ensidlist <- colnames(gex)[-1]
  setwd(base.dir)
  for(i in 1:length(ensidlist)){
    cat(i,"/",length(ensidlist),"\n")
    gene <- ensidlist[i]
    geneinfo <- gtf[match(gene, rownames(gtf)),]
    c <- geneinfo$chr
    start <- geneinfo$start - 1e6 ### 1Mb lower bound for cis-eQTLS
    end <- geneinfo$end + 1e6 ### 1Mb upper bound for cis-eQTLs
    chrsnps <- subset(bim, bim[,1]==c) ### pull snps on same chr
    lower_bound <- subset(chrsnps, chrsnps[,2]>=start)
    cissnps <- subset(lower_bound, chrsnps[,2]<=end) ### pull cis-SNP info
    if (nrow(cissnps)){
      snplist <- cissnps[,3:ncol(cissnps)]
      write.table(snplist, file= ge.dir %&% "tmp." %&% tis %&% ".geno" %&% gene, quote=F,col.names=F,row.names=F)
      geneexp <- cbind(as.numeric(gex[,i+1]))
      write.table(geneexp, file= pheno.dir %&% "tmp.pheno." %&% gene, col.names=F, row.names = F, quote=F) #output pheno for gemma
    }

  }
}

fn_run_gemma_grm(tis)

