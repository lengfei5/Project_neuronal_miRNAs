##########################################################################
##########################################################################
## Project: Chiara's neuron calss-specific miRNA expression
## Script purpose: to check the functional enrichment for miRNAs in two different neuron groups
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Wed May 30 09:57:06 2018
##########################################################################
##########################################################################
#RdataDir = paste0("../results/tables_for_decomvolution/Rdata/")
dataDir = "../results/decomvolution_results/functional_enrich"
resDir = "../results/decomvolution_results/functional_enrich"
if(!dir.exists(resDir)) dir.create(resDir)

library(openxlsx)
groups = read.xlsx(paste0(dataDir, "/miRNA_signature_4_targets_conservation.xlsx"), 
                   sheet = 1, colNames = TRUE, skipEmptyRows = TRUE, skipEmptyCols = TRUE)
groups = groups[which(!is.na(groups$conservation)),]
groups$gene = sapply(groups$gene, function(x) gsub('r-', "R-", x), USE.NAMES = FALSE)

target.scan.db = read.delim("../../../../Databases/TargetScan_Worm_6.2/Summary_Counts.txt",
                            header = TRUE, sep = "\t")
target.db = target.scan.db[which(target.scan.db$Total.num.conserved.sites>=2), ]

family.info = read.delim("../../../../Databases/TargetScan_Worm_6.2/miR_Family_Info.txt",
                          header = TRUE, sep = "\t")
gps = unique(groups$neuron.group)


Functional.Enrichment.Analysis = TRUE
if(Functional.Enrichment.Analysis){
  library(clusterProfiler)
  library(org.Ce.eg.db)
  keytypes(org.Ce.eg.db)
  
  fdr.cutoff = 0.05
  pvalue.cutoff = 0.01
  fun.enrich = "enrichKEGG"
  #fun.enrich = "enrichGO"
  
}

for(n in 1:length(gps))
{
  # n = 2
  targets = c()
  jj = which(groups$neuron.group == gps[n])
  #mm = match(target.db$Representative.miRNA, gr)
  #targets = 
  for(j in jj) 
  {
    #kk = grep(groups$gene[j], target.db$Representative.miRNA)
    #nn = grep(groups$gene[j], family.info$MiRBase.ID)
    nn = which(family.info$MiRBase.ID == paste0("cel-", groups$gene[j]))
    if(length(nn) == 1) {
      kk = which(target.db$miRNA.family ==  as.character(family.info$Seed.m8[nn]))
      if(length(kk)>0) {
        cat(groups$gene[j], 'familly -- ',  as.character(family.info$miR.family[nn]), "--", as.character(unique(target.db$Representative.miRNA[kk])),  
            length(kk), "targets",   "\n")
        targets = c(targets, unique(as.character(target.db$Gene.Symbol[kk])))
      }else{
        cat(groups$gene[j], 'familly -- ',  as.character(family.info$miR.family[nn]), "--", as.character(unique(target.db$Representative.miRNA[kk])),  
            length(kk), "targets",   "\n")
        #cat("No targets found \n")
      }
    }else{
      cat(groups$gene[j],  " ERROR : nb of family found", length(nn), "\n")
    }
  }
  
  targets = unique(targets)
  
  #write.table(targets, file = paste0(resDir, "Targets_for_", gps[n], ".txt"), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  if(Functional.Enrichment.Analysis)
  {
    save.file = paste0(resDir, '/FunctionalEnrichment_analysis_fdr_', fdr.cutoff, '_pval_', pvalue.cutoff, "_", gps[n], '.xlsx')
    save.pdf = paste0(resDir, '/FunctionalEnrichment_analysis_fdr_', fdr.cutoff, '_pval_', pvalue.cutoff, "_", gps[n], '.pdf')
    
    ids = bitr(targets, fromType="SYMBOL", toType=c("ENTREZID", "PMID", "WORMBASE",  "UNIPROT"), OrgDb="org.Ce.eg.db", drop = FALSE);
    #ids.kegg = bitr(ids$WORMBASE, fromType = "WORMBASE", toType = "kegg", OrgDb="org.Ce.eg.db", drop = FALSE);
        
    ego <- enrichGO(gene         = ids$ENTREZID,
                     OrgDb         = org.Ce.eg.db,
                     keyType       = 'ENTREZID',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = pvalue.cutoff, 
                    qvalueCutoff  = fdr.cutoff,
                     readable      = TRUE)
    #ego = as.data.frame(ego)
    
    head(ego)
    #search_kegg_organism('cel', by='kegg_code')
    ekegg <- enrichKEGG(gene         = ids$UNIPROT,
                     keyType = "uniprot",
                     organism     = 'cel',
                     pvalueCutoff = pvalue.cutoff,
                     qvalueCutoff  = fdr.cutoff)
    #ekegg = as.data.frame(ekegg)
    head(ekegg)
    
    pdf(save.pdf, width = 10, height = 8)
    head(as.data.frame(ego))
    write.xlsx(as.data.frame(ego), file= save.file,
               col.names = TRUE, row.names = TRUE, sheetName = "GO_term")
    
    barplot(ego, drop=TRUE, showCategory=20)
    #dotplot(, showCategory=max(as.numeric(table(as.data.frame(ck)[,1]))))
    
    if(nrow(ekegg)>0){
      barplot(ekegg, drop=TRUE, showCategory=20)
      write.xlsx(as.data.frame(ekegg), file= save.file,
                 col.names = TRUE, row.names = TRUE, sheetName = "kegg")
      
    }
   
    dev.off()
    
  }
  
}
