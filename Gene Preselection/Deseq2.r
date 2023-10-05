

readexprs <- function(path, lbl){
    #read gene expression matrices
    df = read.csv(path , sep=',')
    rownames(df) = df$X
    df = df[, -1]
    
    names = paste0(lbl, 1:ncol(df))
    colnames(df) = names
    
    return(df)
}

path_r = readlines("Enter file path for resistance patients: ")
path_s = readlines("Enter file path for sensitive patients: ")

R = readexprs(path_r, "R")
S = readexprs(path_s, "S")

prep = function(case, control){
    #Prepare column data for DESeq2 
    library(dplyr)
    
    all = cbind(case, control)
    
    coldata_late = data.frame(Sample=colnames(case), Group=rep(c('C'), each=ncol(case)))
    coldata_early = data.frame(Sample=colnames(control), Group=rep(c('N'), each=ncol(control)))
    
    coldata = rbind(coldata_late, coldata_early)
    rownames(coldata) = coldata$Sample
    
    return(coldata)
}

coldf = prep(R,S)

de = function(sensitive, resistance, coldata){
    #Perform differential gene expression analysis
    library(DESeq2)
    library(dplyr)
    library(devtools)
    library(BiocParallel)
    library(data.table)
    
    register(MulticoreParam())
    
    dataset = cbind(resistance, sensitive)
    
    dds = DESeqDataSetFromMatrix(countData = round(dataset), 
                                 colData = coldata,
                                 design = ~ Group)
    
    dds$Group = relevel(dds$Group, ref = "R")
    
    keep = rowSums(counts(dds)) >= 10
    dds = dds[keep,]
        
    dds_Deseq = DESeq(dds, parallel=TRUE, fitType = "local")    
    res = results(dds_Deseq, contrast=c("Group", "R", "S"), alpha=0.05, parallel=TRUE)
    resOrdered = res[order(res$padj), ]
    
    all_genes = as.character(rownames(resOrdered))
    all_genes = sub("\\.\\d+", "", all_genes)
    rownames(resOrdered) = all_genes
    
    up = resOrdered[which(resOrdered$log2FoldChange > 1), ]
    up = up[which(up$padj < 0.05), ]
    
    down = resOrdered[which(resOrdered$log2FoldChange < -1), ]
    down = down[which(down$padj < 0.05), ]

    up_genes = resistance[rownames(up), ]
    down_genes = resistance[rownames(down), ]
    all = rbind(up_genes, down_genes)
    
    all_exprs = varistran::vst(all)
    
    return(all_exprs)

}

df_exprs = de(S,R,coldf)


TOM_Matrices = function(dataset){

  library(WGCNA)

  enableWGCNAThreads(nThreads=12)
  allowWGCNAThreads()
  options(stringsAsFactors = FALSE)

  dataset = t(dataset)
  adj = abs(bicor(dataset, use = 'pairwise.complete.obs'))^7
  dissTOM = TOMdist(adj, TOMType = "signed")
    
  hierTOMa = hclust(as.dist(dissTOM), method='average')

  Gene_Modules = labels2colors(cutreeDynamic(hierTOMa, method='tree', cutHeight=0.99))
  Gene_Clusters = labels2colors(cutreeDynamic(hierTOMa, distM= dissTOM , cutHeight = 0.99,
                                               deepSplit=4, pamRespectsDendro = FALSE))
    
  options(repr.plot.width=5,repr.plot.height=5,repr.plot.res=200)
  plotDendroAndColors(hierTOMa,
                      colors = data.frame(Gene_Clusters),
                      dendroLabels = FALSE,
                      cex.axis = 1.2)

  return(Gene_Modules)

}


Mods = TOM_Matrices(df_exprs)

eigenetic_network = function(dataset, colorh1){

  enableWGCNAThreads(nThreads=32)
  allowWGCNAThreads()
  options(stringsAsFactors = FALSE)

  ADJ1 = abs(bicor(t(dataset), use = 'all.obs'))^7

  colors = unique(colorh1)
  Alldegrees1 = intramodularConnectivity(ADJ1, colorh1)

  datME = moduleEigengenes(t(dataset),colorh1)$eigengenes
  MET = orderMEs(cbind(datME))
  datKME = signedKME(t(dataset), datME, outputColumnName='')

  return(datKME)

}

eig = eigenetic_network(df_exprs, Mods)


enrichment = function(dataset, colorh1, datKME){
    
    library(clusterProfiler)
    library(dplyr)
    library(org.Hs.eg.db)
    library(enrichplot)
    
    intModules = table(colorh1)
    intModules = as.data.frame(intModules)
    intModules =intModules$colorh1
    intModules = as.character(intModules)
    
    dat = data.frame()
    dat_new = data.frame()
    
    colrs = c()
    newclors = c()
    
    dataset = t(dataset)
    
    for (color in intModules){
        
        color =  color
        FilterGenes = abs(subset(datKME, select=c(color))) > 0.85
        genes = dimnames(data.frame(dataset))[[2]][FilterGenes]
        
        dat = cbind.fill(dat, genes, fill = NA)
        colrs = append(color, colrs)
    
    }
    
    dat = dat[,seq(1,ncol(dat),2)]
    colnames(dat) = colrs
    dat = as.data.frame(dat)
    
    dat = dat[,!names(dat) %in% c("grey")]
    colrs = colnames(dat)
    
    for (j in 1:ncol(dat)){
        
        gene = dat[, j]
        
        if (all(is.na(gene)) == FALSE){
            
            eg = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
            genes = eg$ENTREZID
            dat_new = cbind.fill(dat_new, genes, fill = NA)
        }
    }
    
    dat_new = dat_new[,seq(1, ncol(dat_new), 2)]
    
    colnames(dat_new) = newclors    
    dat_new = as.data.frame(dat_new)
    
    names = paste0('C', 1:ncol(dat_new))
    colnames(dat_new) = names
        
    return(dat_new)


}

enrich = enrichment(df_exprs, Mods, eig)

clusterKEGG = function(dat_new){
    
    dat_new = dat_new[,!names(dat_new) %in% c("grey")]
    ck = compareCluster(geneCluster = dat_new, fun = enrichKEGG, 
                        pvalueCutoff = 0.05)
    
    return(ck)

}

options(repr.plot.width=8,repr.plot.height=7,repr.plot.res=250)
enrichplot::dotplot(clusterKEGG(enrich), showCategory=5)
print("Analysis Complete!")