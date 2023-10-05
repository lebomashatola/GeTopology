
readexprs <- function(path, lbl){
    #read gene expression matrices
    df = read.csv(path , sep=',')
    rownames(df) = df$X
    df = df[, -1]
    
    names = paste0(lbl, 1:ncol(df))
    colnames(df) = names
    
    return(df)
}

print("Gene Preselection")
print("")
path_r = readlines("Enter file path for resistance patients: ")
path_s = readlines("Enter file path for sensitive patients: ")

R = readexprs(path_r, "R")
S = readexprs(path_s, "S")

limma = function(sensitive, resistance){
    #Perform differential gene expression analysis using Limma
    library(limma)
    library(data.table)
    library(dplyr)
    library(EnhancedVolcano)
    
    design = model.matrix(~ 0 + factor(c(rep(1, ncol(resistance)), rep(2, ncol(sensitive)))))
    colnames(design) = c('R', 'S')
    
    dataset = cbind(resistance, sensitive) 
    dataset = normalizeBetweenArrays(dataset, method="quantile")
    
    fit = lmFit(dataset, design, method='ls') 
    
    contr = makeContrasts(R - S, levels = colnames(coef(fit)))
    tmp = contrasts.fit(fit, contr)
    
    fit2 = eBayes(tmp, robust=TRUE) 
    
    plotSA(fit2, main="Gene-level")
        
    results = decideTests(fit2, method='hierarchical', adjust.method='BH', p.value=0.05)
    sig = topTable(fit2, n=Inf, adjust='BH', coef=1, sort.by='P', p.value=0.05)
    
    up = subset(sig, logFC > 1)
    down = subset(sig, logFC < -1)
    
    up = resistance[rownames(up) , ]
    down = resistance[rownames(down) , ]
    all = rbind(up, down)

    return(all)

}

deg = limma(R,S)

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

mods = TOM_Matrices(deg)

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

eig = eigenetic_network(deg, mods)

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

enrich = enrichment(deg, mods, eig)

clusterKEGG = function(dat_new){
    
    dat_new = dat_new[,!names(dat_new) %in% c("grey")]
    ck = compareCluster(geneCluster = dat_new, fun = enrichKEGG, 
                        pvalueCutoff = 0.05)
    
    return(ck)

}


options(repr.plot.width=8,repr.plot.height=7,repr.plot.res=250)
enrichplot::dotplot(clusterKEGG(enrich), showCategory=5)