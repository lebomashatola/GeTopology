

readexprs = function(path, lbl){
    #read gene expression matrices
    df = read.csv(path , sep=',')
    rownames(df) = df$X
    df = df[, -1]
    
    names = paste0(lbl, 1:ncol(df))
    colnames(df) = names
    
    return(df)
}

TOM_Matrices = function(dataset, dir, name){
    #Compute Hub Genes and Preselect Dataset for the Construction of the Signed-TOM matrix
    library(WGCNA)
    enableWGCNAThreads(nThreads=12)
    allowWGCNAThreads()
    options(stringsAsFactors = FALSE)
    
    dat = c()
    
    adj = abs(bicor(t(dataset), use = 'pairwise.complete.obs'))^7
    dissTOM = TOMdist(adj, TOMType = "signed")
    hierTOMa = hclust(as.dist(dissTOM), method='average')
    
    Gene_Modules = labels2colors(cutreeDynamic(hierTOMa, method='tree', cutHeight=0.99))
    colors = unique(Gene_Modules)
    
    datME = moduleEigengenes(t(dataset),Gene_Modules)$eigengenes
    datKME = signedKME(t(dataset), datME, outputColumnName='')
        
    intModules = table(Gene_Modules)
    intModules = as.data.frame(intModules)
    intModules =intModules$Gene_Modules
    intModules = as.character(intModules)
    
    for (color in intModules){
        
        FilterGenes = abs(subset(datKME, select=c(color))) > 0.85
        genes = dimnames(data.frame(t(dataset)))[[2]][FilterGenes]
        
        dat = append(dat, genes)
    
    }
    
    dat = unique(dat)
    dataset = dataset[dat, ]
    
    adj = abs(bicor(t(dataset), use = 'pairwise.complete.obs'))^power
    dissTOM = TOMdist(adj, TOMType = "signed")
    
    write.csv(dissTOM, paste(dir, '/', name, '_signed_TOM',sep=''))

}

print("Constructing the Signed-TOM")
print("")
path_dir_in = readline("Enter input file path: ")
path_name = readline("Enter phenotype of patients: ")

df = read_exprs(path_dir_in, path_name)
path_dir_out = readline("Enter output file path: ")
TOM_Matrices(resistance, path_dir_out, path_name)
print("Process Complete!")