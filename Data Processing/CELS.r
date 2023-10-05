

library(biomaRt)

preprocess = function(file_in, lbl){
    
    library(oligo)
    library(gcrma)
    library(data.table)
    library(biomaRt)
    
    read = list.celfiles(file_in, full.name=TRUE)
    read = read.celfiles(read)
    read = oligo::rma(read)
    read = exprs(read)
    
    probeIDs = rownames(read)
    read = as.data.frame(read)
    
    read = setDT(read, keep.rownames = TRUE)[]
    colnames(read)[1] = 'affy_hg_u133_plus_2'
    
    mart = useMart(biomart='ENSEMBL_MART_ENSEMBL',
                dataset = 'hsapiens_gene_ensembl',
                host = 'https://useast.ensembl.org')
    
    df = getBM(attributes= c('affy_hg_u133_plus_2', 'hgnc_symbol'),
               filters = 'affy_hg_u133_plus_2',
               values = probeIDs,
               mart = mart)
    
    df_new = merge(df, read, by='affy_hg_u133_plus_2')
    df_new = df_new[,-1]
    
    df_new = df_new[!duplicated(df_new$hgnc_symbol), ]
    rownames(df_new) = df_new$hgnc_symbol
    
    df_new = df_new[,-1]
    write.csv(df_new, paste(file_in, '/', lbl, '.csv', sep=''))
    
}

file_in = readlines("Enter input file path for microarray CEL files: ")
lbl = readlines("Enter label for samples: ")

preprocess(file_in, lbl)