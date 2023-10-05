import os 

while True:
    
    print(" GeTopology (v0.1) \n Welcome to Gene Preselection Analysis!")
    print("###################################")
    opt = int(input('Select option: \n 1. Analyse Microarray Data \n 2. Analyse RNA-Seq Data \n 3. Exit \ n : ')) 
    
    if opt == 1:
        os.system('Rscript Limma.r')

    if opt == 2:
        os.system('Rscript Deseq2.r')

    if opt == 3:
        break