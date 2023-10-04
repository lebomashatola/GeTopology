'''

Preprocess HTSEq count data from https://portal.gdc.cancer.gov

'''

import os
import pandas as pd

class Expression:

    def __init__(self, HTSEqCounts, name):
        """
        Initialise input variables including file path for all HTSeq count folders and phenotype name
        """ 
        self.HTSEqCounts = HTSEqCounts
        self.name = name

    def read_files_create_dataframe(self):
        """
        Open each folder and read HTSeq txt file, select protein-coding genes with TPM counts 
        :return:Gene Expression Matrix and write to compressed PICKLE file 
        """ 
        self.df = pd.DataFrame()
        
        for file in os.listdir(self.HTSEqCounts):
            
            self.filename = os.fsdecode(file)
            self.read_data = self.HTSEqCounts + '/' + self.filename
            
            try:
                self.files = os.listdir(self.read_data)
            except:
                continue
                    
            for i in self.files:
                
                if i.endswith(".tsv"):
                    
                    self.df_files = self.read_data + '/' + i
                    self.df_ = pd.read_csv(self.df_files, sep='\t', header=1)
                    self.df_ = self.df_.iloc[4:] 
                    self.df_ = self.df_[self.df_['gene_type'].str.match('protein_coding')]
                    self.df_ = self.df_.set_index('gene_name')
                    self.df_ = self.df_[['tpm']] #Alternative options include tpm unstranded; fpkm_unstranded 
                    self.df = pd.concat([self.df, self.df_], axis=1)
    
        self.df = self.df.T
        self.df = pd.DataFrame(self.df)
    
        self.output_directory = '/directory/' + self.name + '.pkl'
        self.df.to_pickle(self.output_directory, compression='infer', protocol=5, storage_options=None)

if __name__ == '__main__':

    resistance = Expression("/dir/Resistance", "Resistance")
    resistance.read_files_create_dataframe()

    sensitive = Expression("/dir/Sensitive", "Sensitive")
    sensitive.read_files_create_dataframe()
