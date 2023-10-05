'''
Paired-end RNA-Sequencing Processing, including:
- QC 
- Adaptor sequence trimming
- Sequence alignment 
- Aligned sequences QC
- Quantification of RNA-Seq counts  
'''

import os
import glob
import pandas as pd

class Process:

    def __init__(self, fastq, trimmed, genome, fasta, bam, bamQC, count):
        """
        Initialise the file paths to read and save outputs from 
        """
        self.fastq = fastq
        self.trimmed = trimmed
        self.genome = genome
        self.fasta = fasta
        self.bam = bam
        self.bamQC = bamQC
        self.count = count

    def fastqc(self):
        """
        Read fastq files and perform QC using FastQC-0.12.0
        :returns:FastQC report saved as HTML files
        """
        os.chdir(self.fastq)
        dirs_ = os.listdir(self.fastq)
        dirs_.pop(0)
        
        for i in dirs_:
            
            os.chdir(self.fastq + i)
            seq = [f for f in glob.glob("*.fastq.gz")]
            os.system('fastqc -t 12 ' + self.fastq + i + '/' + seq[0] + ' ' + self.fastq + i + '/' + seq[1])

    def trimmomatics(self):
        """
        Trim adapter sequences using Trimmomatics-0.39
        :return:outputs trimmed fastq files 
        """
        dirs_ = os.listdir(self.fastq)
        dirs_.pop(0)
        
        for i in dirs_:
            os.chdir(self.fastq + i)
            seq = [f for f in glob.glob("*.fastq.gz")]
            
            os.chdir('/Trimmomatic-0.39')
            os.system('java -jar trimmomatic-0.39.jar PE -phred33 ' + self.fastq + i + '/' + seq[0] + ' ' + self.fastq + i + '/' + seq[1] + ' ' +
                self.trimmed  + i + '/' + i + '_forward_paired.fastq.gz ' + self.trimmed + i + '/' + i + '_forward_unpaired.fastq.gz ' +
                self.trimmed + i + '/' + i + '_reverse_paired.fastq.gz ' + self.trimmed + i + '/' + i + '_reverse_unpaired.fastq.gz ' +
                'ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36')

    def finish(self):
        """
        Remove unpaired files generated after adptor sequence trimming
        :return:trimmed paired.fastq.gz files 
        """
        dirs_ = os.listdir(self.trimmed)
        dirs_.pop(0)
        
        for i in dirs_:
            os.system('rm -rf ' + self.trimmed + i + '/*_unpaired*')

    def build_star(self):
        """
        Assmeble STAR genome using reference hg38.fa file
        """
        os.system('STAR --runThreadN 12 --runMode genomeGenerate --genomeDir ' + self.genome + 
        ' --genomeFastaFiles ' + self.fasta) #Genome/hg38.fa

     def align_star(self):
        """
        Sequence alignment to reference genome
        :return:Bam files using STAR aligner
        """
        dirs_ = os.listdir(self.trimmed) #OR self.fastq
        dirs_.pop(0)
        
        for i in dirs_:
            os.chdir(self.trimmed + i)
            seq = [f for f in glob.glob("*.fastq.gz")]
            
            os.system('STAR --runThreadN 12 --genomeDir ' + self.genome + ' --readFilesIn ' + seq[0] + ' ' + seq[1] +
                  ' --outFileNamePrefix ' + self.bam + i + '/' + i +
                  ' --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat')   

    def align_QC(self):
         """
        Sequence alignment QC on BAM files
        :return:QC statistics and plots using samtools
        """
        dirs_ = os.listdir(self.trimmed)
        dirs_.pop(0)
        
        for i in dirs_:
            os.system('samtools stats ' + self.bam + '/' + i + '/' + i + '.bam > ' + self.bamQC + i + '/' + i + '_out.stats')
            os.system('plot-bamstats -p ' + self.bamQC + i + '/' + i + '_outQC ' + self.bamQC + i + '/' + i + '_out.stats')

    def count(self):
        """
        Quantification of RNA-Seq counts 
        :return:HTSEQ counts (.txt) files 
        """
        dirs_ = os.listdir(self.bam)
        dirs_.pop(0)
        
        for i in dirs_:
            os.system('htseq-count -f ' + self.bam + i + '~/gencode.v38.annotation.gff3 ' + self.count + i + '.txt')

if __name__ == '__main__':

    exp = Process("/dir/fastq", "/dir/trimmed", "/dir/genome", "/dir/fasta", "/dir/bam", "/dir/bamQC", "/dir/count")

    exp.fastqc()
    exp.trimmomatics()
    exp.finish()
    exp.build_star()
    exp.align_star()
    exp.align_QC()
    exp.count()

    while True:
        
        print(" GeTopology (v0.1) \n Welcome to RNA-Seq Processing!")
        print("###################################")
        opt = int(input('Select option: \n 1. FastQC \n 2. Trim Adaptor Seq \n 3. Run STAR \n 4. Run Alignment QC \n 5. Get Counts \n 6. Exit \n : '))

        if opt == 1:
            exp.fastqc()
            print('Process Complete! \n')

        if opt == 2:
            exp.trimmomatics()
            exp.finish()
            print('Process Complete! \n')
    
        if opt == 3:
            exp.build_star()
            exp.align_star()
            print('Process Complete! \n')

        if opt == 4:
            exp.align_QC()
            print('Process Complete! \n')

        if opt == 5:
            exp.count()
            print('Process Complete! \n')

        if opt == 6:
            break
