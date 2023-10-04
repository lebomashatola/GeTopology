# GeTopology - Bimodal Phenotype Prediction

## Description

GeTopology processes RNA-Seq files (.fastq), microarray files (.CEL) and tissue slide images (.SVS), performs shape encoding of gene expression data from TDA, and performs predictive learning tasks using a bimodal approach with convolutional and multilayer perceptronic neural networks to predict phenotypes. 

## Getting Started

### Programming Languages Required

* Python >= 3.9 
* R >= 4.0

### Required Software

* Linux/MacOS
* FastQC 0.12.0 [https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* Trimmomatics 0.36 [http://www.usadellab.org/cms/?page=trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* Star 2.7.11a [https://github.com/alexdobin/STAR](https://github.com/alexdobin/STAR)
* HTSeq 2.0.3 [https://htseq.readthedocs.io/en/master/install.html](https://htseq.readthedocs.io/en/master/install.html)

### Installing Dependencies


```
python3 -m pip install -r requirements.txt
```

```
Rscript requirements.r
```

### Executing program
## Data Processing
* Image Processing - Image Segmentation and Tiling
```
code blocks for commands
```
* Microarray Processing - RMA normalization of CEL files, forming a Gene Expression Matrix
```
code blocks for commands
```

* RNA-Seq Processing - Combining HTSeq counts forming a Gene Expression Matrix
```
code blocks for commands
```
## Gene Preselection

* Limma - WGCNA
```
code blocks for commands
```
* DESeq2 - WGCNA 
```
code blocks for commands
```

## Data Encoding

* Topological Data Analysis 
```
code blocks for commands
```

## Predictive Learning 

* Topological Summaries Classification - Multilayer Perceptron
```
code blocks for commands
```

* Image Classification - Convolutional Neural Networks
```
code blocks for commands
```

* Bimodal Classification (MLP & CNN)
```
code blocks for commands
```

## Authors

Contributors names and contact info

Lebohang Mashatola  [681452@students.wits.ac.za](681452@students.wits.ac.za)

Mandeep kaur [Mandeep.Kaur@wits.ac.za](Mandeep.Kaur@wits.ac.za)

## Version History

* 0.1
    * Initial Release


