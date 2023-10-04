# GeTopology - Bimodal Phenotype Prediction

## Description

GeTopology processes RNA-Seq files (.fastq), microarray files (.CEL) and tissue slide images (.SVS), performs shape encoding on gene expression data suing TDA, and performs a predictive learning task using a bimodal approach using convolutional and multilayer perceptron neural networks to predict phenotype. 

## Getting Started

### Programming Languages Required

* Python >= 3.9 
* R >= 4.0
* Unix OS 

### Installing Dependencies

```
python -m pip install -r requirements.txt
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
* Microarray Processing - RMA normalization CEL files, forming a Gene Expression Matrix
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

* Topological Summaries Classification - Multilayer Perceptron Layers
```
code blocks for commands
```

* Image Classification - Convolutional Neural Networks
```
code blocks for commands
```

* Bimodal Classification 
```
code blocks for commands
```

## Authors

Contributors names and contact info

Lebohang Mashatola  
[681452@students.wits.ac.za](681452@students.wits.ac.za)

## Version History

* 0.1
    * Initial Release


