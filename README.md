# MoChAP-Modification-of-ChIP-Analysis-Pipeline
MoChAP is a one-step ChIP analysis pipeline that begins by aligning raw sequencing reads using Smith-Waterman "local alignment." The pipeline then performs motif analysis on biological sequences, leveraging data from a SAM file to assess the presence of a conserved motif. Visualization of motif locations, retrieval of sequences from an online source, comparison of motif locations, identification of common locations, and calculation and visualization of a motif matrix are integral components.

To enhance accuracy, MoChAP incorporates a deep convolutional neural network model trained on curated ChIP-Seq peak data. This model effectively distinguishes bona fide peaks from false ones in the visualized data. Additionally, MoChAP utilizes MetaScope for the visualization of ChIP-seq results, providing a comprehensive and efficient analysis of ChIP data.



### Python packages
Install with pip
```
Required Python packages:
numpy
scipy
tensorflow
Keras
pandas
matplotlib
os
requests
```

### Installation
Install my-project with npm
```
 Colab or jupyter notebook 
     Using Colab in local runtime : https://research.google.com/colaboratory/local-runtimes.html 
Local computer : Download .py file ( python 3 )
```
    
## Pipeline
![Workflow](https://drive.google.com/uc?id=1LL8q6Lu6Es0ECQH4OL5IlKwI_IffXrR9)




