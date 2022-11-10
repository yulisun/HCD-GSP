# HCD-GSP
Graph Signal Processing for Heterogeneous Change Detection.

## Introduction

This paper provides a new strategy for the heterogeneous change detection (HCD) problem: solving HCD
from the perspective of graph signal processing (GSP). We construct a graph to represent the structure of each image,
and treat each image as a graph signal defined on the graph. In this way, we convert the HCD into a GSP problem: a
comparison of the responses of signals on systems defined on the graphs, which attempts to find structural differences and
signal differences due to the changes between heterogeneous images. 

In the previous version, this paper consists of two parts: part I, calculating the structural differences based on Vertex Domain Filtering (arXiv:2208.01881); 
part II, calculating the signal differences based on Spectral Domain Analysis (arXiv:2208.01905).

===================================================

## Datasets and Graph Cut algorithm

dataset#2 is download from Professor Michele Volpi's webpage at https://sites.google.com/site/michelevolpiresearch/home.

dataset#5 is download from Dr. Luigi Tommaso Luppino's webpage (https://sites.google.com/view/luppino/data) and it was downsampled to 875*500 as shown in our paper.

The graphCut algorithm is download from Professor Anton Osokin's webpage at https://github.com/aosokin/graphCutMex_BoykovKolmogorov.

If you use these resources, please cite their relevant papers.

===================================================

## Citation

If you use this code for your research, please cite our paper. Thank you!

@ARTICLE{9477152,
  author={Sun, Yuli and Lei, Lin and Guan, Dongdong and Kuang, Gangyao and Liu, Li},  
  journal={IEEE Transactions on Geoscience and Remote Sensing},   
  title={Graph Signal Processing for Heterogeneous Change Detection},   
  year={2022}}  
  
## Running

Unzip the Zip files (GC) and run the VDF and SDA demo files (tested in Matlab 2016a)! 

If you have any queries, please do not hesitate to contact me (sunyuli@mail.ustc.edu.cn).
