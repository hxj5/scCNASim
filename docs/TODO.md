# TODO List

## Input


## Output


## Implementation
### cs
- Save memory - split adata by features and save splitted adata into
  files for each thread (batch).
- Write a wrapper function for cs module to make it an independent tool for
  count(RDR)-based CNA simulatioin, to take CNA profile and clone annotation
  files as input.
- Add an option about `total library size` to enable simulation of various
  read depth.
  Refer to the `total number of reads` in scDesign2.


## Docs


## Tests
- Prepare notebooks & scripts to test each module.


## Scripts
- Write python version of benchmarking pipeline for CNA detection, mainly for
  inferCNV and Numbat.


## Discussion
### How to process some special UMIs in feature counting?
1. Within one UMI, reads overlap different features, e.g., due to errors in
   UMI collapsing, or features overlapping each other.

2. Within one UMI, some reads overlap one feature while others do not overlap
   any features.

3. Within one UMI, all reads overlap one feature, while some pass the
   --minINCLUDE filtering and others do not.

### The noise in the reference cells in inferCNV heatmap.
Noise: both signals of copy gain and copy loss shown in a strip of genes.
(1) process (e.g., remove) the special UMIs/reads listed above that the noise
    arises from.
(2) process (e.g., remove) the special features (e.g., overlapping features) 
    that the noise arises from.
