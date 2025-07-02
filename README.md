# genespots

The R script identifies **genomic hotspots** - regions significantly enriched in genes of interest (GOIs) - using gene annotations from a GTF file and a provided list of GOIs (differentially expressed genes, gene clusters, etc). 

Script workflow:
- The script first maps GOIs to their genomic cordinates using the provided GTF file.
- Scans each chromosome using a **sliding window (default: 250 kb)** to estimate gene density.
- Peaks in gene density are identified as candidate hotspots.
- Statistical significance of these hotspots is tested generating a null distribution by randomly permuting the GOI labels across all genes (default: 100 permutations).
- Empirical p-values and FDRs (default: < 0.05) are computed for each hotspot.
- Ouputs: table of significant hotspots and plot showing observed densities, null densities and hotspot locations.



## Installation

This script can run in R or Rstudio in local.
Alternatively, to clone this repository use:
```bash
git clone https://github.com/ccarloscr/genespots.git
```


## Sample Inputs

The Input folder contains sample files for both the table of GOIs and the GTF file.
Note that the GTF file should be decompressed before use.



## Configuration

Multiple parameters can be customized by directly editing the script:
- **gtfdir**: path to the GTF annotation file.
- **clusterdir**: path to the table of genes of interest.
- **bandwidth**: size of the sliding window.
- **n_perm**: number of permutations for null distribution.
- **min_density_threshold**: fraction of max density used to filter hotspot peaks.
- **min_genes**: minimum number GOIs overlapping a sliding window to consider a valid hotspot.
- **filter_chromosomes**: logical flag to restrict plots to chromosomes with sig. hotspots.


