# genespots

genespots is an R script for identifying genomic hotspots—regions significantly enriched in genes of interest (GOIs)—using gene annotations from a GTF file and a provided list of GOIs (e.g., differentially expressed genes or gene clusters).


## Workflow overview
- Maps GOIs to their genomic coordinates using the provided GTF file.
- Scans each chromosome using a sliding window (default: 250 kb) to estimate gene density.
- Identifies peaks in gene density as candidate hotspots.
- Assesses statistical significance by generating a null distribution via permutation of GOI labels (default: 100 permutations).
- Computes empirical p-values and adjusts for multiple testing using FDR (default threshold: < 0.05).
- Outputs:
  - A table of significant hotspots.
  - A PDF plot showing observed densities, null densities, and hotspot locations.


## Installation

You can run this script locally in R or RStudio.
To clone the repository:
```bash
git clone https://github.com/ccarloscr/genespots.git
```


## Sample Inputs

The [Input/](Input/) folder contains example files:
- A GTF annotation file [dmel-all-r6.62.gtf](Input/dmel-all-r6.62.gtf)
- A table of GOIs [defs_degs.txt](Input/defs_degs.txt)
Make sure the GTF file is decompressed before use.



## Configuration

You can customize parameters by editing the script directly:
- **gtfdir**: path to the GTF annotation file.
- **clusterdir**: path to the table of genes of interest.
- **bandwidth**: size of the sliding window.
- **n_perm**: number of permutations for null distribution.
- **min_density_threshold**: fraction of max density used to filter hotspot peaks.
- **min_genes**: minimum number GOIs overlapping a sliding window to consider a valid hotspot.
- **filter_chromosomes**: logical flag to restrict plots to chromosomes with sig. hotspots.


## License

This project is licensed under the MIT License.
