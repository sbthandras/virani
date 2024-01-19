### Scope  

The scope of this repository is to provide and ANI (Average Nucleotide Identity) measurement framework,
that not only takes into account the length difference between genomes but also works with alignments 
as measurement units. More exactly, this tool tries to find the best alignment for each part of the genome,
and if two alignments overlap, it chooses the better one over the one with less percentage identity, only on the 
overlapping section, but keeps the worse on sections when it is the only alignment.

### Dependencies  

- ncbi-blast
- R, specifically with these libraries:
  - tidyverse
  - foreach
  - doParallel
  - grid
  - dplyr

### Running

1. Clone the repo
   ```sh
   git clone https://github.com/sbthandras/virani
   ```
2. create the **genomesdir** and **results*
3. In the **genomesdir** create an **input_folder** to place your genomes in .fna nucleotide format
4. Rscript scripts/clean.R **threads** **genomesdir/input_folder** **results/output_folder**
