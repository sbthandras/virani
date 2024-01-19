### Scope  

The scope of this repository is to provide and ANI (Average Nucleotide Identity) measurement framework,
that not only takes into account the length difference between genomes but also works with alignments 
as measurement units. More exactly, this tool tries to find the best alignment for each part of the genome,
and if two alignments overlap, it chooses the better one over the one with less percentage identity, only on the 
overlapping section, but keeps the worse on sections when it is the only alignment.

### Dependencies  

- ncbi-blast
- R, specifically with these libraries
  - tidyverse
  - foreach
  - doParallel
  - grid

### Installation

_Below is an example of how you can instruct your audience on installing and setting up your app. This template doesn't rely on any external dependencies or services._

1. Get a free API Key at [https://example.com](https://example.com)
2. Clone the repo
   ```sh
   git clone https://github.com/your_username_/Project-Name.git
   ```
3. Install NPM packages
   ```sh
   npm install
   ```
4. Enter your API in `config.js`
   ```js
   const API_KEY = 'ENTER YOUR API';
   ```
