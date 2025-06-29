# penrose-dimer-monte-carlo
Monte Carlo code to calculate correlators of a quantum monomer-dimer model on the quasicrystalline Penrose tiling. The code calculates dimer-dimer, vison-vison correlators as well as open Wilson lines and closed Wilson loops corresponding to monomers and visons. The results from this code were used in the preprint: https://arxiv.org/abs/2503.15588

# Getting started

First, clone this repository and extract the data.zip file. 

There are four parts to the code:

1. Generate the Penrose tiling using lattice_generation.ipynb, and save all the requried information about the graph such as the location of the vertices and the adjacency matrix. This data has been already generated and stored in ./data/tiling_input

2. Once the tiling is generated, we generate many vison and monomer loops/paths and use them to do Monte Carlo simulation using main-monte-carlo.ipynb. The resulting correlators are stored in ./data/R2_mc and ./data/R3_mc for the two regions used in the paper.

3. Next, we use the data from these folders to analyze the data (perform path averaging and other post-processing) using data-analysis.ipynb. The analyzed data is stored in ./data/plot_data

4. Finally, ./data/plot_data is used by ./plot.ipynb to generate the plots included in the paper, which are saved in ./images/paper/.
