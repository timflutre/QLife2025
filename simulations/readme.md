Simulations generated for the QLife workshop (march 2025)

Creation of the conda environment
$ conda create -n qlife python=3.9.18

For the simulations, we want to use msprime (coalescent simulations) for the neutral burn-in phase, and then use a forward-in-time simulator (simuPOP, SLiM) to add selection (stabilizing selection)
For now, did the installation in the slim environment (where slim was already installed)
$ conda activate qlife
$ conda install -c conda-forge simupop
$ conda install -c conda-forge msprime
$ conda install -c conda-forge scipy

Version of the installed packages
simupop 1.1.17
msprime 1.3.3
tskit 0.6.0
scipy 1.10.1




