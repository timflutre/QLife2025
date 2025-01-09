Simulations generated for the QLife workshop (march 2025)       

### Installation of the necessary packages

Creation of the conda environment        
$ conda create -n qlife_simu python=3.9.18              

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


### Usage

Run the simulator: $ python src/main.py                                 

A certain number of parameters can be given as argument in the command line, under the format -parameter value.                          
The parameters that can be changed are the following (default value and expected type in parenthesis):                               
- savedFolder (~/default; string): folder where the outputs of the simulator will be saved                          
- optim (0; float): optimum of the fitness function                                 
- varW (100; float): width of the fitness function                         
- G (10; integer): number of generations                      
- N (100; integer): sample size (forward in-time selection)                       
- Npop (10000; integer): population size (coalescent burn-in)                        
- h2 (0.5; float): heritability                         
- nTrait (1; integer): number of characters                             
- nChr (1; integer): number of chromosomes                               
- Lchr (1000; integer): number of sites per chromosome                       
- rho (10e-8; float): recombination rate (per site)                         

Eg of usage with parameters: $ python src/main.py -optim 0 -varW 100                                  



