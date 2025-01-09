Simulations generated for the QLife workshop (march 2025)       

### Installation of the necessary packages

Creation of the conda environment, via a yaml configuration file.  
$ conda env create -f simu_env.yml                           

For the simulations, we use msprime (coalescent simulations) for the neutral burn-in phase and generate the initial population.                   
We then use a forward-in-time simulator (simuPOP) to add selection (stabilizing selection)                                                
$ conda activate qlife_simu                                                                  

### Usage

Run the simulator: $ python src/main.py                                 

A certain number of parameters can be given as argument in the command line, under the format -parameter value.                          
The parameters that can be changed are the following (default value and expected type in parenthesis):                               
- savedFolder (/home/eliset/Desktop/qlife_2025/default; string): folder where the outputs of the simulator will be saved                          
- optim (0; float): optimum of the fitness function                                 
- varW (100; float): width of the fitness function                         
- G (10; integer): number of generations for the forward-in-time step                      
- N (100; integer): sample size (forward in-time selection)                       
- Npop (10000; integer): population size (coalescent burn-in)                        
- h2 (0.5; float): heritability                         
- nTrait (1; integer): number of characters                             
- nChr (1; integer): number of chromosomes                               
- Lchr (1000; integer): number of sites per chromosome                       
- rho (10e-8; float): recombination rate (per site)                         

Eg of usage with parameters (used to generate the example data): $ python QLife2025/simulations/src/main.py -savedFolder '/home/eliset/Desktop/qlife_2025/example_data' -optim 5 -varW 1 -nTrait 2 -nChr 5 -Lchr 200 -rho 0.005                              


### Outputs     

List of the outputs of the simulator:                
- log.txt: list of the parameters used for the simulation
- beta.txt: list of the different marker effects (dim = L x nTrait, L = total number of sites and nTrait = number of traits)
- pedigree.txt: pedigree for the forward-in-time step - all generations included (columns = ind_id father_id mother_id sex affection generation phenotype breeding_value)
- genotype.txt: genotypes (0,1,2) of the individuals of the forward-in-time phase; the first column is the ind_id (dim = N*(G+1) x (L+1), N = population size of the forward-in-time step, G = number of generations of this step)
- IBD.txt: origin of the alleles of the individuals of the forward-in-time phase (via the id of the founder - i.e. g = 0 - with a sign (+ or -) to indicate the chromosome of origin); the first column is the ind_id and the second column is the homologous chromosome (dim = 2*N*(G+1) x (L+2)) 


