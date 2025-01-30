Simulations generated for the QLife workshop (march 2025)      

For the simulations, we use msprime (coalescent simulations) for the neutral burn-in phase and generate coalescence trees.                                     
Mutations are then overlayed to these trees, to generate the genotypes of the founders (initial population).                                      
We then use a forward-in-time simulator (simuPOP) to do multiple round of selection, using the previously generate founders as base population.                                    

Individual are selected based on their probability to have offsprings, i.e their fitness.                                                         
The fitness of an individual is calculated using its phenotype, via an exponential fitness function: exp(-(phenotype - optimum) ** 2 / varW).                                   
This function means that the farther the phenotype of an individual is from the optimum phenotype, the lower the fitness. This decrease of the fitness also depends on the width varW of the function: the smaller varW is, the stronger the selection will be.                   

The phenotype P is calculated via adding a genetic contribution G to an environmental contribution E.                             
The genetic contribution (or breeding value) is the part of the phenotype that will be transmitted to the next generation while the environmental contribution is a noise drawn from a normal distribution.                               
The heritability $h^2$ will determine the part of the phenotype that is due to the genetic contribution (if $h^2 = 1$, then the phenotype is equal to the breeding value; on the contrary, if $h^2 = 0$, then the phenotype is equal to the environmental noise).
Finally, the genetic contribution is obtained by adding the genotypes x at all SNPs, weighted by the SNP effects $\beta$.                                    
To summarize: 
- $G = \sum_{i = SNP} x_i \beta_i$
- $E \sim N(0, \sigma^{2}_{E})$
- $P = G + E$
- $h^2 = \sigma^{2}_{G} / \sigma^{2}_{P} = \sigma^{2}_{G} / (\sigma^{2}_{G} + \sigma^{2}_{E})$                                

### Installation of the necessary packages

Creation of the conda environment, via a yaml configuration file.  
```
conda env create -f simu_env.yml                           
```

```
conda activate qlife_simu                                                                  
```

### Usage

Run the simulator: 
```
python src/main.py                                 
```

A certain number of parameters can be given as argument in the command line, under the format -parameter value.                          
The parameters that can be changed are the following (default value and expected type in parenthesis):                               
- savedFolder (default; string): folder where the outputs of the simulator will be saved                          
- optim (0; float): optimum of the fitness function                                 
- varW (100; float): width of the fitness function                         
- G (10; integer): number of generations for the forward-in-time step                      
- N (100; integer): sample size (forward in-time selection)                       
- Npop (10000; integer): population size (coalescent burn-in)                        
- h2 ([0.5]; float): heritability (need to have more than one element in nTrait > 1; only the first trait will be under selection; full pleiotropy between the trait; cannot be <= 0)                         
- nTrait (1; integer): number of traits                             
- nChr (1; integer): number of chromosomes                               
- Lchr (1000; integer): number of sites per chromosome                       
- rho (1e-7; float): recombination rate (per site)
- mu (1e-4; float): mutation rate (per site)
- proportionQTL (1.0; float): proportion of the SNPs that are QTLs (need to be between 0 and 1)
- varEffect (1.0; float): variance of the distribution of the QTL effects
- corTrait (0.5; float): correlation between the effect of two trait                         

Some remarks on the parameters:
- the heritability must be > 0
- the number of heritabilities to give must be equal to nTrait
- if nTrait > 1, then only the first trait will be under selection
- we consider full pleiotropy between the traits: if a QTL has an effect for one trait, it will have an effect for the other traits (similarly, a neutral site will be neutral for all traits)
- if we want no LD between two markers, put them in two different chromosomes
- the parameter is only relevant when nTrait > 1. Note that, if there are more than 2 traits and do not want the same correlation between the effects of the different trait, need to directly modify the variance-covariance matrix in the script
- if the mutation rate is too low, we may have less variants (polymorphic sites) than Lchr, after the overlay of mutations on the tree (coalescent burn-in phase). In this case, we have added monomorphic (fixed) sites to have the total number of sites equal to Lchr per chromosome. Note that these sites have alternative alleles = NA

Eg of usage with parameters (used to generate the example data): 
```
python QLife2025/simulations/src/main.py -savedFolder 'example_data_1' -optim 5 -varW 1 -nTrait 2 -h2 0.5,0.9 -nChr 5 -Lchr 200 -rho 0.005
```

Some more examples of simulation of specific situations:                                     
- no selection 
```
python QLife2025/simulations/src/main.py -savedFolder 'example_data_1' -optim 0 -varW 100  
```

- selection
```
python QLife2025/simulations/src/main.py -savedFolder 'example_data_2' -optim 5 -varW 1  
```

- no linkage (free recombination): put one site per chromosome
```
python QLife2025/simulations/src/main.py -savedFolder 'example_data_3' -nChr 1000 -Lchr 1  
```

- two traits (each with their own heritability)
```
python QLife2025/simulations/src/main.py -savedFolder 'example_data_4' -nTrait 2 -h2 0.5,0.9
```

- only 10% of the sites are QTLs
```
python QLife2025/simulations/src/main.py -savedFolder 'example_data_5' -proportionQTL 0.1
```



### Outputs     

List of the outputs of the simulator:                
- log.txt (sep = '\t'): list of the parameters used for the simulation
- SNP_INFO.txt (sep = '\t'): information about the SNPs (SNP ID, chromosome, genetic position, reference allele, alternative allele, effect on the first trait, effect on the other traits - if more than 1, 1 column per trait)         
the alternative allele can be NA if the site is monomorphic (which can happen depending on the mutation rate chosen for the overlay of mutations after the coalescent burn-in)
- pedigree.txt (sep = '\t'): pedigree for the forward-in-time step - all generations included (ID of the individual, father ID, mother ID, sex, affection, generation, phenotype trait 1, breeding value trait 1, phenotype and breeding value for the other traits if more than 1, fitness)
- genotype.txt: genotypes (0,1,2) of the individuals of the forward-in-time phase; the first column is the ind_id (dim = N*(G+1) x (L+1), N = population size of the forward-in-time step, G = number of generations of this step) and the first line is the SNP ID
- IBD.txt: origin of the alleles of the individuals of the forward-in-time phase (via the id of the founder - i.e. g = 0 - with a sign (+ or -) to indicate the chromosome of origin); the first column is the ind_id and the second column is the homologous chromosome (dim = 2*N*(G+1) x (L+2)); the first line is the SNP ID 


