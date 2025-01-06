## neutral burn-in using msprime
## add selection via simupop

# $ cd ~/Desktop/qlife_2025/simu
# $ conda activate qlife
# $ python

#-------------------- import packages

from simuOpt import setOptions
setOptions(alleleType = 'lineage', quiet=True)

import msprime
import simuPOP as sim
import numpy as np
from scipy.stats import norm
import os

#-------------------- functions

## function to calculate the fitness based on the phenotypic value
class FisherGeometricalModel:
    def __init__(self, optW, varW):
        self.loc = optW
        self.scale = np.sqrt(varW)
        ## value of the fitness at the optimum 
        self.max_fitness = norm.pdf(self.loc, loc = self.loc, scale = self.scale)
    ## -------------------------------------
    ## -------------------------------------
    def fitness(self, phenotype):
        return norm.pdf(phenotype, loc = self.loc, scale = self.scale) / self.max_fitness

## function to set the generation of an individual 
## (when the ind was generated)
## will be used for the pedigree
## in particular when randomize the parents for a given generation
def setGen(gen):
    return gen[0] + 1

## initialize the genetic values, for the founders
def initGeneticValue(beta, pop):
    u = [sum([beta[locus] * (ind.allele(locus,0) + ind.allele(locus,1)) for locus in range(pop.totNumLoci())]) for ind in pop.individuals()]
    return list(u)

## calculate the genetic values when the population evolve
def setGeneticValue(pop, off, dad, mom, param):
    varE = param[0]
    beta = param[1]
    genetVal = sum([beta[locus] * (off.allele(locus,0) + off.allele(locus,1)) for locus in range(pop.totNumLoci())])
    off.setInfo(genetVal, pop.infoIdx('breedingValue'),)
    ## phenotype: add an environmental noise to the breeding value
    genetVal += np.random.normal(loc = 0, scale = np.sqrt(varE))
    off.setInfo(genetVal, pop.infoIdx('phenotype'),)
    return True

#-------------------- parameters

savedFolder = "/home/eliset/Desktop/qlife_2025/20250106_0"
if not os.path.exists(savedFolder):
    os.mkdir(savedFolder)
    os.mkdir(savedFolder + "/freq")
    os.mkdir(savedFolder + "/count")
    os.mkdir(savedFolder + "/pedigree")
    os.mkdir(savedFolder + "/genotype")
    os.mkdir(savedFolder + "/fit")
    os.mkdir(savedFolder + "/vcf")
    os.mkdir(savedFolder + "/beta")




optim = 5
varW = 1
h2 = 0.5
G = 10
rep = 1
N = 100
L = 10000 ## L = sequence length
rho = 10/L
mu = 10/L ## mutation rate



#-------------------- neutral burn-in

## generate the ancestry of N diploid individuals (i.e. 2N sequences)
## recombination rate: per-site value
ts = msprime.sim_ancestry(samples = N, sequence_length = L, population_size = 10000, recombination_rate = rho, random_seed = 1234)

## add (neutral) mutations to a tree sequence
mut_model = msprime.BinaryMutationModel() ## binary model to only have biallelic variants (0 or 1, 0 being the ancestral allele)
mts = msprime.sim_mutations(ts, rate = mu, random_seed = 5678, model = mut_model)

## print the variants
#for var in mts.variants():
#    print(var.site.position, var.alleles, var.genotypes, sep="\t")

## save the position in a list
pos = [var.site.position for var in mts.variants()]
## can then use it if we want to take into account the level of recombination between two variants based on the fact that their distance (in bp) will be different


## save the variants in a VCF file
with open(savedFolder + "/vcf/vcf" + str(rep) + ".txt", "w") as vcf_file:
    mts.write_vcf(vcf_file, allow_position_zero = True)

## export the genotype matrix as a numpy array
tmp = mts.genotype_matrix()
## format: tmp = [variant1, ..., variantLp]
## with variant = [ind1 chr1, ind1, chr2, ..., indN chr1, indN chr2]

Lp = tmp.shape[0]

## want the genotype matrix in this format:
## founders = [ind1, ..., indN]
## with ind = [variant1 chr1, ..., variantLp chr1, variant1 chr1, ..., variantLp chr2]
founders = [[] for i in range(N)]
for ind in range(N):
    founders[ind] = [[] for i in range(2*Lp)]
    for locus in range(Lp):
        founders[ind][2*locus:(2*locus+2)] = tmp[locus][2*ind:(2*ind+2)]


#-------------------- forward-in-time simulations, with selection

## consider that all the variants can be QTLs
## with a given distribution of effect (normal distribution?)
beta = list(norm.rvs(size = Lp, loc = 0, scale = np.sqrt(0.001)))
## ideally, the betas would be saved in the vcf, to get all the information about the variants at the same place
with open(savedFolder + "/beta/beta" + str(rep) + ".txt", "w") as out:
    for i in beta:
        out.write(str(i) + "\n")

## initialize the population
## with the same population size as the one defined by msprime
## and use the number of variants for the number of loci
## could also add the monomorphic loci, but they wouldn't add any information...
pop = sim.Population(size = N, loci = Lp, ancGen = -1, infoFields = ['ind_id', 'father_id', 'mother_id', 'gen', 'phenotype', 'breedingValue', 'fitness', 'avgFitness'])
sim.IdTagger().reset(1) ## just to make sure IDs starts from 1
sim.tagID(pop)
pop.setIndInfo(0, 'gen')

## add the genotypes defined by msprime
## check that the founder population that we get after initialization is the same as the one we got from msprime
founders2 = []
for i,ind in enumerate(pop.individuals()):
    ind.setGenotype(founders[i])
    founders2.append([ind.info('ind_id'), founders[i]])

selModel = FisherGeometricalModel(optW = optim, varW = varW)

## initialize the genetic values
geneticValue = initGeneticValue(beta, pop)
sim.initInfo(pop, geneticValue, infoFields = 'breedingValue')

## initialize the variances in the population
## ideally, would fix varP and h2 and determine varA and varE
## need to draw the betas in a given distribution to have the wanted varA
## i.e. have varP and h2 -> determine varA -> determine the distribution for betas -> and finally calculate the breeding values
varA = np.var(geneticValue)
varE = varA*(1-h2)/h2
varP = varE + varA

## add an environmental noise to get the phenotype
envNoise = list(norm.rvs(size = N, loc = 0, scale = np.sqrt(varE)))
phenotype = [sum(x) for x in zip(geneticValue, envNoise)]
sim.initInfo(pop, phenotype, infoFields = 'phenotype')

## add varA and N as parameters of the population
## varA will decrease in time due to drift as varA = varA(1-1/2N)
pop.dvars().varA = varA
pop.dvars().N = N


## evolution of the population
pop.evolve(
    initOps = [
        sim.InitLineage(mode = sim.FROM_INFO_SIGNED), ## save the lineage of each allele
        sim.InitSex(),
        sim.PedigreeTagger(outputFields = ['gen', 'phenotype', 'breedingValue'], output='>>' + savedFolder + '/pedigree/pedigree' + str(rep) + '.txt'), ## will need the first generation for the randomization
        sim.Stat(alleleFreq = sim.ALL_AVAIL, popSize = True),
        sim.PyEval(r"'-1\t%d\t' % subPopSize[0]", step = 1, output = '>>' + savedFolder + '/freq/freq' + str(rep) + ".txt"),
        sim.PyEval(r"'\t'.join(['%.4f' % alleleFreq[x][1] for x in range(len(alleleFreq))]) + '\n'", step = 1, output = '>>' + savedFolder + '/freq/freq' + str(rep) + ".txt"),
        sim.PyEval(r"'-1\t%d\t' % subPopSize[0]", step = 1, output = '>>' + savedFolder + '/count/count' + str(rep) + ".txt"),
        sim.PyEval(r"'\t'.join(['%.4f' % alleleNum[x][1] for x in range(len(alleleNum))]) + '\n'", step = 1, output = '>>' + savedFolder + '/count/count' + str(rep) + ".txt"),
    ],
    preOps = [
        sim.PySelector(func = selModel.fitness),
        ## save the fitness of the parents (offsprings of the generation t - 1)
        ## need to do it now, after having evaluated their fitness via PySelector
        ## if were to do it in postOps with avgFitness: evaluate the average fitness of the PARENTS (not what we want!!!)
        ## for the last generation (offspring), apply the selector outside of pop.evolve to get their fitness
        sim.PyEval(r"'%d\t' %gen", step = 1, output = '>>' + savedFolder + '/fit/fit' + str(rep) + ".txt"),
        sim.PyEval(r"'\t'.join(['%.4f' % x for x in fitness]) + '\n'", stmts = "fitness = pop.indInfo('fitness')", exposePop = 'pop', step = 1, output = '>>' + savedFolder + '/fit/fit' + str(rep) + ".txt"),
    ],
    matingScheme = sim.HomoMating(
        chooser = sim.RandomParentsChooser(True, selectionField = 'fitness'),
        generator = sim.OffspringGenerator(
            [
                sim.IdTagger(),
                sim.PyTagger(setGen), ## set the generation of the new ind to that of its father (that was incremented by 1 in the preOps step)
                sim.Recombinator(rates = rho), 
                sim.PyOperator(setGeneticValue, param = [varE, beta]), 
                sim.SummaryTagger(mode = sim.MEAN, infoFields = ['fitness', 'avgFitness']),  
                sim.PedigreeTagger(outputFields = ['gen', 'phenotype', 'breedingValue'], output='>>' + savedFolder + '/pedigree/pedigree' + str(rep) + '.txt'),
            ]),
    ),
    postOps = [
        sim.Stat(alleleFreq = sim.ALL_AVAIL, popSize = True),
        sim.PyEval(r"'%d\t%d\t' % (gen,subPopSize[0])", step = 1, output = '>>' + savedFolder + '/freq/freq' + str(rep) + ".txt"),
        sim.PyEval(r"'\t'.join(['%.4f' % alleleFreq[x][1] for x in range(len(alleleFreq))]) + '\n'", step = 1, output = '>>' + savedFolder + '/freq/freq' + str(rep) + ".txt"),
        sim.PyEval(r"'%d\t%d\t' % (gen,subPopSize[0])", step = 1, output = '>>' + savedFolder + '/count/count' + str(rep) + ".txt"),
        sim.PyEval(r"'\t'.join(['%.4f' % alleleNum[x][1] for x in range(len(alleleNum))]) + '\n'", step = 1, output = '>>' + savedFolder + '/count/count' + str(rep) + ".txt"),
        ## the additive variance decrease at each generation due to the selection
        ## as a function of the population size
        sim.PyExec("varA = varA*(1 - 0.5/N)"), 
    ],
    gen = 10
)

## apply the selector to the offsprings of the last generation to get their fitness
lastFit = []
for ind in pop.individuals():
    lastFit.append(selModel.fitness(ind.phenotype))

with open(savedFolder + '/fit/fit' + str(rep) + '.txt', 'a') as f:
    f.write(str(G+1) + '\t')
    f.write('\t'.join(['%.4f' % x for x in lastFit]) + '\n')

## dump the whole population (all generations) in a file
## genotypes included, as well as the pedigree informations
## for now do that, but will want to save the genotypes only in a specific file
sim.dump(pop, structure = False, ancGens = sim.ALL_AVAIL, output = ">>" + savedFolder + "/genotype/genotype_" + str(rep) + ".txt")

## save the origins (for the lineage tracing)
for ind in pop.individuals():
    lineage = ind.lineage()
    N = len(lineage) // 2
    with open(savedFolder + '/genotype/origins_' + str(rep) + '.txt', 'a') as f:
        towrite = lineage[:N]
        f.write(str(ind.ind_id) + '\t1\t' + '\t'.join(['%d' % x for x in towrite]) + '\n')
        towrite = lineage[N:]
        f.write(str(ind.ind_id) + '\t2\t' + '\t'.join(['%d' % x for x in towrite]) + '\n')
















