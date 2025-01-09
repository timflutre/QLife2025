#!/usr/bin/python

## neutral burn-in using msprime
## add selection via simupop

# $ cd ~/Desktop/qlife_2025/simu
# $ conda activate qlife
# $ python


#-------------------- functions

## function that read the arguments that were passed from the command line
## and put them in a dicttionary by option name
def getopts(argv, default):
    opts = {}
    while argv:
        if argv[0][0] == '-': ## find "-name value" pairs
            opts[argv[0][1:]] = argv[1] ## dict key is "-name arg"
            argv = argv[2:]
        else:
            argv = argv[1:]
    for i,j in default.items():
        if i not in opts: ## add the missing parameters and give them the default values
            opts[i] = default[i]
    return opts

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

## for a given generation, write the lineage in a file
## lineage = list of size 2*N*L
## -> separate by homologous chromosomes and by individuals
def saveLineage(pop, param):
    savedFolder = param[0]
    for ind in pop.individuals():
        lineage = ind.lineage()
        M = len(lineage) // 2
        with open(savedFolder + '/IBD.txt', 'a') as f:
            towrite = lineage[:M]
            f.write(str(ind.ind_id) + '\thom1\t' + '\t'.join(['%d' % x for x in towrite]) + '\n')
            towrite = lineage[M:]
            f.write(str(ind.ind_id) + '\thom2\t' + '\t'.join(['%d' % x for x in towrite]) + '\n')
    return True

#-------------------- parameters

if __name__ == '__main__':
    
    from simuOpt import setOptions
    setOptions(alleleType = 'lineage', quiet=True)

    import msprime
    import simuPOP as sim
    import numpy as np
    from scipy.stats import norm
    from scipy.stats import multivariate_normal
    import os
    import sys

    from simuPOP.utils import export
    from simuPOP.utils import Exporter

    default = {'savedFolder' : "~/default",
               'optim' : 5,
               'varW' : 1,
               'h2': 0.5,
               'G' : 10,
               'N' : 100,
               'Npop' : 10000,
               'nTrait' : 1,
               'nChr' : 1,
               'Lchr' : 1000,
               'rho' : 10e-8}

    parameters = getopts(sys.argv, default)

    SEED = sim.getRNG().seed()

    savedFolder = parameters['savedFolder']
    rep = 1

    optim = float(parameters['optim'])
    varW = float(parameters['varW'])

    G = int(parameters['G'])
    N = int(parameters['N'])
    Npop = int(parameters['Npop'])

    h2 = float(parameters['h2'])
    nTrait = int(parameters['nTrait'])
    nChr = int(parameters['nChr']) ## the chromosomes are independant
    Lchr = int(parameters['Lchr']) ## Lchr = sequence length (number of sites), per chromosome // consider that all chromosomes have the same number of sites
    rho = float(parameters['rho']) ## if want free recombination, there will be one QTL per chromosome and Lchr = 1
    mu = 10/Lchr

    varEffect = 1
    if nTrait > 1:
        covTrait = [[varEffect, 0.5], [0.5, varEffect]]

    if h2 <= 0:
        sys.exit('Null heritability is not allowed')

    if not os.path.exists(savedFolder):
        os.mkdir(savedFolder)


    with open(savedFolder + '/log.txt', 'a') as f:
        f.write('parameter\tvalue\n')
        f.write('seed\t%d\n' % SEED)
        for key, value in parameters.items():  
            f.write('%s\t%s\n' % (key, str(value)))

    #-------------------- neutral burn-in

    ## generate the ancestry of N diploid individuals (i.e. 2N sequences)
    ## recombination rate: per-site value
    ## each chromosomes are independant -> one replicate = one chromosome
    ts = msprime.sim_ancestry(samples = N, sequence_length = Lchr, population_size = Npop, recombination_rate = rho, num_replicates = nChr)

    ## add (neutral) mutations to a tree sequence
    mut_model = msprime.BinaryMutationModel() ## binary model to only have biallelic variants (0 or 1, 0 being the ancestral allele)
    genoList = np.empty(shape = 0)
    for i, tchr in enumerate(ts):
        mts = msprime.sim_mutations(tchr, rate = mu, model = mut_model)
        genoList = np.append(genoList, mts.genotype_matrix()) ## export the genotype matrix as a numpy array
        ## mts.genotype_matrix() format: [chr1, ..., nChr]
        ## with chr = array([variant1, ..., variant Lchr])
        ## with variant = [ind1 homologous1, ind1 hom2, ..., indN hom1, indN hom2]
        ## => genoList format = [chr1 loc1 ind1 hom1, chr1 loc1 ind1 hom2, chr1 loc1 ind2, chr1 loc2, chr2 ...]
        # print([var.genotypes.tolist() for var in mts.variants()])
        with open(savedFolder + '/log.txt', 'a') as f:
            f.write('number of variants on chromosome %d: %d \n' % (i, len([var.position for var in mts.variants()])))


    ## want the genotype matrix in this format:
    ## founders = [ind1, ..., indN]
    ## with ind = [variant1 hom1, ..., variantL hom1, variant1 hom2, ..., variantL hom2]
    L = nChr * Lchr
    founders = [[] for i in range(N)]
    for ind in range(N):
        founders[ind] = [[] for i in range(2*L)]
        ix1 = list(range(ind*2, 2*L*N, 2*N)) ## index of the variants in genoList for the first homologous chromosome
        ix2 = list(range(ind*2+1, 2*L*N, 2*N)) ## index of the variants in genoList for the second homologous chromosome
        for loc in range(L):
            founders[ind][loc] = genoList[ix1[loc]]
            founders[ind][loc+L] = genoList[ix2[loc]] ## allele on the homologous chromosome (diploids)


    #-------------------- forward-in-time simulations, with selection

    ## consider that all the variants can be QTLs
    ## with a given distribution of effect (normal distribution?)
    if nTrait == 1:
        beta = list(norm.rvs(size = L, loc = 0, scale = np.sqrt(varEffect)))
        with open(savedFolder + "/beta.txt", "w") as out:
            for i in beta:
                out.write(str(i) + "\n")

    if nTrait > 1:
        ## if there are more than one trait, use a multivariate normal to define correlated beta
        ## use the betas of the first trait (arbitrary decision) to calculate the trait that will be used to define the fitness in simuPOP
        ## to not calculate the breeding value and the phenotype of the other traits (can be calculated as posteriori as we have the betas)
        betaAll = multivariate_normal.rvs(size = L, mean = [0]*nTrait, cov = covTrait)
        with open(savedFolder + "/beta.txt", "w") as out:
            for i in betaAll:
                out.write('\t'.join(['%.4f' % x for x in list(i)]) + "\n")
        ## to define the trait associated in simupop, use the first trait
        beta = [i[0] for i in betaAll]

    ## for now the betas are drawn from a normal distribution, but will also need to include:
    ## - a given proportion are neutral with the others being major effects
    ## - multi-trait: correlation between effects; effects dranwn from a multi-normal distribution

    ## initialize the population
    ## with the same population size as the one defined by msprime
    ## and use the number of variants for the number of loci
    ## could also add the monomorphic loci, but they wouldn't add any information...
    pop = sim.Population(size = N, loci = [Lchr]*nChr, ancGen = -1, infoFields = ['ind_id', 'father_id', 'mother_id', 'gen', 'phenotype', 'breedingValue', 'fitness', 'avgFitness'])
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

    # sim.dump(pop, ancGens = sim.ALL_AVAIL, output = ">>" + savedFolder + "/initialPop.txt")

    ## evolution of the population
    pop.evolve(
        initOps = [
            sim.InitLineage(mode = sim.FROM_INFO_SIGNED), ## save the lineage of each allele
            sim.InitSex(),
            sim.PedigreeTagger(outputFields = ['gen', 'phenotype', 'breedingValue'], output='>>' + savedFolder + '/pedigree.txt'), ## will need the first generation for the randomization
            sim.Stat(alleleFreq = sim.ALL_AVAIL, popSize = True),
            sim.PyOperator(saveLineage, param = [savedFolder]),
            Exporter(format = "csv", output = '>>' + savedFolder + '/genotype.txt', infoFields = ["ind_id"], header = False, genoFormatter = {(0,0):0, (0,1):1, (1,0):1, (1,1):2}, sexFormatter = None, affectionFormatter = None, delimiter = "\t"),
        ],
        preOps = [
            sim.PySelector(func = selModel.fitness),
            ## save the fitness of the parents (offsprings of the generation t - 1)
            ## need to do it now, after having evaluated their fitness via PySelector
            ## if were to do it in postOps with avgFitness: evaluate the average fitness of the PARENTS (not what we want!!!)
            ## for the last generation (offspring), apply the selector outside of pop.evolve to get their fitness
            #sim.PyEval(r"'%d\t' %gen", step = 1, output = '>>' + savedFolder + '/fit/fit' + str(rep) + ".txt"),
            #sim.PyEval(r"'\t'.join(['%.4f' % x for x in fitness]) + '\n'", stmts = "fitness = pop.indInfo('fitness')", exposePop = 'pop', step = 1, output = '>>' + savedFolder + '/fit/fit' + str(rep) + ".txt"),
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
                    sim.PedigreeTagger(outputFields = ['gen', 'phenotype', 'breedingValue'], output='>>' + savedFolder + '/pedigree.txt'),
                ]),
        ),
        postOps = [
            sim.Stat(alleleFreq = sim.ALL_AVAIL, popSize = True),
            ## the additive variance decrease at each generation due to the selection
            ## as a function of the population size
            sim.PyExec("varA = varA*(1 - 0.5/N)"),
            sim.PyOperator(saveLineage, param = [savedFolder]),
            Exporter(format = "csv", output = '>>' + savedFolder + '/genotype.txt', infoFields = ["ind_id"], header = False, genoFormatter = {(0,0):0, (0,1):1, (1,0):1, (1,1):2}, sexFormatter = None, affectionFormatter = None, delimiter = "\t"),
        ],
        gen = G
    )



    ## apply the selector to the offsprings of the last generation to get their fitness
    #lastFit = []
    #for ind in pop.individuals():
    #    lastFit.append(selModel.fitness(ind.phenotype))
    #with open(savedFolder + '/fit/fit' + str(rep) + '.txt', 'a') as f:
    #    f.write(str(G+1) + '\t')
    #    f.write('\t'.join(['%.4f' % x for x in lastFit]) + '\n')

    ## dump the whole population (all generations) in a file
    ## genotypes included, as well as the pedigree informations
    ## for now do that, but will want to save the genotypes only in a specific file
    # sim.dump(pop, ancGens = sim.ALL_AVAIL, output = ">>" + savedFolder + "/genotype.txt")














