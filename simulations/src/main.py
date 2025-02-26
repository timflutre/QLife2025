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
        with open(savedFolder + '/IBD_allSites.txt', 'a') as f:
            towrite = lineage[:M]
            f.write(str(ind.ind_id) + '\thom1\t' + '\t'.join(['%d' % x for x in towrite]) + '\n')
            towrite = lineage[M:]
            f.write(str(ind.ind_id) + '\thom2\t' + '\t'.join(['%d' % x for x in towrite]) + '\n')
    return True


## function to test the value of the parameters inputed
## if they do not correspond, then exit the simulation
def testParameters(h2, h2_others, nTrait, prop0, nQTLChr, Lchr, varEffect0):
    if h2 <= 0:
        sys.exit('Null heritability is not allowed')
    if len(h2_others) != (nTrait-1):
        sys.exit('Not enough heritabilities defined: each trait need to have a h2')
    if (prop0 > 1.0) or (prop0 < 0.0):
        sys.exit('the proportion of SNPs that are from another distribution need to be between 0 - i.e. all SNPs are either QTLs with effects from the same distribution or neutral - and 1 - i.e. all SNPs have an effect from a gaussian distribution with variance varEffect0')
    if (nQTLChr + round(prop0*(Lchr-nQTLChr))) > Lchr:
        sys.exit('cannot have more QTLs than sites')
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
    import random

    from simuPOP.utils import export
    from simuPOP.utils import Exporter

    default = {'savedFolder' : "default",
               'optim' : 0,
               'varW' : 100,
               'h2': [0.5], ## can have multiple values if more than one trait; it's the first trait that will be selected for
               'G' : 10,
               'N' : 100,
               'Npop' : 10000,
               'nTrait' : 1,
               'nChr' : 1,
               'nQTLChr' : 1, ## number of QTLs per chromosome (whose effect is drawn from N(0,varEffect))
               'Lchr' : 1000,
               'LG' : 100, ## genetic length of one chromosome, in cM
               'mu' : 1e-4,
               'proportion0' : 0.0, ## proportion of QTLs with an effect from a different distribution (effect drawn from N(0,varEffect0))
               'varEffect' : 1.0, ## variance of the distribution of QTL effects (suppose it's the same for all traits if multiple ones)
               'varEffect0' : 0.0, ## consider that we can have two distributioin for the QTLs effects; if this variance is zero then these effects are 0 (neutral sites) -> controlled by proportion0
               'corTrait': 0.5} ## correlation between the QTLs effects of the different traits; consider full pleiotropy (used only for ntrait > 1)

    ## rmk: the variance and covariance for the QTL effect is BEFORE the normalization by the phenotypic variance!!

    parameters = getopts(sys.argv, default)
    
    ## set the seed for simupop and numpy random generator
    SEED = random.randint(0, 2 ** 32)
    SEED = 1234
    random.seed(SEED)
    sim.setRNG(seed = SEED)

    savedFolder = parameters['savedFolder']
    rep = 1

    optim = float(parameters['optim'])
    varW = float(parameters['varW'])

    G = int(parameters['G'])
    N = int(parameters['N'])
    Npop = int(parameters['Npop'])

    ## test if h2 is a list
    ## if not, transform it in a list
    if not isinstance(parameters['h2'], list):
        parameters['h2'] = parameters['h2'].split(sep = ",")

    h2 = float(parameters['h2'][0])
    h2_others = [float(i) for i in parameters['h2'][1:]] ## heritabilities of the other traits; if nTrait = 1, then it will be empty
    nTrait = int(parameters['nTrait'])
    nChr = int(parameters['nChr']) ## the chromosomes are independant
    nQTLChr = int(parameters['nQTLChr']) ## number of QTLs per chromosome // consider that all chromosomes have the same number of QTLs
    Lchr = int(parameters['Lchr']) ## Lchr = sequence length (number of sites), per chromosome // consider that all chromosomes have the same number of sites
    prop0 = float(parameters['proportion0']) ## proportion of markers drawn from N(0, varEffect0) among the markers; if varEffect0 -> these are neutral
    LG = float(parameters['LG'])/100 ## genetic length, in M, of each chromosome
    rho = LG / Lchr ## in M/bp/generation; if want free recombination, there will be one QTL per chromosome and Lchr = 1
    
    mu = float(parameters['mu']) ## mutation rate of the overlay of mutations after the coalescent burn-in phase
    
    
    varEffect0 = float(parameters['varEffect0'])   
    varEffect = float(parameters['varEffect'])
    if nTrait > 1:
        corTrait = float(parameters['corTrait'])
        covTrait = []
        for i in range(nTrait):
            covTrait.append([corTrait * varEffect]*nTrait)
            covTrait[i][i] = varEffect
        parameters['covTrait'] = covTrait
        ## for the other distribution of effects
        covTrait0 = []
        for i in range(nTrait):
            if varEffect0 == 0.0:
                corTrait = 0.0
            covTrait0.append([corTrait * varEffect0]*nTrait)
            covTrait0[i][i] = varEffect0
        parameters['covTrait0'] = covTrait0   


    testParameters(h2, h2_others, nTrait, prop0, nQTLChr, Lchr, varEffect0)

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
    ts = msprime.sim_ancestry(samples = N, sequence_length = Lchr, population_size = Npop, recombination_rate = rho, num_replicates = nChr, random_seed = SEED)

    ## add (neutral) mutations to a tree sequence
    mut_model = msprime.BinaryMutationModel() ## binary model to only have biallelic variants (0 or 1, 0 being the ancestral allele)
    genoList = np.empty(shape = 0)
    variants = [] ## list of positions with variants
    for i, tchr in enumerate(ts):
        mts = msprime.sim_mutations(tchr, rate = mu, model = mut_model)
        var_mts = [var.position for var in mts.variants()]
        variants.append(var_mts)
        tmp = mts.genotype_matrix()
        if len(var_mts) < Lchr: 
            ## if there are less variants than the number of sites given as parameter
            ## add sites that are fixed (allele 0 for all inds)
            tmp = np.array([[0]*(2*N)]*Lchr) ## initialize the genotypes with only 0
            tmp2 = mts.genotype_matrix()
            for jx,jvar in enumerate(var_mts):
                tmp[int(jvar)] = tmp2[jx] ## replace the sites that are polymorphic with their correct genotypes
        genoList = np.append(genoList, tmp) ## export the genotype matrix as a numpy array
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


    ## ID of the SNPs, and their corresponding chromosomes
    snp_id = [str(i+1).rjust(len(str(L)), '0') for i in range(L)]
    chr_id = [[i+1]*Lchr for i in range(nChr)] 
    chr_id = [str(j) for i in chr_id for j in i]
    
    ## reference and alternative alleles for each SNP
    ## randomly drawn
    alleles = np.array(['A', 'T', 'C', 'G'])
    ref = [random.choice(alleles) for i in range(L)]
    alt = [random.choice(alleles[alleles != ref[i]]) for i in range(L)]
    
    ## put NA instead of an alternative allele for the monomorphic sites
    ## i.e. they do not have an alternative allele
    ixVar = [int(k + i*Lchr) for i,j in enumerate(variants) for k in j] ## get the index of the variants in the concatenated list of alleles
    alt = [k if i in ixVar else 'NA' for i,j in enumerate(alt) for k in j]
    
    ## genetic coordinates of each SNP
    ## consider a uniform recombination map, with the same recombination rate between all SNPs
    ## genetic positions = cumulative sum of the recombination rate
    genPos = [ np.array([rho]*Lchr).cumsum() for i in range(nChr)] 
    genPos = [j for i in genPos for j in i]
    
    ## rmk: the snp id, chr id, ref and alt alleles, and the genetic positions will be saved in the same file as the SNPs effects
    
    
    #-------------------- forward-in-time simulations, with selection
    
    ## define which markers are QTLs and which ones are neutral
    ## in the case of proportionQTL < 1
    ## consider that all chromosomes have the same number of QTLs
    ## that are randomly placed along the chromosomes
    if Lchr >= 10:
        L0 = round(prop0*(Lchr-nQTLChr)) ## number of QTL drawn from distribution 0 (varEffect0) per chromosome (proportion of the total number of sites, minus the QTLs from the other distribution)
        qtl = [random.sample([0]*(Lchr-(L0+nQTLChr)) + [1]*nQTLChr + [2]*L0, Lchr) for i in range(nChr)] ## 0 = neutral sites, 1 = QTL, 2 = QTL second distribution of effects
        ## if the number of QTLs is 1
        ## force the QTL to be at a SNP
        if nQTLChr == 1:
            qtl = [random.sample([0]*(Lchr-L0)+ [2]*L0, Lchr) for i in range(nChr)]
            for ichr in range(nChr):
                qtl[ichr][int(random.choice(variants[ichr]))] = 1
        qtl = [j for i in qtl for j in i]
    
    ## if there are less than 10 markers per chromosome 
    ## consider all chromosomes together
    ## (if less than 10, then the proportion of QTLs cannot be correctly represented on a per-chromosome basis)
    if Lchr < 10:
        L0 = round(prop0*(L-nQTLChr*Lchr))
        qtl = random.sample([0]*(L-(L0+nQTLChr)) + [1]*nQTLChr*Lchr + [2]*L0, L)
        ## !! NEED TO ADD CONDITIONS TO CHECK THAT THE QTL IS AT A SNP IF nQTL = 1
    
    ## consider that all the variants can be QTLs
    ## with a given distribution of effect (normal distribution?)
    ## rmk: to define beta, it must not have been defined previously (in the specific scenario)
    if nTrait == 1: 
        beta = list(norm.rvs(size = L, loc = 0, scale = np.sqrt(varEffect)))
        beta0 = list(norm.rvs(size = L, loc = 0, scale = np.sqrt(varEffect0)))
        for ix,i in enumerate(qtl):
            if i == 0:
                beta[ix] = 0.0
            if i == 2:
                beta[ix] = beta0[ix]
    
    if nTrait > 1:
        ## if there are more than one trait, use a multivariate normal to define correlated beta
        ## use the betas of the first trait (arbitrary decision) to calculate the trait that will be used to define the fitness in simuPOP
        ## to not calculate the breeding value and the phenotype of the other traits (can be calculated as posteriori as we have the betas)
        ## in the case there are two different distribution for the QTLs effects:
        ## there will also be two distribution for the other traits (consider full pleiotropy)
        betaAll = multivariate_normal.rvs(size = L, mean = [0]*nTrait, cov = covTrait)
        betaAll0 = multivariate_normal.rvs(size = L, mean = [0]*nTrait, cov = covTrait0)
        for ix,i in enumerate(qtl):
            if i == 0:
                betaAll[ix] = 0.0
            if i == 2: 
                betaAll[ix] = betaAll0[ix]
        ## to define the trait associated in simupop, use the first trait
        beta = [i[0] for i in betaAll]
    
    
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

    ## calculate the genetic values before the standardization (to get varP = 1)
    geneticValue = initGeneticValue(beta, pop)
    varA = np.var(geneticValue)

    ## want to standardize the phenotypic variance to 1
    ## thus, we will standardize the beta by the realized phenotypic standard deviation (ie sqrt(varA/h2))
    if nTrait == 1:
        beta = [i  / np.sqrt(varA/h2) for i in beta]
        with open(savedFolder + "/SNP_INFO_allSites.txt", "w") as out:
            out.write("snp_id\tchr_id\tgen_pos\tREF\tALT\tbeta_trait_1" + '\n')
            for ix,i in enumerate(beta):
                listInfo = [snp_id[ix], chr_id[ix], str(genPos[ix]), ref[ix], alt[ix]]
                out.write('\t'.join(['%s' % x for x in listInfo]) + '\t' + str(i) + "\n")
    
    if nTrait > 1:
        betaAll = betaAll / np.sqrt(varA/h2)
        beta = [i[0] for i in betaAll]
        ## if more than one trait, divide all betas by realized varP 
        with open(savedFolder + "/SNP_INFO_allSites.txt", "w") as out:
            out.write("snp_id\tchr_id\tgen_pos\tREF\tALT\t" + '\t'.join(['beta_trait_' + str(i+1) for i in range(nTrait)]) + '\n')
            for ix,i in enumerate(betaAll):
                listInfo = [snp_id[ix], chr_id[ix], str(genPos[ix]), ref[ix], alt[ix]]
                out.write('\t'.join(['%s' % x for x in listInfo]) + '\t' + '\t'.join(['%.4f' % x for x in list(i)]) + "\n")

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
    
    ## save the header of the genotype file
    with open(savedFolder + "/genotype_allSites.txt", "w") as f:
        f.write("ind_id" + '\t' + '\t'.join(['%s' %x for x in snp_id]) + '\n')

    ## save the header of the IBD file
    with open(savedFolder + "/IBD_allSites.txt", "w") as f:
        f.write("ind_id\thomologous_chr" + '\t' + '\t'.join(['%s' %x for x in snp_id]) + '\n')

    ## evolution of the population
    pop.evolve(
        initOps = [
            sim.InitLineage(mode = sim.FROM_INFO_SIGNED), ## save the lineage of each allele
            sim.InitSex(),
            sim.PedigreeTagger(outputFields = ['gen', 'phenotype', 'breedingValue'], output='>>' + savedFolder + '/pedigree.txt'), ## will need the first generation for the randomization
            sim.Stat(alleleFreq = sim.ALL_AVAIL, popSize = True),
            sim.PyOperator(saveLineage, param = [savedFolder]),
            Exporter(format = "csv", output = '>>' + savedFolder + '/genotype_allSites.txt', infoFields = ["ind_id"], header = False, genoFormatter = {(0,0):0, (0,1):1, (1,0):1, (1,1):2}, sexFormatter = None, affectionFormatter = None, delimiter = "\t"),
        ],
        preOps = [
            sim.PySelector(func = selModel.fitness),
            ## save the fitness of the parents (offsprings of the generation t - 1)
            ## need to do it now, after having evaluated their fitness via PySelector
            ## if were to do it in postOps with avgFitness: evaluate the average fitness of the PARENTS (not what we want!!!)
            ## for the last generation (offspring), apply the selector outside of pop.evolve to get their fitness
            sim.PyEval(r"'%d\t' %gen", step = 1, output = '>>' + savedFolder + '/fitness.txt'),
            sim.PyEval(r"'\t'.join(['%.4f' % x for x in fitness]) + '\n'", stmts = "fitness = pop.indInfo('fitness')", exposePop = 'pop', step = 1, output = '>>' + savedFolder + '/fitness.txt'),
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
            Exporter(format = "csv", output = '>>' + savedFolder + '/genotype_allSites.txt', infoFields = ["ind_id"], header = False, genoFormatter = {(0,0):0, (0,1):1, (1,0):1, (1,1):2}, sexFormatter = None, affectionFormatter = None, delimiter = "\t"),
        ],
        gen = G
    )



    ## apply the selector to the offsprings of the last generation to get their fitness
    lastFit = []
    for ind in pop.individuals():
        lastFit.append(selModel.fitness(ind.phenotype))
    with open(savedFolder + '/fitness.txt', 'a') as f:
        f.write(str(G+1) + '\t')
        f.write('\t'.join(['%.4f' % x for x in lastFit]) + '\n')
    
    ## calculate the breeding value and phenotype of the other trait(s) (i.e. the ones not under selection)
    ## each trait will have a different heritability
    ## need the genotype matrix
    if nTrait > 1:
        geno = np.loadtxt(savedFolder + '/genotype_allSites.txt', skiprows = 1)
        bv = []
        pheno = []
        for i in range(nTrait-1):
            ibeta = [j[i+1] for j in betaAll]
            bv.append([sum([ibeta[iloc] * locus  for iloc,locus in enumerate(ind[1:])]) for ind in geno])
            ## calculate the genetic variance of the initial population 
            ## -> will be used to calculate the environmental variance
            ## and then the phenotype
            u = bv[i][:N]
            varA = np.var(u)
            varE = varA*(1-h2_others[i])/h2_others[i] ## varE is constant with time
            varP = varE + varA
            envNoise = list(norm.rvs(size = len(bv[i]), loc = 0, scale = np.sqrt(varE)))
            phenotype = [sum(x) for x in zip(bv[i], envNoise)]
            pheno.append(phenotype)
    
    ## dump the whole population (all generations) in a file
    ## genotypes included, as well as the pedigree informations
    ## for now do that, but will want to save the genotypes only in a specific file
    # sim.dump(pop, ancGens = sim.ALL_AVAIL, output = ">>" + savedFolder + "/genotype.txt")

    ## add the fitness as a column of the pedigree file
    fitness = np.loadtxt(savedFolder + '/fitness.txt')
    fitness = fitness.tolist()
    fitness = [j for i in fitness for j in i[1:]] ## the first element is the generation

    pedigree = []
    with open(savedFolder + '/pedigree.txt', 'r') as f:
        for ix, i in enumerate(f):
            iline = i.split('\n')[0].split(' ')
            for j in range(nTrait-1):
                iline.append(str(pheno[j][ix]))
                iline.append(str(bv[j][ix]))
            iline.append(str(fitness[ix]))
            pedigree.append(iline)

    ## overwrite the pedigree with the completed one
    with open(savedFolder + '/pedigree.txt', 'w') as f:
        lineHeader = ['ind_id', 'father_id', 'mother_id', 'sex', 'affection', 'generation'] 
        lineTrait = [['pheno' + str(i+1), 'bv' + str(i+1)] for i in range(nTrait)]
        lineTrait = [j for i in lineTrait for j in i]
        lineHeader = lineHeader + lineTrait + ['fitness']
        f.write(' '.join(['%s' %x for x in lineHeader]) + '\n')
        for i in pedigree:
            f.write(' '.join(i) + '\n')

    ## remove the - now useless - fitness file
    os.remove(savedFolder + '/fitness.txt')
    
    ## remove the fixed sites from SNP_INFO.txt, IBD.txt and genotype.txt
    ## still keep the files with all sites ?
    sites = []
    header_sites = []
    with open(savedFolder + '/SNP_INFO_allSites.txt', 'r') as f:
        for ix, i in enumerate(f):
            iline = i.split('\n')[0].split('\t')
            if ix == 0:
                header_sites.append(iline)
            else:
                sites.append(iline)
    
    geno = []
    header_geno = []
    with open(savedFolder + '/genotype_allSites.txt', 'r') as f:
        for ix, i in enumerate(f):
            iline = i.split('\n')[0].split('\t')
            if ix == 0:
                header_geno.append(iline)
            else:
                geno.append(iline)
                
    ibd = []
    header_ibd = []
    with open(savedFolder + '/IBD_allSites.txt', 'r') as f:
        for ix, i in enumerate(f):
            iline = i.split('\n')[0].split('\t')
            if ix == 0:
                header_ibd.append(iline)
            else:
                ibd.append(iline)
    
    sites = np.array(sites)
    geno = np.array(geno)
    ibd = np.array(ibd)
    header_ibd = np.array(header_ibd)
    fixed = np.where(sites[:, 4] == 'NA')[0]
    
    snp = np.delete(sites, fixed, axis = 0)
    genoSNP = np.delete(geno, fixed + 1, axis = 1) ## the first column = ind id
    ibdSNP = np.delete(ibd, fixed + 2, axis = 1) ## the first column = ind id, the second one = homologous chromosome
    header_genoSNP = np.delete(header_geno, fixed + 1, axis = 1) ## the first column = ind id
    header_ibdSNP = np.delete(header_ibd, fixed + 2, axis = 1) ## the first column = ind id
    
    ## save the information about the SNPs only
    with open(savedFolder + '/SNP_INFO.txt', 'w') as f:
        f.write('\t'.join(['%s' %x for x in header_sites[0]]) + '\n')
        for i in snp:
            f.write('\t'.join(i) + '\n')
    
    with open(savedFolder + '/genotype.txt', 'w') as f:
        f.write('\t'.join(['%s' %x for x in header_genoSNP[0]]) + '\n')
        for i in genoSNP:
            f.write('\t'.join(i) + '\n')
    
    with open(savedFolder + '/IBD.txt', 'w') as f:
        f.write('\t'.join(['%s' %x for x in header_ibdSNP[0]]) + '\n')
        for i in ibdSNP:
            f.write('\t'.join(i) + '\n')
        
    print('simulation ' + savedFolder + ' done')





