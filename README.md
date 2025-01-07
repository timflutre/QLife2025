---
title: "QLife 2025"
authors: "E. Tourrette, T. Flutre, B. Servin"
---


<!-- pandoc README.md -t html -s -o README.html --toc -->



# Genomic prediction for quanti. traits

```
quarto render genpred.qmd --to html
```

* goal: show how to answer the question "how to replace pedigree-based kinship with SNP-based kinship?"

Simulations:

* 1st: show bias of GEBV when there is selection (train on 10 generations and predict the 11th, no LD, with vs without selection)

* 2nd: influence of recombination (limited and heterogeneous recomb)

On suppose dans cet atelier que les traits ont une archi genet polygeniq infintesimale.

1. comprendre le déterminisme génétique des traits (var genet add et h2)

    * notion de parenté :
    
        * coancestry : ref pop
    
        * matrice de parenté calculée à partir du pedigree (notée "A") : base pop = founders
        
        * genomic relationship matrices ("GRM" notées "K" ici) : ce sont des matrices de var-cov génétiques

                * VanRaden, Yang, Speed (LDAK), Schreck (R/sommer)
                
        * lien entre les différentes notions de parenté : Toro (2011), Legarra (2016), meta-founders
    
    * multi-trait: matrice G, kronecker ; R/sommer, R/MCMCglmm
    
2. prédire les BV d'individus non-phénotypés

    * Henderson MMEs; mentionner le ssGBLUP en perspective
    
    * multi-trait: montrer gain de précision de prédiction
    
    * gBLUP = snpBLUP = rrBLUP



# GWAS for quanti. traits

```
quarto render gwas.qmd --to html
```

* goal: show how to answer the question "how to identify QTLs and model their effects?"

Simulations:

* some SNPs have no effect, no selection

1. SNP-by-SNP GWAS:

    * tests d'hyp: FWER (Bonferroni), FDR

    * R/gaston (implements GEMMA), mais astuce que pour un effet aléatoire
    
    * R/MM4LMM, rapide aussi pour plusieurs effets aléatoires
    
    * ouvrir sur Q et K ; parler de LOCO ; parler de Rio et al (2020)
    
    * imputation et fine mapping : on n'en parle qu'en perspective

    * méta-analyse de différentes cohortes : metal

2. multi-SNP GWAS:

    * single-SNP (en fixe) vs rrBLUP (en aléatoire)
    
    * modéliser la distrib des effets des SNP :
    
        * mélange de Normals : BSLMM, BayesR
        
        * sélection de variables : one-stage R/mlmm.gwas, two-stage R/ash, one-stage spike-and-slab (R/varbvs), lasso/elastic-net (R/glmnet)

3. multi-response GWAS:

    * responses: traits vs envts/cohortes
    
    * joint model: multivariate elastic-net (R/glmnet), R/spring (Chiquet et al)

    * méta-analyse:
    
        * tests d'hyp : R/qch
        
        * modelo effects: R/mash, mantra

* discuss and show how to perform two-stage analysis by deregressing BLUPs

* how to model SNP effects by contrasting backsolving with SNP-by-SNP GWAS

* concretely, we plan to use 2 simulations, one with an infinitesimal architecture and the other with a subset of QTLs with larger effects (>=10 traits) 
