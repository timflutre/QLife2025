---
title: "QLife 2025 (PSL/ENS)"
authors: "E. Tourrette, T. Flutre, B. Servin"
---


<!-- pandoc README.md -t html -s -o README.html --toc -->



# Genomic prediction for quanti. traits

See the file `genpred.qmd`:

1. open it in RStudio 

2. run the simulations in the terminal

3. click on `Render`

It is also possible to execute it in a terminal with this command-line:
```
quarto render genpred.qmd --to html
```

Simulations:

* 1st: show bias of GEBV when there is selection (train on 10 generations and predict the 11th, no LD, with vs without selection)

* 2nd: influence of recombination (limited and heterogeneous recomb)



# GWAS for quanti. traits

Same but with the file `gwas.qmd`.
