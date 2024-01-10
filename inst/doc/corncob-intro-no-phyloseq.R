## ----echo = FALSE-------------------------------------------------------------
knitr::opts_chunk$set(fig.width=8, fig.height=4) 

## ----eval = FALSE-------------------------------------------------------------
#  remotes::install_github("statdivlab/corncob")

## ----message = FALSE----------------------------------------------------------
library(corncob)
library(magrittr)
data(soil_phylo_sample)
data(soil_phylo_otu)
data(soil_phylo_taxa)

## -----------------------------------------------------------------------------
head(soil_phylo_sample)

## -----------------------------------------------------------------------------
soil_phylo_otu[1:5, 1:5]

## -----------------------------------------------------------------------------
soil_phylo_taxa[1:3, ]

## -----------------------------------------------------------------------------
data(soil_phylum_small_sample)
sample_data <- soil_phylum_small_sample
data(soil_phylum_small_otu)
data <- soil_phylum_small_otu

## -----------------------------------------------------------------------------
pro_data <- cbind(sample_data, 
                  W = unlist(data["Proteobacteria", ]),
                  M = colSums(data))

## -----------------------------------------------------------------------------
corncob <- bbdml(formula = cbind(W, M - W) ~ 1,
             phi.formula = ~ 1,
             data = pro_data)

## -----------------------------------------------------------------------------
plot(corncob, B = 50)

## -----------------------------------------------------------------------------
plot(corncob, total = TRUE, B = 50)

## -----------------------------------------------------------------------------
plot(corncob, total = TRUE, color = "DayAmdmt", B = 50)

## -----------------------------------------------------------------------------
plot(corncob, color = "DayAmdmt", B = 50)

## -----------------------------------------------------------------------------
corncob_da <- bbdml(formula = cbind(W, M - W) ~ DayAmdmt,
             phi.formula = ~ DayAmdmt,
             data = pro_data)

## -----------------------------------------------------------------------------
plot(corncob_da, color = "DayAmdmt", total = TRUE, B = 50)

## -----------------------------------------------------------------------------
plot(corncob_da, color = "DayAmdmt", B = 50)

## -----------------------------------------------------------------------------
lrtest(mod_null = corncob, mod = corncob_da)

## -----------------------------------------------------------------------------
summary(corncob_da)

## -----------------------------------------------------------------------------
set.seed(1)
da_analysis <- differentialTest(formula = ~ DayAmdmt,
                                 phi.formula = ~ DayAmdmt,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ DayAmdmt,
                                 test = "Wald", boot = FALSE,
                                 data = data,
                                 sample_data = sample_data,
                                 taxa_are_rows = TRUE, 
                                 fdr_cutoff = 0.05)

## -----------------------------------------------------------------------------
da_analysis

## -----------------------------------------------------------------------------
da_analysis$significant_taxa

## -----------------------------------------------------------------------------
set.seed(1)
dv_analysis <- differentialTest(formula = ~ DayAmdmt,
                                 phi.formula = ~ DayAmdmt,
                                 formula_null = ~ DayAmdmt,
                                 phi.formula_null = ~ 1,
                                 test = "LRT", boot = FALSE,
                                 data = data,
                                 sample_data = sample_data,
                                 taxa_are_rows = TRUE, 
                                 fdr_cutoff = 0.05)
dv_analysis$significant_taxa

## -----------------------------------------------------------------------------
da_analysis$p[1:5]

## -----------------------------------------------------------------------------
da_analysis$p_fdr[1:5]

## -----------------------------------------------------------------------------
plot(da_analysis)

## -----------------------------------------------------------------------------
which(is.na(da_analysis$p)) %>% names

## -----------------------------------------------------------------------------
data["GN04", ]

## -----------------------------------------------------------------------------
ex1 <- differentialTest(formula = ~ Day,
                        phi.formula = ~ 1,
                        formula_null = ~ 1,
                        phi.formula_null = ~ 1,
                        data = data,
                        taxa_are_rows = TRUE,
                        sample_data = sample_data, 
                        test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05)
plot(ex1)

## -----------------------------------------------------------------------------
ex2 <- differentialTest(formula = ~ Day,
                        phi.formula = ~ Day,
                        formula_null = ~ 1,
                        phi.formula_null = ~ Day,
                        data = data,
                        taxa_are_rows = TRUE,
                        sample_data = sample_data, 
                        test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05)
plot(ex2)

## -----------------------------------------------------------------------------
ex3 <- differentialTest(formula = ~ Day,
                        phi.formula = ~ Day,
                        formula_null = ~ 1,
                        phi.formula_null = ~ 1,
                        data = data,
                        taxa_are_rows = TRUE,
                        sample_data = sample_data, 
                        test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05)
plot(ex3)

