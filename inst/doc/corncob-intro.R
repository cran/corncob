## ----echo = FALSE-------------------------------------------------------------
knitr::opts_chunk$set(fig.width=8, fig.height=4) 

## ----results = 'hide'---------------------------------------------------------
phy <- requireNamespace("phyloseq", quietly = TRUE) == TRUE

## ----echo = FALSE-------------------------------------------------------------
print(paste0("phyloseq is installed: ", phy))

## ----eval = FALSE-------------------------------------------------------------
# remotes::install_github("statdivlab/corncob")

## ----message = FALSE, eval = phy----------------------------------------------
library(corncob)
library(phyloseq)
library(magrittr)
data(soil_phylo_sample)
data(soil_phylo_otu)
data(soil_phylo_taxa)
soil_phylo <- phyloseq::phyloseq(phyloseq::sample_data(soil_phylo_sample),
                                phyloseq::otu_table(soil_phylo_otu, taxa_are_rows = TRUE),
                                phyloseq::tax_table(soil_phylo_taxa))

## ----eval = phy---------------------------------------------------------------
soil_phylo

## ----eval = phy---------------------------------------------------------------
otu_table(soil_phylo)[1:3, 1:3]

## ----eval = phy---------------------------------------------------------------
sample_data(soil_phylo)[1:3, ]

## ----eval = phy---------------------------------------------------------------
tax_table(soil_phylo)[1:3, ]

## ----eval = phy---------------------------------------------------------------
soil <- soil_phylo %>% 
            phyloseq::subset_samples(DayAmdmt %in% c(11,21)) %>%
            phyloseq::tax_glom("Phylum") 

## ----eval = phy---------------------------------------------------------------
soil

## ----eval = phy---------------------------------------------------------------
tax_table(soil)[1:5, ]

## ----eval = phy---------------------------------------------------------------
corncob <- bbdml(formula = OTU.1 ~ 1,
             phi.formula = ~ 1,
             data = soil)

## ----eval = phy---------------------------------------------------------------
plot(corncob, B = 50)

## ----eval = phy---------------------------------------------------------------
plot(corncob, total = TRUE, B = 50)

## ----eval = phy---------------------------------------------------------------
plot(corncob, total = TRUE, color = "DayAmdmt", B = 50)

## ----eval = phy---------------------------------------------------------------
plot(corncob, color = "DayAmdmt", B = 50)

## -----------------------------------------------------------------------------
df <- plot(corncob, color = "DayAmdmt", B = 50, data_only = TRUE)
head(df)

## ----eval = phy---------------------------------------------------------------
corncob_da <- bbdml(formula = OTU.1 ~ DayAmdmt,
             phi.formula = ~ DayAmdmt,
             data = soil)

## ----eval = phy---------------------------------------------------------------
plot(corncob_da, color = "DayAmdmt", total = TRUE, B = 50)

## ----eval = phy---------------------------------------------------------------
plot(corncob_da, color = "DayAmdmt", B = 50)

## ----eval = phy---------------------------------------------------------------
lrtest(mod_null = corncob, mod = corncob_da)

## ----eval = phy---------------------------------------------------------------
summary(corncob_da)

## ----eval = phy---------------------------------------------------------------
set.seed(1)
da_analysis <- differentialTest(formula = ~ DayAmdmt,
                                 phi.formula = ~ DayAmdmt,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ DayAmdmt,
                                 test = "Wald", boot = FALSE,
                                 data = soil,
                                 fdr_cutoff = 0.05)

## ----eval = phy---------------------------------------------------------------
da_analysis

## ----eval = phy---------------------------------------------------------------
da_analysis$significant_taxa

## ----eval = phy---------------------------------------------------------------
set.seed(1)
dv_analysis <- differentialTest(formula = ~ DayAmdmt,
                                 phi.formula = ~ DayAmdmt,
                                 formula_null = ~ DayAmdmt,
                                 phi.formula_null = ~ 1,
                                 data = soil,
                                 test = "LRT", boot = FALSE,
                                 fdr_cutoff = 0.05)
dv_analysis$significant_taxa

## ----eval = phy---------------------------------------------------------------
otu_to_taxonomy(OTU = da_analysis$significant_taxa, data = soil)

## ----eval = phy---------------------------------------------------------------
otu_to_taxonomy(OTU = dv_analysis$significant_taxa, data = soil)

## ----eval = phy---------------------------------------------------------------
da_analysis$p[1:5]

## ----eval = phy---------------------------------------------------------------
da_analysis$p_fdr[1:5]

## ----eval = phy---------------------------------------------------------------
plot(da_analysis)

## -----------------------------------------------------------------------------
df <- plot(da_analysis, data_only = TRUE)

# we can easily remove special characters used in our formatting steps
df <- df %>%
  dplyr::mutate(variable = gsub("\nDifferential Abundance", "",
                                variable, fixed = TRUE))

head(df)

## ----eval = phy---------------------------------------------------------------
which(is.na(da_analysis$p)) %>% names

## ----eval = phy---------------------------------------------------------------
otu_to_taxonomy(OTU = "OTU.4206", data = soil)

## ----eval = phy---------------------------------------------------------------
otu_table(soil)["OTU.4206"]

## ----eval = phy---------------------------------------------------------------
check_GN04 <- bbdml(formula = OTU.4206 ~ DayAmdmt,
                 phi.formula = ~ DayAmdmt,
                 data = soil)

## ----eval = phy---------------------------------------------------------------
soil_full <- soil_phylo %>% 
  tax_glom("Phylum") 

## ----eval = phy---------------------------------------------------------------
ex1 <- differentialTest(formula = ~ Day,
                                 phi.formula = ~ 1,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ 1,
                                 data = soil_full,
                                 test = "Wald", boot = FALSE,
                                 fdr_cutoff = 0.05)
plot(ex1)

## ----eval = phy---------------------------------------------------------------
ex2 <- differentialTest(formula = ~ Day,
                                 phi.formula = ~ Day,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ Day,
                                 data = soil_full,
                                 test = "Wald", boot = FALSE,
                                 fdr_cutoff = 0.05)
plot(ex2)

## ----eval = phy---------------------------------------------------------------
ex3 <- differentialTest(formula = ~ Day,
                                 phi.formula = ~ Day,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ 1,
                                 data = soil_full,
                                 test = "Wald", boot = FALSE,
                                 fdr_cutoff = 0.05)
plot(ex3)

## ----eval = phy---------------------------------------------------------------
ex4 <- differentialTest(formula = ~ Day + Amdmt,
                                 phi.formula = ~ Day,
                                 formula_null = ~ Amdmt,
                                 phi.formula_null = ~ 1,
                                 data = soil_full,
                                 test = "Wald", boot = FALSE,
                                 fdr_cutoff = 0.05)
plot(ex4)

## ----eval = phy---------------------------------------------------------------
ex5 <- differentialTest(formula = ~ Day + Amdmt,
                                 phi.formula = ~ Day + Amdmt,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ Day + Amdmt,
                                 data = soil_full,
                                 test = "Wald", boot = FALSE,
                                 fdr_cutoff = 0.05)
plot(ex5)

## ----eval = phy---------------------------------------------------------------
ex6 <- differentialTest(formula = ~ Day,
                                 phi.formula = ~ Day + Amdmt,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ Day,
                                 data = soil_full,
                                 test = "Wald", boot = FALSE,
                                 fdr_cutoff = 0.05)
plot(ex6)

## ----eval = phy---------------------------------------------------------------
data(ibd_phylo_sample)
data(ibd_phylo_otu)
data(ibd_phylo_taxa)
ibd_phylo <-  phyloseq::phyloseq(phyloseq::sample_data(ibd_phylo_sample), 
                                 phyloseq::otu_table(ibd_phylo_otu, taxa_are_rows = TRUE),
                                 phyloseq::tax_table(ibd_phylo_taxa))
ibd <- ibd_phylo %>% 
            phyloseq::tax_glom("Genus") 

## ----eval = phy---------------------------------------------------------------
ex7 <- differentialTest(formula = ~ ibd,
                                 phi.formula = ~ 1,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ 1,
                                 data = ibd,
                                 test = "Wald", boot = FALSE,
                                 fdr_cutoff = 0.05)
plot(ex7)

## ----eval = phy---------------------------------------------------------------
plot(ex7, level = c("Family", "Genus"))

## ----eval = phy---------------------------------------------------------------
ex8 <- differentialTest(formula = ~ ibd,
                                 phi.formula = ~ ibd,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ 1,
                                 data = ibd,
                                 test = "Wald", boot = FALSE,
                                 fdr_cutoff = 0.05)
plot(ex8, level = "Genus")

