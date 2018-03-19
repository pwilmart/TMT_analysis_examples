# comparing PSM, peptide, and protein intensity properties

# Data from this publication:
#D’Angelo, G., Chaerkady, R., Yu, W., Hizal, D.B., Hess, S., Zhao, W., Lekstrom, K., Guo, X., White, W.I., 
# Roskos, L. and Bowen, M.A., 2017. Statistical models for the analysis of isobaric tags multiplexed 
#quantitative proteomics. Journal of proteome research, 16(9), pp.3124-3136.

# load libraries
library(tidyverse)
library(limma)
library(edgeR)

# read in the four datasets and remove the spike-in proteins at the top of the frames
# PSMs
psm_ms2_all <- read.csv("./MS2/MS2_evidence.csv", stringsAsFactors = FALSE)
print("Raw import numbers:")
dim(psm_ms2_all)
psm_ms2 <- psm_ms2_all[471:dim(psm_ms2_all)[1], 2:11]
psm_ms3_all <- read.csv("./MS3/MS3_evidence.csv", stringsAsFactors = FALSE)
dim(psm_ms3_all)
psm_ms3 <- psm_ms3_all[354:dim(psm_ms3_all)[1], 2:11]

# Peptides
pep_ms2_all <- read.csv("./MS2/MS2_peptides.csv", stringsAsFactors = FALSE)
dim(pep_ms2_all)
pep_ms2 <- pep_ms2_all[219:dim(pep_ms2_all)[1], 2:11]
pep_ms3_all <- read.csv("./MS3/MS3_peptides.csv", stringsAsFactors = FALSE)
dim(pep_ms3_all)
pep_ms3 <- pep_ms3_all[183:dim(pep_ms3_all)[1], 2:11]

# Proteins
prot_ms2_all <- read.csv("./MS2/MS2_proteinGroups.csv", stringsAsFactors = FALSE)
dim(prot_ms2_all)
prot_ms2 <- prot_ms2_all[13:dim(prot_ms2_all)[1], 2:11]
prot_ms3_all <- read.csv("./MS3/MS3_proteinGroups.csv", stringsAsFactors = FALSE)
dim(prot_ms3_all)
prot_ms3 <- prot_ms3_all[13:dim(prot_ms3_all)[1], 2:11]

# replace zeros with "10" (something small relative to the smallest non-zero values)
psm_ms2[psm_ms2 <= 1] <- 10
psm_ms3[psm_ms3 <= 1] <- 10
pep_ms2[pep_ms2 <= 1] <- 10
pep_ms3[pep_ms3 <= 1] <- 10

# use something a bit larger for proteins
prot_ms2[prot_ms2 <= 1] <- 50
prot_ms3[prot_ms3 <= 1] <- 50

# check column sums
sum_frame <- data.frame(format(round(colSums(psm_ms2), digits = 0), big.mark = ","),
                        format(round(colSums(pep_ms2), digits = 0), big.mark = ","),
                        format(round(colSums(prot_ms2), digits = 0), big.mark = ","),
                        format(round(colSums(psm_ms3), digits = 0), big.mark = ","),
                        format(round(colSums(pep_ms3), digits = 0), big.mark = ","),
                        format(round(colSums(prot_ms3), digits = 0), big.mark = ","))
colnames(sum_frame) <- c("PSM_MS2", "Peptide_MS2", "Prot_MS2", "PSM_MS3", "Peptide_MS3", "Prot_MS3")
print("Column sums of different aggregation levels of E. coli proteins:")
sum_frame

# use the short names as working names and save at each stage to more descriptive variable names
psm_ms2_raw <- psm_ms2
psm_ms3_raw <- psm_ms3
pep_ms2_raw <- pep_ms2
pep_ms3_raw <- pep_ms3
prot_ms2_raw <- prot_ms2
prot_ms3_raw <- prot_ms3

# compare some intensity (peak heights) distributions with box plots
par(mfrow = c(1, 2))
boxplot(log10(psm_ms2), notch = TRUE, main = "PSMs MS2")
boxplot(log10(psm_ms3), notch = TRUE, main = "PSMs MS3")
boxplot(log10(pep_ms2), notch = TRUE, main = "Peptides MS2")
boxplot(log10(pep_ms3), notch = TRUE, main = "Peptides MS3")
boxplot(log10(prot_ms2), notch = TRUE, main = "Proteins MS2")
boxplot(log10(prot_ms3), notch = TRUE, main = "Proteins MS3")

# compare some intensity (peak heights) distributions
par(mfrow = c(2, 1))
plotDensities(log10(psm_ms2), main = "PSMs MS2", legend = FALSE)
plotDensities(log10(psm_ms3), main = "PSMs MS3", legend = FALSE)
plotDensities(log10(pep_ms2), main = "Peptides MS2", legend = FALSE)
plotDensities(log10(psm_ms3), main = "Peptides MS3", legend = FALSE)
plotDensities(log10(prot_ms2), main = "Proteins MS2", legend = FALSE)
plotDensities(log10(prot_ms3), main = "Proteins MS3", legend = FALSE)
par(mfrow = c(1, 1))

# see what we have in the columns at each level of aggregation
print("PSM_MS2:")
summary(psm_ms2[6])
print("Peptides_MS2:")
summary(pep_ms2[6])
print("Proteins_MS2")
summary(prot_ms2[6])
print("PSM_MS3:")
summary(psm_ms3[6])
print("Peptides_MS3:")
summary(pep_ms3[6])
print("Proteins_MS3")
summary(prot_ms3[6])

# figure out the global scaling values
tar_psm_ms2 <- mean(colSums(psm_ms2))
tar_pep_ms2 <- mean(colSums(pep_ms2))
tar_prot_ms2 <- mean(colSums(prot_ms2))
tar_psm_ms3 <- mean(colSums(psm_ms3))
tar_pep_ms3 <- mean(colSums(pep_ms3))
tar_prot_ms3 <- mean(colSums(prot_ms3))

# do the sample loading normalization before the IRS normalization
# there is a different correction factor for each column
facs_psm_ms2 <- tar_psm_ms2 / colSums(psm_ms2)
psm_ms2 <- sweep(psm_ms2, 2, facs_psm_ms2, FUN = "*")
facs_pep_ms2 <- tar_pep_ms2 / colSums(pep_ms2)
pep_ms2 <- sweep(pep_ms2, 2, facs_pep_ms2, FUN = "*")
facs_prot_ms2 <- tar_prot_ms2 / colSums(prot_ms2)
prot_ms2 <- sweep(prot_ms2, 2, facs_prot_ms2, FUN = "*")
facs_psm_ms3 <- tar_psm_ms3 / colSums(psm_ms3)
psm_ms3 <- sweep(psm_ms3, 2, facs_psm_ms3, FUN = "*")
facs_pep_ms3 <- tar_pep_ms3 / colSums(pep_ms3)
pep_ms3 <- sweep(pep_ms3, 2, facs_pep_ms3, FUN = "*")
facs_prot_ms3 <- tar_prot_ms3 / colSums(prot_ms3)
prot_ms3 <- sweep(prot_ms3, 2, facs_prot_ms3, FUN = "*")

# see what the SL normalized data look like
par(mfrow = c(1, 2))
boxplot(log10(psm_ms2), notch = TRUE, main = "SL PSM MS2 data")
boxplot(log10(psm_ms3), notch = TRUE, main = "SL PSM MS3 data")
boxplot(log10(pep_ms2), notch = TRUE, main = "SL Pep MS2 data")
boxplot(log10(pep_ms3), notch = TRUE, main = "SL Pep MS3 data")
boxplot(log10(prot_ms2), notch = TRUE, main = "SL Prot MS2 data")
boxplot(log10(prot_ms3), notch = TRUE, main = "SL Prot MS3 data")

# NOTE: density distributions order samples differently than box plots...
par(mfrow = c(2, 1))
plotDensities(log10(psm_ms2), legend = FALSE, main = "SL PSM MS2 data")
plotDensities(log10(psm_ms3), legend = FALSE, main = "SL PSM MS3 data")
plotDensities(log10(pep_ms2), legend = FALSE, main = "SL Pep MS2 data")
plotDensities(log10(pep_ms3), legend = FALSE, main = "SL Pep MS3 data")
plotDensities(log10(prot_ms2), legend = FALSE, main = "SL Prot MS2 data")
plotDensities(log10(prot_ms3), legend = FALSE, main = "SL Prot MS3 data")
par(mfrow = c(1, 1))

# do TMM on data
psm_ms2_tmm <- calcNormFactors(psm_ms2)
psm_ms2 <- sweep(psm_ms2, 2, psm_ms2_tmm, FUN = "/") # this is data after SL and TMM on original scale
pep_ms2_tmm <- calcNormFactors(pep_ms2)
pep_ms2 <- sweep(pep_ms2, 2, pep_ms2_tmm, FUN = "/") # this is data after SL and TMM on original scale
print("Print the TMM correction factors (MS2 proteins then MS3 proteins:")
(prot_ms2_tmm <- calcNormFactors(prot_ms2))
prot_ms2 <- sweep(prot_ms2, 2, prot_ms2_tmm, FUN = "/") # this is data after SL and TMM on original scale

psm_ms3_tmm <- calcNormFactors(psm_ms3)
psm_ms3 <- sweep(psm_ms2, 2, psm_ms3_tmm, FUN = "/") # this is data after SL and TMM on original scale
pep_ms3_tmm <- calcNormFactors(pep_ms3)
pep_ms3 <- sweep(pep_ms2, 2, pep_ms3_tmm, FUN = "/") # this is data after SL and TMM on original scale
(prot_ms3_tmm <- calcNormFactors(prot_ms3))
prot_ms3 <- sweep(prot_ms2, 2, prot_ms3_tmm, FUN = "/") # this is data after SL and TMM on original scale

# see what the SL/TMM normalized data look like
par(mfrow = c(1, 2))
boxplot(log10(psm_ms2), notch = TRUE, main = "SL/TMM PSM MS2 data")
boxplot(log10(psm_ms3), notch = TRUE, main = "SL/TMM PSM MS3 data")
boxplot(log10(pep_ms2), notch = TRUE, main = "SL/TMM Pep MS2 data")
boxplot(log10(pep_ms3), notch = TRUE, main = "SL/TMM Pep MS3 data")
boxplot(log10(prot_ms2), notch = TRUE, main = "SL/TMM Prot MS2 data")
boxplot(log10(prot_ms3), notch = TRUE, main = "SL/TMM Prot MS3 data")

# NOTE: density distributions order samples differently than box plots...
par(mfrow = c(2, 1))
plotDensities(log10(psm_ms2), legend = FALSE, main = "SL/TMM PSM MS2 data")
plotDensities(log10(psm_ms3), legend = FALSE, main = "SL/TMM PSM MS3 data")
plotDensities(log10(pep_ms2), legend = FALSE, main = "SL/TMM Pep MS2 data")
plotDensities(log10(pep_ms3), legend = FALSE, main = "SL/TMM Pep MS3 data")
plotDensities(log10(prot_ms2), legend = FALSE, main = "SL/TMM Prot MS2 data")
plotDensities(log10(prot_ms3), legend = FALSE, main = "SL/TMM Prot MS3 data")
par(mfrow = c(1, 1))

# save the normalized data frames
psm_ms2_norm <- psm_ms2
psm_ms3_norm <- psm_ms3
pep_ms2_norm <- pep_ms2
pep_ms3_norm <- pep_ms3
prot_ms2_norm <- prot_ms2
prot_ms3_norm <- prot_ms3

# create an average vector for the x-axis
psm_ms2$ref <- rowMeans(psm_ms2)
pep_ms2$ref <- rowMeans(pep_ms2)
prot_ms2$ref <- rowMeans(prot_ms2)
psm_ms3$ref <- rowMeans(psm_ms3)
pep_ms3$ref <- rowMeans(pep_ms3)
prot_ms3$ref <- rowMeans(prot_ms3)

# linear scales (thank you stack exchange - plotting setup in ggplot is very different...)
psm_ms2 %>% 
  gather(channel, intensity, -ref) %>%
  ggplot(aes(ref, intensity)) + geom_point() + facet_wrap(~channel) + ggtitle("PSMs MS2")

pep_ms2 %>% 
  gather(channel, intensity, -ref) %>%
  ggplot(aes(ref, intensity)) + geom_point() + facet_wrap(~channel) + ggtitle("Peptides MS2")

prot_ms2 %>% 
  gather(channel, intensity, -ref) %>%
  ggplot(aes(ref, intensity)) + geom_point() + facet_wrap(~channel) + ggtitle("Proteins MS2")

# plot the linear MS3 method data
psm_ms3 %>% 
  gather(channel, intensity, -ref) %>%
  ggplot(aes(ref, intensity)) + geom_point() + facet_wrap(~channel) + ggtitle("PSMs MS3")

pep_ms3 %>% 
  gather(channel, intensity, -ref) %>%
  ggplot(aes(ref, intensity)) + geom_point() + facet_wrap(~channel) + ggtitle("Peptides MS3")

prot_ms3 %>% 
  gather(channel, intensity, -ref) %>%
  ggplot(aes(ref, intensity)) + geom_point() + facet_wrap(~channel) + ggtitle("Proteins MS3")

# log transform the data
log_psm_ms2 <- log10(psm_ms2)
log_pep_ms2 <- log10(pep_ms2)
log_prot_ms2 <- log10(prot_ms2)
log_psm_ms3 <- log10(psm_ms3)
log_pep_ms3 <- log10(pep_ms3)
log_prot_ms3 <- log10(prot_ms3)

# plot the MS2 method data with the log scales
log_psm_ms2 %>% 
  gather(channel, intensity, -ref) %>%
  ggplot(aes(ref, intensity)) + geom_point() + facet_wrap(~channel) + ggtitle("log PSMs MS2")

log_pep_ms2 %>% 
  gather(channel, intensity, -ref) %>%
  ggplot(aes(ref, intensity)) + geom_point() + facet_wrap(~channel) + ggtitle("log Peptides MS2")

log_prot_ms2 %>% 
  gather(channel, intensity, -ref) %>%
  ggplot(aes(ref, intensity)) + geom_point() + facet_wrap(~channel) + ggtitle("log Proteins MS2")

# finally, the MS3 method data with log scales
log_psm_ms3 %>% 
  gather(channel, intensity, -ref) %>%
  ggplot(aes(ref, intensity)) + geom_point() + facet_wrap(~channel) + ggtitle("log PSMs MS3")

log_pep_ms3 %>% 
  gather(channel, intensity, -ref) %>%
  ggplot(aes(ref, intensity)) + geom_point() + facet_wrap(~channel) + ggtitle("log Peptides MS3")

log_prot_ms3 %>% 
  gather(channel, intensity, -ref) %>%
  ggplot(aes(ref, intensity)) + geom_point() + facet_wrap(~channel) + ggtitle("log Proteins MS3")

# define a function for computing CVs
compute_cv <- function(temp) {
  cv <- 100 * apply(temp, 1, sd) / rowMeans(temp)
}

# compute the CVs and look at the distributions
print("MS2 method (PSMs, peptides, proteins)")
summary(compute_cv(psm_ms2))
summary(compute_cv(pep_ms2))
summary(compute_cv(prot_ms2))

print("MS3 method (PSMs, peptides, proteins)")
summary(compute_cv(psm_ms3))
summary(compute_cv(pep_ms3))
summary(compute_cv(prot_ms3))

# make some box plots, too
par(mfrow = c(2, 3))
boxplot(compute_cv(psm_ms2), ylim = c(0, 50), notch = TRUE, main = "MS2: PSMs")
boxplot(compute_cv(pep_ms2), ylim = c(0, 50), notch = TRUE, main = "MS2: Peptides")
boxplot(compute_cv(prot_ms2), ylim = c(0, 50), notch = TRUE, main = "MS2: Proteins")
boxplot(compute_cv(psm_ms3), ylim = c(0, 50), notch = TRUE, main = "MS3: PSMs")
boxplot(compute_cv(pep_ms3), ylim = c(0, 50), notch = TRUE, main = "MS3: Peptides")
boxplot(compute_cv(prot_ms3), ylim = c(0, 50), notch = TRUE, main = "MS3: Proteins")
par(mfrow = c(1, 1))

# check clustering at the protein level betwen MS2 and MS3 methods
plotMDS(log2(prot_ms2), pch = 4, main = "MS2 Proteins")
plotMDS(log2(prot_ms3), pch = 3, main = "MS3 Proteins")

# remember that we saved the final normalized data in data frames that we could 
# write out if we wanted to get back into static summary files for publication