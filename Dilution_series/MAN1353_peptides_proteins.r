
# comparing peptide intensities to protein intensities
# load libraries
library(tidyverse)
library(limma)

# Set up to select the proper columns and label them
keep <- c(1, 2, 5, 6, 9, 10)
labels = c("a_25", "b_20", "c_15", "d_10", "e_5", "f_2.5")

# read in the three datasets, keep the 6 columns, and set labels
# PSMs
psm_start <- read_csv("psm_tmt.csv")
psm_start <- psm_start[keep]
colnames(psm_start) <- labels
# PAW Grouped peptides
peptide_start <- read_csv("grouped_peptide_summary_TMT_8.csv")
peptide_start <- peptide_start[keep]
colnames(peptide_start) <- labels
# PAW grouped proteins
protein_start <- read_csv("grouped_protein_summary_TMT_8.csv")
protein_start <- protein_start[keep]
colnames(protein_start) <- labels

# filter out rows with any missing data points and see how many rows remain
# this should cull out some low quality PSMs
psms <- psm_start[apply(psm_start, 1, function(x) all(x > 0)), ] 
print("PSMs (before and after):")
dim(psm_start)[1]; dim(psms)[1]
peptides <- peptide_start[apply(peptide_start, 1, function(x) all(x > 0)), ] 
print("Peptides (before and after):")
dim(peptide_start)[1]; dim(peptides)[1]
proteins <- protein_start[apply(protein_start, 1, function(x) all(x > 0)), ] 
print("Proteins (before and after):")
dim(protein_start)[1]; dim(proteins)[1]

# check column sums
format(round(colSums(psms), digits = 0), big.mark = ",")
format(round(colSums(peptides), digits = 0), big.mark = ",")
format(round(colSums(proteins), digits = 0), big.mark = ",")

# compare some intensity (peak heights) distributions
plotDensities(log2(psms), main = "PSMs")

plotDensities(log2(peptides), main = "Peptides")

plotDensities(log2(proteins), main = "Proteins")

# see what we have in the columns at each level of aggregation
print("PSMs:")
summary(psms)
print("Peptides:")
summary(peptides)
print("Proteins")
summary(proteins)

# create an average vector for the x-axis
psms$ref <- rowMeans(psms)
peptides$ref <- rowMeans(peptides)
proteins$ref <- rowMeans(proteins)

# we can simplify plotting if we put data in long form (tidy data)
gpsms <- gather(psms, key = dilution, value = intensity, a_25:f_2.5)
gpeptides <- gather(peptides, key = dilution, value = intensity, a_25:f_2.5)
gproteins <- gather(proteins, key = dilution, value = intensity, a_25:f_2.5)

# check some things
head(gproteins)

# make frames for MA style plots
log_psms <- log2(psms[1:6] / psms$ref)
log_psms$ref <- log10(psms$ref)
log_peptides <- log2(peptides[1:6] / peptides$ref)
log_peptides$ref <- log10(peptides$ref)
log_proteins <- log2(proteins[1:6] / proteins$ref)
log_proteins$ref <- log10(proteins$ref)

# also tidy the log data frames
glog_psms <- gather(log_psms, key = dilution, value = log_ratios, a_25:f_2.5)
glog_peptides <- gather(log_peptides, key = dilution, value = log_ratios, a_25:f_2.5)
glog_proteins <- gather(log_proteins, key = dilution, value = log_ratios, a_25:f_2.5)

# compute the ratios of each dilution channel to the reference
# save as log values for horizontal lines in the MA plots
calc_ratios <- colMeans(psms)
calc_ratios <- log2(calc_ratios[1:6] / calc_ratios[7])

# check some things
round(calc_ratios, 2)

# full scale, linear axes
ggplot(data = gpsms, aes(x = ref, y = intensity)) +
  geom_point(aes(color = dilution, shape = dilution)) +
  geom_smooth(aes(color = dilution), method = "lm") +
  ggtitle("raw PSMs")

# expanded scale
ggplot(data = gpsms, aes(x = ref, y = intensity)) +
  geom_point(aes(color = dilution, shape = dilution)) + 
  geom_smooth(aes(color = dilution), method = "lm") +
  xlim(c(0, 500000)) + ylim(c(0, 500000)) +
  ggtitle("raw PSMs (expanded scale)")

ggplot(data = gpeptides, aes(x = ref, y = intensity)) +
  geom_point(aes(color = dilution, shape = dilution)) + 
  geom_smooth(aes(color = dilution), method = "lm") +
  ggtitle("Peptides")

# expanded scale
ggplot(data = gpeptides, aes(x = ref, y = intensity)) +
  geom_point(aes(color = dilution, shape = dilution)) + 
  geom_smooth(aes(color = dilution), method = "lm") +
  xlim(c(0, 500000)) + ylim(c(0, 500000)) +
  ggtitle("Peptides (expanded scale)")

ggplot(data = gproteins, aes(x = ref, y = intensity)) +
  geom_point(aes(color = dilution, shape = dilution)) + 
  geom_smooth(aes(color = dilution), method = "lm") +
  ggtitle("Proteins")

# expanded scale
ggplot(data = gproteins, aes(x = ref, y = intensity)) +
  geom_point(aes(color = dilution, shape = dilution)) + 
  geom_smooth(aes(color = dilution), method = "lm") +
  xlim(c(0, 1000000)) + ylim(c(0, 1000000)) +
  ggtitle("Proteins (expanded scale)")

# MA style plot for PSMs
ggplot(data = glog_psms, aes(x = ref, y = log_ratios)) +
  geom_point(aes(color = dilution, shape = dilution)) +
  geom_hline(yintercept = calc_ratios) +
  xlab("log ref intensity") + ggtitle("raw PSMs (MA plot)")

# MA style plot for peptides
ggplot(data = glog_peptides, aes(x = ref, y = log_ratios)) +
  geom_point(aes(color = dilution, shape = dilution)) +
  geom_hline(yintercept = calc_ratios) +
  xlab("log ref intensity") + ggtitle("Peptides (MA plot)")

# MA style plot for proteins
ggplot(data = glog_proteins, aes(x = ref, y = log_ratios)) +
  geom_point(aes(color = dilution, shape = dilution)) +
  geom_hline(yintercept = calc_ratios) +
  xlab("log ref intensity") + ggtitle("Proteins (MA plot)")

# load edgeR and then put data into DGEList objects
library(edgeR)
y_psms  <- DGEList(counts = psms[1:6], group = factor(labels))
y_peptides  <- DGEList(counts = peptides[1:6], group = factor(labels))
y_proteins  <- DGEList(counts = proteins[1:6], group = factor(labels))
y_psms$samples
# y_peptides$samples
# y_proteins$samples

# these are the original column sums (they will match the library sizes)
round(colSums(psms), 0)

# dividing by library sizes makes each sample sum to 1.0
colSums(sweep(psms[1:6], 2, y_psms$samples$lib.size, FUN = "/"))

# transforming to the CPM scales makes each column sum to one million
colSums(cpm(y_psms))
sum(colSums(cpm(y_psms)))

# compute TMM factors - they get added to $samples
y_psms_2 <- calcNormFactors(y_psms)
round(y_psms_2$samples$norm.factors, 6)

# see what column sums we get after a cpm function call now that there are TMM factors
round(colSums(cpm(y_psms_2)), 0)
round(sum(colSums(cpm(y_psms_2))), 0)

# combine the library sizes and TMM factors to see what we get
real_factors <- 1 / (y_psms_2$samples$lib.size * y_psms_2$samples$norm.factors)
by_hand_cpm <- 1000000 * sweep(psms[1:6], 2, real_factors, FUN = "*")
round(colSums(by_hand_cpm), 0)

# define a function for scaling data and computing CVs
compute_cv <- function(temp) {
#  middle <- c(25, 20, 15, 10, 5, 2.5) # the "known" fators
  middle <- c(22.26, 17.60, 13.56, 9.43, 4.50, 2.50) # the measured factors
  average <- mean(middle)
  temp[1] <- temp[1] * average/middle[1]
  temp[2] <- temp[2] * average/middle[2]
  temp[3] <- temp[3] * average/middle[3]
  temp[4] <- temp[4] * average/middle[4]
  temp[5] <- temp[5] * average/middle[5]
  temp[6] <- temp[6] * average/middle[6]
  cv <- 100 * apply(temp, 1, sd) / rowMeans(temp)
}

# compute the CVs and look at the distribution summary numbers
summary(compute_cv(psms[1:6]))
summary(compute_cv(peptides[1:6]))
summary(compute_cv(proteins[1:6]))

# we can also let edgeR scale the data for us
# combine the library size and TMM factors (i.e. normalize the data)
y_psms_norm <- calcNormFactors(y_psms)
y_peptides_norm <- calcNormFactors(y_peptides)
y_proteins_norm <- calcNormFactors(y_proteins)

# we can compute the CVs on the cpm function returned objects
cv_psms <- 100 * apply(cpm(y_psms_norm), 1, sd) / rowMeans(cpm(y_psms_norm))
cv_peptides <- 100 * apply(cpm(y_peptides_norm), 1, sd) / rowMeans(cpm(y_peptides_norm))
cv_proteins <- 100 * apply(cpm(y_proteins_norm), 1, sd) / rowMeans(cpm(y_proteins_norm))

# check the distribution summary numbers
summary(cv_psms)
summary(cv_peptides)
summary(cv_proteins)

# see if we have made the dilution series more "the same"
plotMDS(y_proteins_norm)

# compare some box plots of the CV distributions
par(mfrow = c(1, 3))
boxplot(compute_cv(psms[1:6]), notch = TRUE, ylim = c(0, 100), main = "PSMs (13.6%)")
boxplot(compute_cv(peptides[1:6]), notch = TRUE, ylim = c(0, 100), main = "Peptides (11.1%)")
boxplot(compute_cv(proteins[1:6]), notch = TRUE, ylim = c(0, 100), main = "Proteins (6.5%)")
par(mfrow = c(1, 1))

# log the R session
sessionInfo()
