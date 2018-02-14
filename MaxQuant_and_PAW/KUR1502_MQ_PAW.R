# load libraries
library(tidyverse) 
library(limma) 
library(edgeR) 

# read the data files (saved as CSV exports from XLSX files)
data_MQ <- read_csv("MQ_prepped_data.csv")
dim(data_MQ)
data_PAW <- read_csv("PAW_prepped_data.csv")
dim(data_PAW)

# save the annotation column and remove from MQ frame
anno_MQ <- data_MQ[8]
data_MQ <- data_MQ[1:7]
row.names(data_MQ) <- anno_MQ$Accession
head(data_MQ)
head(row.names(data_MQ))

# same for PAW
anno_PAW <- data_PAW[8]
data_PAW <- data_PAW[1:7]
row.names(data_PAW) <- anno_PAW$Accession
head(data_PAW)
head(row.names(data_PAW))

# take care of the zeros in MQ frame
data_MQ[data_MQ <= 1] <- 5

# let's see what the starting MQ data look like
# do this in a 2x2 plot frame
par(mfrow = c(2, 2))
boxplot(log2(data_MQ), col = c(rep("red", 3), rep("blue", 4)), 
        notch = TRUE, main = "RAW MQ data")
# NOTE: density distributions order samples differently than box plots...
plotDensities(log2(data_MQ), col = c(rep("blue", 4), rep("red", 3)), main = "Raw MQ data")

# PAW data
boxplot(log2(data_PAW), col = c(rep("red", 3), rep("blue", 4)), 
        notch = TRUE, main = "RAW PAW data")
plotDensities(log2(data_PAW), col = c(rep("blue", 4), rep("red", 3)), main = "Raw PAW data")
par(mfrow = c(1, 1))

# check the column totals to see if they are equal
print("MQ:")
format(round(colSums(data_MQ), digits = 0), big.mark = ",")
print("PAW:")
format(round(colSums(data_PAW), digits = 0), big.mark = ",")

# figure out the global scaling values for sample loading normalizations
target_MQ <- mean(colSums(data_MQ))
target_PAW <- mean(colSums(data_PAW))

# do the sample loading normalization before other normalizations
# there is a different correction factor for each column
norm_facs_MQ <- target_MQ / colSums(data_MQ)
data_MQ_sl <- sweep(data_MQ, 2, norm_facs_MQ, FUN = "*")
norm_facs_PAW <- target_PAW / colSums(data_PAW)
data_PAW_sl <- sweep(data_PAW, 2, norm_facs_PAW, FUN = "*")

par(mfrow = c(2, 2))
# see what the SL normalized data look like
boxplot(log2(data_MQ_sl), col = c(rep("red", 3), rep("blue", 4)), 
        notch = TRUE, main = "SL MQ data")
# NOTE: density distributions order samples differently than box plots...
plotDensities(log2(data_MQ_sl), col = c(rep("blue", 4), rep("red", 3)), main = "SL MQ data")

# PAW data
boxplot(log2(data_PAW_sl), col = c(rep("red", 3), rep("blue", 4)), 
        notch = TRUE, main = "SL PAW data")
plotDensities(log2(data_PAW_sl), col = c(rep("blue", 4), rep("red", 3)), main = "SL PAW data")
par(mfrow = c(1, 1))

# check the columnn totals
format(round(colSums(data_MQ_sl), digits = 0), big.mark = ",")
format(round(colSums(data_PAW_sl), digits = 0), big.mark = ",")

# do TMM on MQ data
MQ_tmm <- calcNormFactors(data_MQ_sl)
data_MQ_tmm <- sweep(data_MQ_sl, 2, MQ_tmm, FUN = "/") # this is data after SL and TMM on original scale

# also on PAW data
PAW_tmm <- calcNormFactors(data_PAW_sl)
data_PAW_tmm <- sweep(data_PAW_sl, 2, PAW_tmm, FUN = "/") # this is data after SL and TMM on original scale

par(mfrow = c(2, 2))
boxplot(log2(data_MQ_tmm), col = c(rep("red", 3), rep("blue", 4)), 
        notch = TRUE, main = "SL/TMM MQ data")
plotDensities(log2(data_MQ_tmm), col = c(rep("blue", 4), rep("red", 3)), main = "SL/TMM MQ data")

# PAW data
bp_mq <- boxplot(log2(data_PAW_tmm), col = c(rep("red", 3), rep("blue", 4)), 
                 notch = TRUE, main = "SL/TMM PAW data")
pd_mq <- plotDensities(log2(data_PAW_tmm), col = c(rep("blue", 4), rep("red", 3)), main = "SL/TMM PAW data")
par(mfrow = c(1, 1))

# check final column totals after TMM
format(round(colSums(data_MQ_tmm), digits = 0), big.mark = ",")
format(round(colSums(data_PAW_tmm), digits = 0), big.mark = ",")

# see how things cluster after we have gotten the boxplots and desity plots looking nice
plotMDS(log2(data_MQ_tmm), col = c(rep("red", 3), rep("blue", 4)), main = "SL/TMM MQ")
plotMDS(log2(data_PAW_tmm), col = c(rep("red", 3), rep("blue", 4)), main = "SL/TMM PAW")

# function computes CVs per time point
make_CVs <- function(df) {
  # separate by samples
  media <- df[1:3]
  exo <- df[4:7]
  
  media$ave <- rowMeans(media)
  media$sd <- apply(media[1:3], 1, sd)
  media$cv <- 100 * media$sd / media$ave
  exo$ave <- rowMeans(exo)
  exo$sd <- apply(exo[1:4], 1, sd)
  exo$cv <- 100 * exo$sd / exo$ave  
  ave_df <- data.frame(media$ave, exo$ave)
  sd_df <- data.frame(media$sd, exo$sd)
  cv_df <- data.frame(media$cv, exo$cv)
  return(list(ave_df, sd_df, cv_df))
}

# get CVs and averages
list_MQ_sl <- make_CVs(data_MQ_sl)
list_PAW_sl <- make_CVs(data_PAW_sl)
list_MQ_tmm <- make_CVs(data_MQ_tmm)
list_PAW_tmm <- make_CVs(data_PAW_tmm)

# compare CV distributions
par(mfrow = c(2, 2))
boxplot(list_MQ_sl[[3]], notch = TRUE, main = "MQ/SL CVs", ylim = c(0, 100))
boxplot(list_MQ_tmm[[3]], notch = TRUE, main = "MQ/SL/TMM CVs", ylim = c(0, 100))
boxplot(list_PAW_sl[[3]], notch = TRUE, main = "PAW/SL CVs", ylim = c(0, 100))
boxplot(list_PAW_tmm[[3]], notch = TRUE, main = "PAW/SL/TMM CVs", ylim = c(0, 100))
par(mfrow = c(1, 1))

# print out the average median CVs
print("MQ (%) (SL then SL/TMM):")
(MQ_sl_med_cv <- round(mean(apply(list_MQ_sl[[3]], 2, median)), 2))
(MQ_tmm_med_cv <- round(mean(apply(list_MQ_tmm[[3]], 2, median)), 2))
print("PAW (%) (SL then SL/TMM):")
(PAW_sl_med_cv <- round(mean(apply(list_PAW_sl[[3]], 2, median)),2))
(PAW_tmm_med_cv <- round(mean(apply(list_PAW_tmm[[3]], 2, median)), 2))

# compare biological replicates to each other by condition
library(psych)
pairs.panels(log2(data_MQ_tmm[1:3]), lm = TRUE, main = "MQ SL/TMM Media")
pairs.panels(log2(data_PAW_tmm[1:3]), lm = TRUE, main = "PAW SL/TMM Media")
pairs.panels(log2(data_MQ_tmm[4:7]), lm = TRUE, main = "MQ SL/TMM Exosome")
pairs.panels(log2(data_PAW_tmm[4:7]), lm = TRUE, main = "PAW SL/TMM Exosome")

# compare within (MQ or PAW) to between (MQ vs PAW) biological replicates
df_MQ <- cbind(anno_MQ, data_MQ_tmm)
df_PAW <- cbind(anno_PAW, data_PAW_tmm)
df_both <- merge(df_MQ, df_PAW, by = "Accession")
dim(df_both)
# head(df_both)

media <- df_both[c(2, 3, 4, 9, 10, 11)]
pairs.panels(log2(media), lm = TRUE, main = "MQ vs PAW, Media")
exo <- df_both[c(5, 6, 7, 8, 12, 13, 14, 15)]
pairs.panels(log2(exo), lm = TRUE, main = "MQ vs PAW, Exosome")

# check what averages of each condition are like between MQ and PAW
library(ggExtra)

# add marginal distrubution histograms to basic correlation plot (good starting point)
ave_MQ <- data.frame(media = rowMeans(data_MQ_tmm[1:3]), exosome = rowMeans(data_MQ_tmm[4:7]))
ggplot()
corr_plot <- ggplot(ave_MQ, aes(x = log10(media), y = log10(exosome))) +
  geom_point() + ggtitle("Media vs Exo: MQ")
ggMarginal(corr_plot, type = "histogram")

ave_PAW <- data.frame(media = rowMeans(data_PAW_tmm[1:3]), exosome = rowMeans(data_PAW_tmm[4:7]))
ggplot()
corr_plot <- ggplot(ave_PAW, aes(x = log10(media), y = log10(exosome))) +
  geom_point() + ggtitle("Media vs Exosome: PAW")
ggMarginal(corr_plot, type = "histogram")

# let's do DE testing with edgeR
# set up the sample mapping
group <- c(rep("media", 3), rep("exo", 4))

# make group into factors and set the order
group <- factor(group, levels = c("media", "exo"))
str(group)

# create a DGEList object with our data
y_MQ <- DGEList(counts = data_MQ_sl, group = group)
y_MQ <- calcNormFactors(y_MQ)
y_MQ <- estimateDisp(y_MQ)

# y_MQ is a list: y_MQ$counts is the data, and y_MQ$samples has interesting content
y_MQ$samples
plotBCV(y_MQ, main = "Biological variation MQ")

# the exact test object has columns like fold-change, CPM, and p-values
et_MQ <- exactTest(y_MQ, pair = c("media", "exo"))
summary(decideTestsDGE(et_MQ)) # this counts up, down, and unchanged genes (here it is proteins)

# the topTags function adds the BH FDR values to an exactTest data frame. Make sure not to change row order!
tt_MQ <- topTags(et_MQ, n = 10000, sort.by = "none")
tt_MQ <- tt_MQ$table # tt_sl is a list. We just need the data frame table

# add the default value as a new column
tt_MQ$candidate <- "no"
tt_MQ[which(tt_MQ$FDR <= 0.10 & tt_MQ$FDR > 0.05), dim(tt_MQ)[2]] <- "low"
tt_MQ[which(tt_MQ$FDR <= 0.05 & tt_MQ$FDR > 0.01), dim(tt_MQ)[2]] <- "med"
tt_MQ[which(tt_MQ$FDR <= 0.01), dim(tt_MQ)[2]] <- "high"
tt_MQ$candidate <- factor(tt_MQ$candidate, levels = c("high", "med",  "low", "no"))

# what does tt_MQ look like?
head(tt_MQ)

# what does the test p-value distribution look like?
ggplot(tt_MQ, aes(PValue)) + 
  geom_histogram(bins = 100, fill = "white", color = "black") + 
  geom_hline(yintercept = mean(hist(tt_MQ$PValue, breaks = 100, plot = FALSE)$counts[26:100])) +
  ggtitle("MQ/edgeR p-value distribution")

# create a DGEList object with our data
y_PAW <- DGEList(counts = data_PAW_sl, group = group)
y_PAW <- calcNormFactors(y_PAW)
y_PAW <- estimateDisp(y_PAW)

# y_PAW is a list: y_PAW$counts is the data, and y_PAW$samples has interesting content
y_PAW$samples
plotBCV(y_PAW, main = "Biological variation PAW")

# the exact test object has columns like fold-change, CPM, and p-values
et_PAW <- exactTest(y_PAW, pair = c("media", "exo"))
summary(decideTestsDGE(et_PAW)) # this counts up, down, and unchanged genes (here it is proteins)

# the topTags function adds the BH FDR values to an exactTest data frame. Make sure not to change row order!
tt_PAW <- topTags(et_PAW, n = 10000, sort.by = "none")
tt_PAW <- tt_PAW$table # tt_PAW is a list. We just need the data frame table

# add the default value as a new column
tt_PAW$candidate <- "no"
tt_PAW[which(tt_PAW$FDR <= 0.10 & tt_PAW$FDR > 0.05), dim(tt_PAW)[2]] <- "low"
tt_PAW[which(tt_PAW$FDR <= 0.05 & tt_PAW$FDR > 0.01), dim(tt_PAW)[2]] <- "med"
tt_PAW[which(tt_PAW$FDR <= 0.01), dim(tt_PAW)[2]] <- "high"
tt_PAW$candidate <- factor(tt_PAW$candidate, levels = c("high", "med",  "low", "no"))

# what does tt_PAW look like?
head(tt_PAW)

# what does the test p-value distribution look like?
ggplot(tt_PAW, aes(PValue)) + 
  geom_histogram(bins = 100, fill = "white", color = "black") + 
  geom_hline(yintercept = mean(hist(tt_PAW$PValue, breaks = 100, plot = FALSE)$counts[26:100])) +
  ggtitle("PAW/edgeR p-value distribution")

# We need one extra library
library(scales)

# function for MA plots
pw_ma_plot <- function(frame, x, y, f, title) {
  # frame = data frame with data
  # x, y are the string names of the x and y columns
  # f is the factor for faceting, title is a string for plot titles
  # make the main MA plot
  temp <- data.frame(log2((frame[x] + frame[y])/2), log2(frame[y] / frame[x]), frame[f])
  colnames(temp) <- c("Ave", "FC", "candidate")
  first  <- ggplot(temp, aes(x = Ave, y = FC)) +
    geom_point(aes(color = candidate, shape = candidate)) +
    scale_y_continuous(paste0("logFC (", x, "/", y, ")")) +
    scale_x_continuous("Ave_intensity") +
    ggtitle(title) + 
    geom_hline(yintercept = 0.0, color = "black") + # one-to-one line
    geom_hline(yintercept = 1.0, color = "black", linetype = "dotted") + # 2-fold up
    geom_hline(yintercept = -1.0, color = "black", linetype = "dotted") # 2-fold down
  
  # make separate MA plots
  second <- ggplot(temp, aes(x = Ave, y = FC)) +
    geom_point(aes(color = candidate, shape = candidate)) +
    scale_y_continuous(paste0("logFC (", x, "/", y, ")")) +
    scale_x_continuous("Ave_intensity") +
    geom_hline(yintercept = 0.0, color = "black") + # one-to-one line
    geom_hline(yintercept = 1.0, color = "black", linetype = "dotted") + # 2-fold up
    geom_hline(yintercept = -1.0, color = "black", linetype = "dotted") + # 2-fold down
    facet_wrap(~ candidate) +
    ggtitle(paste(title, "(separated)", sep=" "))
  # plots inside functions do not automatically display
  print(first)
  print(second)
}

pw_scatter_plot <- function(frame, X, Y, f, title) {
  # frame = data frame with data
  # X, Y are the string names of the x and y columns
  # f is the factor for faceting, title is a string for plot titles
  # make the combined candidate corelation plot
  first <- ggplot(frame, aes_string(X, Y)) +
    geom_point(aes_string(color = f, shape = f)) +
    scale_y_log10() +
    scale_x_log10() +
    ggtitle(title) + 
    geom_abline(intercept = 0.0, slope = 1.0, color = "black") + # one-to-one line
    geom_abline(intercept = 0.301, slope = 1.0, color = "black", linetype = "dotted") + # 2-fold up
    geom_abline(intercept = -0.301, slope = 1.0, color = "black", linetype = "dotted") # 2-fold down
  
  # make separate corelation plots
  second <- ggplot(frame, aes_string(X, Y)) +
    geom_point(aes_string(color = f, shape = f)) +
    scale_y_log10() +
    scale_x_log10() +
    geom_abline(intercept = 0.0, slope = 1.0, color = "black") + # one-to-one line
    geom_abline(intercept = 0.301, slope = 1.0, color = "black", linetype = "dotted") + # 2-fold up
    geom_abline(intercept = -0.301, slope = 1.0, color = "black", linetype = "dotted") + # 2-fold down
    facet_wrap(~ candidate) +
    ggtitle(paste(title, "(separated)", sep=" ")) 
  
  print(first)
  print(second)
}

# MaxQuant first:
# for plotting results, we will use the average intensities for the samples
ave_MQ$candidate <- tt_MQ$candidate
volcano_MQ <- data.frame(log2(ave_MQ$exosome / ave_MQ$media), log10(tt_MQ$FDR)*(-1), ave_MQ$candidate)
colnames(volcano_MQ) <- c("FoldChange", "FDR", "candidate")
head(volcano_MQ)

# start with MA plot
pw_ma_plot(ave_MQ, "media", "exosome", "candidate", "MQ data")
# now the scatter plot
pw_scatter_plot(ave_MQ, "media", "exosome", "candidate", "MQ data")

# make a volcano plot
ggplot(volcano_MQ, aes(x = FoldChange, y = FDR)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  xlab("Fold-Change (Log2)") +
  ylab("-Log10 FDR") +
  ylim(c(NA, 40)) +
  xlim(c(-5, 10)) +
  ggtitle("MQ Volcano Plot")

# PAW pipeline next:
# for plotting results, we will use the average intensities for the samples
ave_PAW$candidate <- tt_PAW$candidate
volcano_PAW <- data.frame(log2(ave_PAW$exosome / ave_PAW$media), log10(tt_PAW$FDR)*(-1), ave_PAW$candidate)
colnames(volcano_PAW) <- c("FoldChange", "FDR", "candidate")
head(volcano_PAW)

# start with MA plot
pw_ma_plot(ave_PAW, "media", "exosome", "candidate", "PAW data")
# now the scatter plot
pw_scatter_plot(ave_PAW, "media", "exosome", "candidate", "PAW data")

# make a volcano plot
ggplot(volcano_PAW, aes(x = FoldChange, y = FDR)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  xlab("Fold-Change (Log2)") +
  ylab("-Log10 FDR") +
  ylim(c(NA, 40)) +
  xlim(c(-5, 10)) +
  ggtitle("PAW Volcano Plot")

# see how many proteins are in each DE category
summary(ave_MQ$candidate)
summary(ave_PAW$candidate)

# collect up some of the data frames and write to disk
final_MQ_frame <- cbind(anno_MQ, data_MQ_tmm, tt_MQ)
write.csv(final_MQ_frame, file = "final_MQ_frame.csv")

# compare the edgeR testing to a more basic t-test
# do the t-test on log transformed intensities to be safe
ttest_MQ <- log2(data_MQ_tmm)
# add average ratio columns (non-logged ratios), fold-change column, and row names
ttest_MQ$ave_media <- rowMeans(data_MQ_tmm[1:3])
ttest_MQ$ave_exo  <- rowMeans(data_MQ_tmm[4:7])
ttest_MQ$logFC <- log2(ttest_MQ$ave_exo / ttest_MQ$ave_media)
row.names(ttest_MQ) <- anno_MQ$Accession

# apply the basic two-sample t-test (we will pool variance)
t.result <- apply(ttest_MQ, 1, function(x) t.test(x[1:3], x[4:7], var.equal = TRUE))
# extract the p-value column from the t-test thingy 
ttest_MQ$p_value <- unlist(lapply(t.result, function(x) x$p.value))
# do a Benjamini-Hochberg multiple testing correction
ttest_MQ$fdr <- p.adjust(ttest_MQ$p_value, method = "BH")

# add a DE candidate status column
ttest_MQ$candidate <- "no"
ttest_MQ[which(ttest_MQ$fdr <= 0.10 & ttest_MQ$fdr > 0.05), dim(ttest_MQ)[2]] <- "low"
ttest_MQ[which(ttest_MQ$fdr <= 0.05 & ttest_MQ$fdr > 0.01), dim(ttest_MQ)[2]] <- "med"
ttest_MQ[which(ttest_MQ$fdr <= 0.01), dim(ttest_MQ)[2]] <- "high"
ttest_MQ$candidate <- factor(ttest_MQ$candidate, levels = c("high", "med",  "low", "no"))
head(ttest_MQ)

# count up, down and the rest (FDR less than 0.05)
all <- dim(ttest_MQ)[1]
up <- dim(ttest_MQ[(ttest_MQ$fdr <= 0.05) & (ttest_MQ$logFC > 0.0), ])[1]
down <- dim(ttest_MQ[(ttest_MQ$fdr <= 0.05) & (ttest_MQ$logFC <= 0.0), ])[1]
print("This is like the decideTest in edgeR - 5% FDR cut:")
up 
all - up - down
down
print("Candidate Counts:")
summary(ttest_MQ$candidate)

# what does the test p-value distribution look like?
ggplot(ttest_MQ, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "white", color = "black") + 
  geom_hline(yintercept = mean(hist(tt_PAW$PValue, breaks = 100, plot = FALSE)$counts[26:100])) +
  ggtitle("MQ data with t-test p-value distribution")

# start with MA plot
pw_ma_plot(ttest_MQ, "ave_media", "ave_exo", "candidate", "MQ t-test data")
# now the scatter plot
pw_scatter_plot(ttest_MQ, "ave_media", "ave_exo", "candidate", "MQ t-test data")

# and the volcano plot
volcano_tt <- data.frame(log2(ttest_MQ$ave_exo / ttest_MQ$ave_media), log10(ttest_MQ$fdr)*(-1), ttest_MQ$candidate)
colnames(volcano_tt) <- c("FoldChange", "FDR", "candidate")
ggplot(volcano_tt, aes(x = FoldChange, y = FDR)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  xlab("Fold-Change (Log2)") +
  ylab("-Log10 FDR") +
  ylim(c(NA, 5)) +
  ggtitle("MQ t-test Volcano Plot")

# replot the edgeR volcano plot on same y-scale for comparison
ggplot(volcano_MQ, aes(x = FoldChange, y = FDR)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  xlab("Fold-Change (Log2)") +
  ylab("-Log10 FDR") +
  ylim(c(NA, 5)) +
  ggtitle("MQ edgeR Volcano Plot")