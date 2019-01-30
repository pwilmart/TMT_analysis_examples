# TMT_analysis_examples
## Examples of TMT data analyses using Jupyter notebooks and R
#### Phillip Wilmarth
#### Oregon Health & Science University
#### Proteomic Shared Resource
#### 2018, 2019

## Other repositories that may be helpful:
* [Multi-TMT experiments and IRS normalization](https://github.com/pwilmart/IRS_normalization.git)

## Folders and descriptions:
### ([Jupyter notebooks](http://jupyter.org) are now viewable via the links below)

### (1) [MaxQuant_and_PAW](https://github.com/pwilmart/MaxQuant_and_PAW.git)
#### [Notebook HTML file](https://pwilmart.github.io/TMT_analysis_examples/KUR1502_MQ_PAW.html)

Comparison of 3 control versus 4 treatment mouse cell culture data. The data are SPS MS3 on a Thermo Fusion using TMT 10-plex. Data were analyzed with two pipelines: MaxQuant (v1.5.7.4) and an OHSU in-house pipeline (PAW). Anaysis started with protein reports from both (files are in the repository). Details are provided for prepping the data for analysis with R. R analysis script and Jupyter notebook used for analysis.

> Huan, J., Hornick, N.I., Goloviznina, N.A., Kamimae-Lanning, A.N., David, L.L., Wilmarth, P.A., Mori, T., Chevillet, J.R., Narla, A., Roberts Jr, C.T. and Loriaux, M.M., 2015. Coordinate regulation of residual bone marrow function by paracrine trafficking of AML exosomes. Leukemia, 29(12), p.2285.

### (2) [Dilution_series](https://github.com/pwilmart/Dilution_series)
#### [Notebook HTML file](https://pwilmart.github.io/TMT_analysis_examples/MAN1353_peptides_proteins.html)

**Updated January 1, 2019.** Analysis of a dilution series to compare the properties of reporter ions at the PSM, the peptide, and the protein levels. This looks at the advantages of aggregating TMT reporter ions into protein intensities. It also explores data normalization in edgeR with TMM. Updated with better R scripting and more use of ggplot.

### (3) [Multiple_TMT_MQ](https://github.com/pwilmart/Multiple_TMT_MQ.git)
#### [Notebook HTML file](https://pwilmart.github.io/TMT_analysis_examples/multiple_TMT_MQ.html)

Analysis of the mouse lens development data with MaxQuant. This is a three TMT experiment, and how to match the data between TMT experiments is demonstrated. This focuses on normalization methods. Statsitical testing is not explored since that was done in the other repository referenced above.

> Khan, S.Y., Ali, M., Kabir, F., Renuse, S., Na, C.H., Talbot, C.C., Hackett, S.F. and Riazuddin, S.A., 2018. Proteome Profiling of Developing Murine Lens Through Mass Spectrometry. Investigative Ophthalmology & Visual Science, 59(1), pp.100-107.

 ### (4) [JPR-201712_MS2-MS3_PSM-peptide-protein](https://github.com/pwilmart/JPR-201712_MS2-MS3_PSM-peptide-protein)
 #### [Notebook HTML file](https://pwilmart.github.io/TMT_analysis_examples/MS2MS3_peptides_proteins.html)

 Explores TMT data aggregation at PSM, peptide, and protein levels for the same 10 replicates of an E. coli digest. Also compares MS2-based reporter ions to MS3-based reporter ions. Analysis was done with MaxQuant. Data is from this publication:

 > D’Angelo, G., Chaerkady, R., Yu, W., Hizal, D.B., Hess, S., Zhao, W., Lekstrom, K., Guo, X., White, W.I., Roskos, L. and Bowen, M.A., 2017. Statistical models for the analysis of isobaric tags multiplexed quantitative proteomics. Journal of proteome research, 16(9), pp.3124-3136.

### (5) [Gygi Lab Yeast triple knockout](https://github.com/pwilmart/Yeast_triple_KO_TMT)
#### [Notebook HTML file](https://pwilmart.github.io/TMT_analysis_examples/Triple_KO.html)

Re-analysis of yeast triple knockout TMT data from the Gygi lab.

> Paulo, J.A., O’Connell, J.D. and Gygi, S.P., 2016. A triple knockout (TKO) proteomics standard for diagnosing ion interference in isobaric labeling experiments. Journal of the American Society for Mass Spectrometry, 27(10), pp.1620-1625.

### (6) [Plubell_2017_PAW](https://github.com/pwilmart/Plubell_2017_PAW.git)
#### [Notebook HTML file](https://pwilmart.github.io/TMT_analysis_examples/auto_finder_PAW.html)

Re-analysis of the data in the original IRS paper using PAW/Comet, and other workflows. The original experiment was four 10-plex TMT labelings with the two internal reference pooled standards randomly assigned different channels in each plex. The first notebook explores how to make sure that those two channels in each plex are correctly determined before doing the IRS normalization.

> Plubell, D.L., Wilmarth, P.A., Zhao, Y., Fenton, A.M., Minnier, J., Reddy, A.P., Klimek, J., Yang, X., David, L.L. and Pamir, N., 2017. Extended multiplexing of tandem mass tags (TMT) labeling reveals age and high fat diet specific proteome changes in mouse epididymal adipose tissue. Molecular & Cellular Proteomics, 16(5), pp.873-890.

### (7) [IRS_validation](https://github.com/pwilmart/IRS_validation.git)
#### [auto_finder HTML file](https://pwilmart.github.io/TMT_analysis_examples/auto_finder_BIND-473.html)
#### [IRS_validation HTML file](https://pwilmart.github.io/TMT_analysis_examples/IRS_validation.html)

Thorough testing and validation of the IRS method using reference channel data from a 77-channel TMT experiment.

### (8) [Yeast_CarbonSources](https://github.com/pwilmart/Yeast_CarbonSources.git)
#### [CarbonSources_part-1 HTML file](https://pwilmart.github.io/TMT_analysis_examples/CarbonSources_part-1.html)

This is analysis of a public dataset ([PRIDE PXD002875](https://www.ebi.ac.uk/pride/archive/projects/PXD002875)) from Paulo, O'Connell, Gaun, and Gygi processed with the PAW pipeline using Comet. The part-1 notebook covers basic TMT data sanity checks, normalization, and basic statistical testing with edgeR. Part-2 will explore how much the numbers of differential candidates can vary with statistical test choices.

> Paulo, J.A., O’Connell, J.D., Gaun, A. and Gygi, S.P., 2015. Proteome-wide quantitative multiplexed profiling of protein expression: carbon-source dependency in Saccharomyces cerevisiae. Molecular biology of the cell, 26(22), pp.4063-4074.
