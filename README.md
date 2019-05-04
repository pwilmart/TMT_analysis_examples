# TMT_analysis_examples
## Examples of TMT data analyses using Jupyter notebooks and R
#### Phillip Wilmarth
#### Oregon Health & Science University, PSR Core
#### 2018, 2019

## Other repositories that may be helpful:
* [Multi-TMT experiments and IRS normalization](https://github.com/pwilmart/IRS_normalization.git)

## Folders and descriptions:
### (rendered HTML of [Jupyter notebooks](http://jupyter.org) are in the links below)

### (1) [MaxQuant_and_PAW repository](https://github.com/pwilmart/MaxQuant_and_PAW.git)
#### [PAW analysis Notebook HTML file](https://pwilmart.github.io/TMT_analysis_examples/KUR1502_PAW.html)
#### [Compare to t-test Notebook HTML file](https://pwilmart.github.io/TMT_analysis_examples/KUR1502_PAW_t-test.html)
#### [Compare to limma Notebook HTML file](https://pwilmart.github.io/TMT_analysis_examples/KUR1502_PAW_limma.html)
#### [Compare to limma-voom Notebook HTML file](https://pwilmart.github.io/TMT_analysis_examples/KUR1502_PAW_limma-voom.html)
#### [MQ analysis Notebook HTML file](https://pwilmart.github.io/TMT_analysis_examples/KUR1502_MQ.html)
##### [older Notebook HTML file](https://pwilmart.github.io/TMT_analysis_examples/KUR1502_PAW.html)

Comparisons of 3 control versus 4 treatment mouse cell culture data. The data are from SPS MS3 acquisition on a Thermo Fusion using TMT 10-plex. Data were analyzed with two pipelines: MaxQuant (v1.6.5.0) and an OHSU in-house pipeline (Comet/PAW). Each pipeline was analyzed separately in similar notebook layouts. Both notebooks can be opened side-by-side for an easy head-to-head comparison.

There is also a comparison of the edgeR statistical testing to a two-sample t-test and to limma (with and without using `voom` variance modeling) to clearly demonstrate the improved performance of newer statistical tools.

Analyses started with protein reports from each pipeline (files are in the repository folders). R analysis scripts and Jupyter notebooks were used for analyses. The data is from the publication below.

> Huan, J., Hornick, N.I., Goloviznina, N.A., Kamimae-Lanning, A.N., David, L.L., Wilmarth, P.A., Mori, T., Chevillet, J.R., Narla, A., Roberts Jr, C.T. and Loriaux, M.M., 2015. Coordinate regulation of residual bone marrow function by paracrine trafficking of AML exosomes. Leukemia, 29(12), p.2285.

### (2) [Dilution_series repository](https://github.com/pwilmart/Dilution_series)
#### [Notebook HTML file](https://pwilmart.github.io/TMT_analysis_examples/MAN1353_peptides_proteins.html)

**Updated January 1, 2019.** Analysis of a dilution series to compare the properties of reporter ions at the PSM, the peptide, and the protein levels. This looks at the advantages of aggregating TMT reporter ions into protein intensities. It also explores data normalization in edgeR with TMM. Updated with better R scripting and more use of ggplot.

### (3) [Multiple_TMT_MQ repository](https://github.com/pwilmart/Multiple_TMT_MQ.git)
#### [Notebook HTML file](https://pwilmart.github.io/TMT_analysis_examples/multiple_TMT_MQ.html)

Analysis of the mouse lens development data with MaxQuant. This is a three TMT experiment, and how to match the data between TMT experiments is demonstrated. This focuses on normalization methods. Statsitical testing is not explored since that was done in the other repository referenced above.

> Khan, S.Y., Ali, M., Kabir, F., Renuse, S., Na, C.H., Talbot, C.C., Hackett, S.F. and Riazuddin, S.A., 2018. Proteome Profiling of Developing Murine Lens Through Mass Spectrometry. Investigative Ophthalmology & Visual Science, 59(1), pp.100-107.

### (4) [JPR-201712_MS2-MS3 repository](https://github.com/pwilmart/JPR-201712_MS2-MS3)
#### [First analysis HTML file](https://pwilmart.github.io/TMT_analysis_examples/MS2MS3_peptides_proteins.html)

Looks at PSM, peptide, and protein level data. Analysis preformed with MaxQuant (v1.5.7.4).

#### [Second analysis HTML file](https://pwilmart.github.io/TMT_analysis_examples/JPR-2017_E-coli_MS2-MS3.html)

Compares MS2 reporter ion data to SPS MS3 reporter ion data. Analysis done with PAW pipeline.

#### [serum analysis HTML file](https://pwilmart.github.io/TMT_analysis_examples/JPR-2017_serum.html)

Analysis of depleted serum labeled with 10-plex TMT and processed on Q-Exactive. Analysis done with PAW pipeline.

Data is from this publication:

> D’Angelo, G., Chaerkady, R., Yu, W., Hizal, D.B., Hess, S., Zhao, W., Lekstrom, K., Guo, X., White, W.I., Roskos, L. and Bowen, M.A., 2017. Statistical models for the analysis of isobaric tags multiplexed quantitative proteomics. Journal of proteome research, 16(9), pp.3124-3136.

### (5) [Gygi Lab Yeast triple knockout repository](https://github.com/pwilmart/Yeast_triple_KO_TMT)
#### [Notebook HTML file](https://pwilmart.github.io/TMT_analysis_examples/Triple_KO.html)

Re-analysis of yeast triple knockout TMT data from the Gygi lab.

> Paulo, J.A., O’Connell, J.D. and Gygi, S.P., 2016. A triple knockout (TKO) proteomics standard for diagnosing ion interference in isobaric labeling experiments. Journal of the American Society for Mass Spectrometry, 27(10), pp.1620-1625.

### (6) [Plubell_2017_PAW repository](https://github.com/pwilmart/Plubell_2017_PAW.git)
#### [Notebook HTML file](https://pwilmart.github.io/TMT_analysis_examples/auto_finder_PAW.html)

Re-analysis of the data in the original IRS paper using PAW/Comet, and other workflows. The original experiment was four 10-plex TMT labelings with the two internal reference pooled standards randomly assigned different channels in each plex. The first notebook explores how to make sure that those two channels in each plex are correctly determined before doing the IRS normalization.

> Plubell, D.L., Wilmarth, P.A., Zhao, Y., Fenton, A.M., Minnier, J., Reddy, A.P., Klimek, J., Yang, X., David, L.L. and Pamir, N., 2017. Extended multiplexing of tandem mass tags (TMT) labeling reveals age and high fat diet specific proteome changes in mouse epididymal adipose tissue. Molecular & Cellular Proteomics, 16(5), pp.873-890.

### (7) [IRS_validation repository](https://github.com/pwilmart/IRS_validation.git)
#### [auto_finder HTML file](https://pwilmart.github.io/TMT_analysis_examples/auto_finder_BIND-473.html)
#### [IRS_validation HTML file](https://pwilmart.github.io/TMT_analysis_examples/IRS_validation.html)

Thorough testing and validation of the IRS method using reference channel data from a 77-channel TMT experiment.

### (8) [Yeast_CarbonSources repository](https://github.com/pwilmart/Yeast_CarbonSources.git)
#### [CarbonSources_part-1 HTML file](https://pwilmart.github.io/TMT_analysis_examples/CarbonSources_part-1.html)
#### [CarbonSources_MQ HTML file](https://pwilmart.github.io/TMT_analysis_examples/CarbonSources_MQ.html)

This is analysis of a public dataset ([PRIDE PXD002875](https://www.ebi.ac.uk/pride/archive/projects/PXD002875)) from Paulo, O'Connell, Gaun, and Gygi processed with the PAW pipeline using Comet. The part-1 notebook covers basic TMT data sanity checks, normalization, and basic statistical testing with edgeR. Part-2 will explore how much the numbers of differential candidates can vary with statistical test choices.

I added a MaxQuant 1.6.3.3 processing of the same RAW files worked up using a very similar notebook.

> Paulo, J.A., O’Connell, J.D., Gaun, A. and Gygi, S.P., 2015. Proteome-wide quantitative multiplexed profiling of protein expression: carbon-source dependency in Saccharomyces cerevisiae. Molecular biology of the cell, 26(22), pp.4063-4074.

### (9) [BCP-ALL_QE-TMT_Nat-Comm-2019 repository](https://github.com/pwilmart/BCP-ALL_QE-TMT_Nat-Comm-2019.git)
#### HTML files:
##### [balanced study averages](https://pwilmart.github.io/TMT_analysis_examples/Nat-Comm-2019_TMT_QE_averages.html) - slightly better IRS using plex averages

##### [single pooled standard](https://pwilmart.github.io/TMT_analysis_examples/Nat-Comm-2019_TMT_QE_pools.html) - single pooled internal standards have a little more uncertainty

Re-analysis of data from childhood acute lymphoblastic leukemia study in Nat. Comm. April 2019. Demonstrates an independent analysis of the 216 Q-Exactive RAW files where MS2 reporter ions are kept in their natural intensity scale instead of ratio transformations. Natural intensity scales are more informative than ratios and have more options for statistical testing. Data from the publication below.

> Yang, M., Vesterlund, M., Siavelis, I., Moura-Castro, L.H., Castor, A., Fioretos, T., Jafari, R., Lilljebjörn, H., Odom, D.T., Olsson, L. and Ravi, N., 2019. Proteogenomics and Hi-C reveal transcriptional dysregulation in high hyperdiploid childhood acute lymphoblastic leukemia. Nature communications, 10(1), p.1519.
