# TMT_analysis_examples
## Examples of TMT data analyses using R
#### Phillip Wilmarth
#### Oregon Health & Science Universtiy
#### Proteomic Shared Resource
#### Winter/Spring 2018

## Other repositories that may be helpful:
* [Multi-TMT experiments and IRS normalization](https://github.com/pwilmart/IRS_normalization.git)

## Folders and descriptions:
### ([Jupyter notebooks](http://jupyter.org) are now viewable via the links below)
> **[MaxQuant_and_PAW](https://pwilmart.github.io/TMT_analysis_examples/KUR1502_MQ_PAW.html)**: Comparison of 3 control versus 4 treatment mouse cell culture data. The data are SPS MS3 on a Thermo Fusion using TMT 10-plex. Data were analyzed with two pipelines: MaxQuant (v1.5.7.4) and an OHSU in-house pipeline (PAW). Anaysis started with protein reports from both (files are in the repository). Details are provided for prepping the data for analysis with R. R analysis script and Jupyter notebook used for analysis.

> **[Dilution_series](https://pwilmart.github.io/TMT_analysis_examples/MAN1353_peptides_proteins.html)**: Analysis of a dilution series to compare the properties of reporter ions at the PSM, the peptide, and the protein levels. This looks at the advantages of aggregating TMT reporter ions into protein intensities.

> **[Multiple_TMT_MQ](https://pwilmart.github.io/TMT_analysis_examples/multiple_TMT_MQ.html)**: Analysis of the mouse lens development data with MaxQuant. This is a three TMT experiment, and how to match the data between TMT experiments is demonstrated. This focuses on normalization methods. Statsitical testing is not explored since that was done in the other repository referenced above.

> **[JPR-201712_MS2-MS3_PSM-peptide-protein](https://pwilmart.github.io/TMT_analysis_examples/MS2MS3_peptides_proteins.html)**: Explores TMT data aggregation at PSM, peptide, and protein levels for the same 10 replicates of an E. coli digest. Also compares MS2-based reporter ions to MS3-based reporter ions. Analysis was done with MaxQuant.

### [Gygi Lab Yeast triple knockout](https://github.com/pwilmart/Yeast_triple_KO_TMT)
**[HTML file](https://pwilmart.github.io/TMT_analysis_examples/Triple_KO.html)**

Re-analysis of yeast triple knockout TMT data from the Gygi lab.

> Paulo, J.A., O’Connell, J.D. and Gygi, S.P., 2016. A triple knockout (TKO) proteomics standard for diagnosing ion interference in isobaric labeling experiments. Journal of the American Society for Mass Spectrometry, 27(10), pp.1620-1625.
