# TMT_analysis_examples
## Examples of TMT data analyses using R
#### Phillip Wilmarth
#### Oregon Health & Science Universtiy
#### Proteomic Shared Resource
#### Winter/Spring 2018

## Other repositories that may be helpful:
* [Multi-TMT experiments and IRS normalization](https://github.com/pwilmart/IRS_normalization.git)

## Repositories and descriptions:
> **[MaxQuant_and_PAW](https://github.com/pwilmart/MaxQuant_and_PAW.git)**: Comparison of 3 control versus 4 treatment mouse cell culture data. The data are SPS MS3 on a Thermo Fusion using TMT 10-plex. Data were analyzed with two pipelines: MaxQuant (v1.5.7.4) and an OHSU in-house pipeline (PAW). Anaysis started with protein reports from both (files are in the repository). Details are provided for prepping the data for analysis with R. R analysis script and Jupyter notebook used for analysis. View the notebook at [this link](KUR1502_MQ_PAW.html)

> **[Dilution_series](https://github.com/pwilmart/Dilution_series.git)**: Analysis of a dilution series to compare the properties of reporter ions at the PSM, the peptide, and the protein levels. This looks at the advantages of aggregating TMT reporter ions into protein intensities. View the notebook at [this link](MAN1353_peptides_proteins.html).

> **[Multiple_TMT_MQ](https://github.com/pwilmart/Multiple_TMT_MQ.git)**: Analysis of the mouse lens development data with MaxQuant. This is a three TMT experiment, and how to match the data between TMT experiments is demonstrated. This focuses on normalization methods. Statsitical testing is not explored since that was done in the other repository referenced above. View the notebook at [this link](multiple_TMT_MQ.html).

> **[JPR-201712_MS2-MS3_PSM-peptide-protein](https://github.com/pwilmart/JPR-201712_MS2-MS3_PSM-peptide-protein.git)**: Explores TMT data aggregation at PSM, peptide, and protein levels for the same 10 replicates of an E. coli digest. Also compares MS2-based reporter ions to MS3-based reporter ions. Analysis was done with MaxQuant. View the notebook at [this link](MS2MS3_peptides_proteins.html).

> **[Yeast_triple_KO_TMT](https://github.com/pwilmart/Yeast_triple_KO_TMT.git)**: Re-analysis of the yeast triple knock-out data from the Gygi lab. This analysis using the PAW/Comet pipeline uses IRS normalization to better integrate data from the replicate TMT experiments. This data compares newer SPS MS3 methods to older MS2 methods. View the notebook at [this link](Triple_KO.html).

> **[Plubell_2017_PAW](https://github.com/pwilmart/Plubell_2017_PAW.git)**: Re-analysis of the data in the original IRS paper using PAW/Comet, and other workflows. The original experiment was four 10-plex TMT labelings with the two internal reference pooled standards randomly assigned different channels in each plex. The first notebook explores how to make sure that those two channels in each plex are correctly determined before doing the IRS normalization.  View the "auto_finder" notebook at [this link](auto_finder_PAW.html).
