# TMT_analysis_examples
## Examples of TMT data analyses using R 
#### Phillip Wilmarth
#### Oregon Health & Science Universtiy
#### Proteomic Shared Resource
#### Winter/Spring 2018

## Other repositories that may be helpful:
* [Multi-TMT experiments and IRS normalization](https://github.com/pwilmart/IRS_normalization.git)

## Folders and descriptions:
> **MaxQuant_and_PAW**: Comparison of 3 control versus 4 treatment mouse cell culture data. The data are SPS MS3 on a Thermo Fusion using TMT 10-plex. Data were analyzed with two pipelines: MaxQuant (v1.5.7.4) and an OHSU in-house pipeline (PAW). Anaysis started with protein reports from both (files are in the repository). Details are provided for prepping the data for analysis with R. R analysis script and Jupyter notebook used for analysis.

> **Dilution_series**: Analysis of a dilution series to compare the properties of reporter ions at the PSM, the peptide, and the protein levels. This looks at the advantages of aggregating TMT reporter ions into protein intensities.
