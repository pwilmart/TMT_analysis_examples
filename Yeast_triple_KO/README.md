# Gygi Lab Yeast Triple Knockouts

Re-analysis of data from this 2016 publication:

> Paulo, J.A., O’Connell, J.D. and Gygi, S.P., 2016. A triple knockout (TKO) proteomics standard for diagnosing ion interference in isobaric labeling experiments. Journal of the American Society for Mass Spectrometry, 27(10), pp.1620-1625.

Data downloaded from PRIDE PXD008009.

The RAW files were downloaded and reprocessed using two pipelines. One is an in-house pipeline (PAW) developed at Oregon Health & Science University (OHSU) by Phil Wilmarth that uses MSConvert and Comet open source tools. The other processing was done with MaxQuant. The analyses were done to compare the data for the three knock out genes to the data presented in the original publication. Different global factors like numbers MS2 scans for each instrument in respective acquisition modes and numbers of identified PSMs for each LC run were tabulated.

Plots of the total protein reporter ion intensities for each knock out gene were generated to see how platforms and acquisition methods performed. Two platforms were compared that could perform the cleaner SPS MS3 acquisition method: a Thermo Fusion and a Thermo Fusion Lumos mass spectrometer. The Lumos was also used to acquire TMT reporter ions in MS2 scans.

The PAW pipeline converts the RAW files into text files using MSConvert from the Proteowizard toolkit. These text files are parsed to produce MS2 peak lists for database searching and to extract to corresponding reporter ion peak heights.

Data from the MS2 scans are extracted for the Comet searches and the reporter peak heights are extracted from the MS3 scans. The pipeline uses the target/decoy method to make score histograms and determine score filtering thresholds. Accurate mass is used to create conditional score histograms where target/decoy delta mass histograms are used to set the accurate mass windows. Basic parsimony principles are used for protein inference and two peptides per protein were required. An additional protein grouping step was used to combine nearly identical peptide sets (often these are housekeeping genes).

The pipeline also has a script for performing internal reference scaling (IRS) normalization. The three repeats of each 9-plex experiment that were done in separate mass spec runs were combined using the IRS method where averages across the 9 channels served as stand ins for proper reference channels. IRS is described in detail in several Jupyter notebook analyses available online.

MaxQuant version 1.6.2.10 was used with mostly default values and the same protein database as used in the PAW analysis. The Fuson and Lumos SPS MS3 data could be processed together in one analysis (the TMT label set was the same, namely, MS3 10-plex). The MS2 Lumos data was done in a separate analysis (MS2 10-plex). Second peptide identification and precursor isolation purity were not used. Mass tolerances and FDR cutoffs were kept at default values. The proteinGroups file was used to prepare the input for R.

There are some supporting files:

### Summary sheets:
* PAW_grouped_protein_summary_TMT_8.xlsx (main PAW protein summary)
* PAW_grouped_peptide_summary_TMT_8.xlsx (main PAW peptide summary)
* MS2_proteinGroups.xlsx (MaxQuant summary file)
* MS3_proteinGroups.xlsx (MaxQuant summary file)

### IRS script, input files, output files, console log
* pandas_TMT_IRS_norm_first9.py (python script for IRS norm)
* labeled_Fusion_PAW_grouped_protein_summary_TMT_8.txt
* labeled_Lumos_PAW_grouped_protein_summary_TMT_8.txt
* labeled_MS2_Lumos_PAW_grouped_protein_summary_TMT_8.txt
* labeled_Fusion_PAW_grouped_protein_summary_TMT_8_IRS_normalized.txt
* labeled_Lumos_PAW_grouped_protein_summary_TMT_8_IRS_normalized.txt
* labeled_MS2_Lumos_PAW_grouped_protein_summary_TMT_8_IRS_normalized.txt
* IRS_norm_log.txt

### Summary stats and knockout bar plots
* raw_map_stats.xlsx (summary stats)
* KO_bar_plots_Intensity.xlsx (bar plots for intensities)
* KO_bar_plots_S-to-N.xlsx (bar plots for S/N ratios)
