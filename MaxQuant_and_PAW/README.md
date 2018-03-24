## Comparison between MaxQuant and PAW pipelines. Also comparison between edgeR statistical testing and basic two-sample t-test.
**Click on the Jupyter notebook file (_KUR1502_MQ_PAW.ipynb_) to see the notebook in your browser. It may take a minute for the page to render and load, so please be patient.**

## Data from this publication:
> Huan, J., Hornick, N.I., Goloviznina, N.A., Kamimae-Lanning, A.N., David, L.L., Wilmarth, P.A., Mori, T., Chevillet, J.R., Narla, A., Roberts Jr, C.T. and Loriaux, M.M., 2015. 
Coordinate regulation of residual bone marrow function by paracrine trafficking of AML exosomes. Leukemia, 29(12), p.2285.

This is a mouse bone marrow cell culture experiment with controls (n=3) and leukemia exosome-dosed cells (n=4). 
It is a single TMT labeling design. MaxQuant (v1.5.7.4) was used with a mouse Swiss-Prot protein database. 10-plex TMT and Reporter 
ion MS3 quant were selected. Reporter ion tolerance was 0.003 Da. Data was taken from the proteinGroups.txt summary file. Full details
are provided on how data was prepped for importing into R. An R script, suitable for use in RStudio, and a Jupyter notebook are 
provided for reproducing the analysis. An HTML report has been generated. Excel spreadsheets are present for both pipelines to 
illustrate what starting data look like. The R analysis results are exported from R and have been incorporated back into a more
final Excel file. This final sheet wold be a starting point for further annotation analysis, or other biological followup.

A popular genomics statistical package (edgeR) is used in the analysis of results from both MaxQuant and PAW. An additional basic t-test 
analysis of the MaxQuant data was also done and compared to the edgeR testing.
