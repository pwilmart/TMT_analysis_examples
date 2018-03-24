# Anaysis of a dilution series
**Click on the Jupyter notebook file (_MAN1353_peptides_proteins.ipynb_) to see the notebook in your browser. It may take a minute for the page to render and load, so please be patient.**

The sample is a mouse brain prep that was digested, split into 6 aliquots, digested, and TMT labeled. A mixture was created with 
relative volumes of 25 to 20 to 15 to 10 to 5 to 2.5. The mixture was run on a Thermo Fusion using the synchronous precursor scan 
MS3 method:

> McAlister, G.C., Nusinow, D.P., Jedrychowski, M.P., Wühr, M., Huttlin, E.L., Erickson, B.K., Rad, R., Haas, W. and Gygi, S.P., 2014. MultiNotch MS3 enables accurate, sensitive, and multiplexed detection of differential expression across cancer cell line proteomes. Analytical chemistry, 86(14), pp.7150-7158.

The data were processed with MSConvert from Proteowizard to extract the MS2 and MS3 scans. The MS2 scans were processed with Comet
to identify the peptides. Accurate mass and target/decoy methods were used to filter PSMs to a 1% FDR. Parsimonious protein inference
and protein grouping were used to create protein and peptide reports. Shared or unique peptide status was determined based on the final 
list of identified proteins. Unique peptides were used for quantification.

> Chambers, M.C., Maclean, B., Burke, R., Amodei, D., Ruderman, D.L., Neumann, S., Gatto, L., Fischer, B., Pratt, B., Egertson, J. 
and Hoff, K., 2012. A cross-platform toolkit for mass spectrometry and proteomics. Nature biotechnology, 30(10), p.918.

> Eng, J.K., Jahan, T.A. and Hoopmann, M.R., 2013. Comet: an open‐source MS/MS sequence database search tool. 
Proteomics, 13(1), pp.22-24.

> Wilmarth, P.A., Riviere, M.A. and David, L.L., 2009. Techniques for accurate protein identification in shotgun proteomic studies of 
human, mouse, bovine, and chicken lenses. Journal of ocular biology, diseases, and informatics, 2(4), pp.223-234.

The analysis looks at the properties (mostly variance) of the reporter ion data at three levels: PSMs, peptides, and proteins. The
PSMs are the raw MS3 scan data. Sequence identities are not known. The PSMs will include non-identifiable scans, contaminants, 
decoys, incorrect scans, and correct scans. The peptides are much cleaner. They are derived from 1% FDR PSMs, contaminants and decoys 
have been removed, and shared peptides have been excluded. Peptides also have some degree of aggregation. Multipule MS2 scans and 
different charge states were combined.

Reporter ions for proteins are aggregated (summed) from the reporter ions of all constituent peptides. Only the unique peptides are
used in the sums. There is a great reduction in the number of proteins compared to the number of peptides or PSMs. This greatly 
simplifies analyses. The analysis also demonstrates that the aggregation improves the consistency of the reporter ion measurements.

An interesting aside is looking at the normalization functions in edgeR. edgeR is a useful RNA-Seq statistical package that I have
used in other work in this and related repositories. The normalization factors that edgeR reports can seem confusing and they are
explained here.

> Robinson, M.D., McCarthy, D.J. and Smyth, G.K., 2010. edgeR: a Bioconductor package for differential expression analysis of 
digital gene expression data. Bioinformatics, 26(1), pp.139-140.

> Robinson, M.D. and Oshlack, A., 2010. A scaling normalization method for differential expression analysis of RNA-seq data. 
Genome biology, 11(3), p.R25.
