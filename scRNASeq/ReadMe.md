# single cell RNA-Sequencing

This is Github repo containing notes / usage examples for analyzing single-cell RNA-Seq (scRNA-Seq) data using various R/Bioconductor packages.

## basic processing and normalization 

* cellity [[bioC link](https://bioconductor.org/packages/release/bioc/html/cellity.html)] - classification of low quality cells in scRNA-Seq data
* scLVM [[GitHub link](https://github.com/PMBio/scLVM)] - Dissects the observed heterogeneity into different sources, thereby allowing for the correction of confounding sources of variation. scLVM was primarily designed to account for cell-cycle induced variations in single-cell RNA-seq data where cell cycle is the primary source of variability.
* [scater](https://github.com/stephaniehicks/bioconductorNotes/blob/master/scRNASeq/basic_processing/scater.Rmd) [[bioC link](https://www.bioconductor.org/packages/release/bioc/html/scater.html)] - This R/Bioconductor package includes (1) QC metrics, (2) transcript quantification using Kallisto, (3) many functions for data visualization.
* SCONE [[bioC link](http://bioconductor.org/packages/devel/bioc/html/scone.html)] - Stands for 'Single-Cell Overview of Normalized Expression'. Includes quality control metrics and normalization. 
* [scran](https://github.com/stephaniehicks/bioconductorNotes/blob/master/scRNASeq/basic_processing/scran.Rmd) [[bioC link](http://bioconductor.org/packages/release/bioc/html/scran.html)] - Basic analyses for scRNA-Seq data. This R/Bioconductor package (1) estimates [pool-based and spike-in based normalization size factors](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0947-7), (2) assigns cells to cell cycle phases, (3) detects highly variable and correlated genes across cells.


## distance, clustering and tree structures 

* clusterExperiment [[bioC devel link](http://bioconductor.org/packages/devel/bioc/html/clusterExperiment.html)] - Functions for running and comparing many different clusterings of single-cell sequencing data. Meant to work with SCONE and slingshot
* SC3 [[bioC link](https://bioconductor.org/packages/release/bioc/html/SC3.html)] - unsupervised clustering of cells
* [sincell](https://github.com/stephaniehicks/bioconductorNotes/blob/master/scRNASeq/distance_clustering/sincell.Rmd) [[bioC link](http://bioconductor.org/packages/release/bioc/html/sincell.html)] - Cell-state similarity and hierarchy for single-cell RNA-Seq data. This R/Bioconductor package (1) assesses cell-to-cell distance and (2) builds a cell-state hierarchy. Assessments of the robustness of the cell-state hierarchy are also provided.

See Seurat below too. 

## cells ordered by time or space

* cellTree [[bioC link](http://bioconductor.org/packages/release/bioc/html/cellTree.html)] - Uses Latent Dirichlet Allocation to build tree structures connecting cells
* [monocle](https://github.com/stephaniehicks/bioconductorNotes/blob/master/scRNASeq/ordered_pseudotime/monocle.Rmd) [[bioC link](https://bioconductor.org/packages/release/bioc/html/monocle.html)] - time-series analysis (but also does differential expression)
* slingshot [[GitHub link](https://github.com/kstreet13/slingshot)] - functions for identifying and characterizing continuous developmental trajectories in single-cell sequencing data
* TSCAN [[GitHub link](https://github.com/zji90/TSCAN)] - Pseudo-time reconstruction and evaluation in scRNA-Seq data
* SPRING [[GitHub link](https://github.com/AllonKleinLab/SPRING/)] - data viz tool for differentiation trajectories. 
    * Allows user to see full complex cell states (similar to using t-SNE on HGV) AND shows user trajectories between states (similar to using diffusion maps)



## differential expression and/or gene set enrichment analysis

* BASiCs [[GitHub link](https://github.com/catavallejos/BASiCS)] - Stands for 'Bayesian Analysis of Single-Cell Sequencing Data'. This R package (1) estimates cell-specific normalization constants as part of the model parameters, (2) unexplained technical variability is quantified based on spike-in genes, (3) identifies high/lowly expressed genes and total variability of the expression counts is decomposed into technical and biological components. Updates to the package can now identify genes that undergo changes in cell-to-cell heterogeneity (beyond comparisons of over all means).
* MAST [[GitHub link](https://github.com/RGLab/MAST)] - Stands for 'Model-based Analysis of Single-cell Transcriptomics'. This R package (1) filters low-quality cells, (2) adaptive thresholding of background noise, (3) tests for univariate differential expression (with adjustment for covariates and cellular detection rate), (4) gene set enrichment analysis, (5) exploration of gene-gene correlations and co-expression. Can be applied to multiplexed qPCR, NanoString, RNA-Seq. Based on a Hurdle model - two part generalized to model expression and can include biological and technical covariates. Borrows strength across genes to estimate gene variances (empirical Bayes)
* scDD [[GitHub link](https://github.com/kdkorthauer/scDD)] - This R package contains a Bayesian method to detect differences in multi-modal gene expression distributions within and between biological conditions. Calculates a Bayes Factor score based on the Product Partition Model. To evaluate significance, condition labels are permuted and new Bayes Factor scores are calculated. Can include covariates by permuting residuals of a linear model including the covariate and using fitted values (but cannot directly incorporate in model like MAST). Can adjust for cellular detection rate (CDR) from Finak et al. if a potential confounder variable.
* SCDE [[bioC link](https://www.bioconductor.org/packages/release/bioc/html/scde.html)] - This R package includes (1) tests for differential expression for two conditions and (2) pathway and gene set overdisperson analysis to identify and characterize putative cell subpopulations (PAGODA). Treats zero-counts a technical dropouts (but impossible to know if they are biological or technical). Method uses a mixture model with a mixing parameter to account for the dropout events, and a location parameter to estimate the mean expression level of expressed cells. A test is provided to detect mean differences of expressed cells, but no guideline is provided to characterize differences due to percentage of zero changes, and extra heterogeneity among expressed cells that gives rise to distributions that are not unimodal is not accommodated (see scDD instead). Cannot adjust for biological covariates or batch/time information. 
* SCPattern [[GitHub link](https://github.com/lengning/SCPattern)] - This R packages uses empirical Bayes to identify differentially expressed genes for ordered scRNA-Seq experiments (e.g. conditions, time, space) vs just two conditions. Considers both zeros and non-zero counts. Based on Kolmogorov-Smirnov statistic. Also classifies expression ?patterns? with probability estimates. Ordered conditions can allows inference for whether zero counts are biological or technical (e.g. gradual changes in percentage of zeros as conditions change may be biological signal vs constant percentage of zeros).

See Monocle, Seurat too.


## spatial reconstruction

* Seurat [[R package link](http://www.satijalab.org/seurat.html)] - This R package can (1) identify highly variable genes, (2) perform dimensionality reduction (PCA, ICA, t-SNE), (3) perform standard unsupervised clustering algorithms (density clustering, hierarchical clustering, k-means) for cell types and states, (4) tests for differentially expressed genes and (5) spatial reconstruction 


## Experimental design and power analyses

* Similar to bulk RNA-Seq, it is important to randomize (or balance) design across as many factors as possible to not introduce artifacts imposed from sample preparation or data collection
	* Potential limiting factors in experimental design that must be balanced with experimental design: 
		* Protocols to isolate cells [see [Saliba et al. (2014)](http://www.ncbi.nlm.nih.gov/pubmed/25053837) and [Kolodziejczyk et al. (2015)](http://www.ncbi.nlm.nih.gov/pubmed/26000846)]
		* Protocols to extract mRNA from isolate cells [see [Grun and van Oudenaarden (2015)](http://www.ncbi.nlm.nih.gov/pubmed/26544934), [Saliba et al. (2014)](http://www.ncbi.nlm.nih.gov/pubmed/25053837) and [Kolodziejczyk et al. (2015)](http://www.ncbi.nlm.nih.gov/pubmed/26000846)]
		* Protocols for amplification methods [see [Dueck et al. (2016)](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-3300-3)]
		* Choice to use synthetic spike-ins or not (or unique molecular identifiers UMIs) [see [Stegle et al. (2015)](http://www.ncbi.nlm.nih.gov/pubmed/25628217)]
			* Spike-ins can take up a HUGE percentage of "read landscape" in sequencing
			* UMIs can remove amplification bias, but are 5' or 3' end biased (bad for isoform or allele-specific expression)
		* Choice of sequencing platform, time and budget constraints. 
* Trade-off between number of cells and sequencing depth
	* Not many papers on determining min number of cells to sequence for a given biological question
	* Most papers focus on min sequencing depths to observe majority of expressed genes (~ 500k to 1 million reads) [see [Shalek et al. (2014)](http://www.ncbi.nlm.nih.gov/pubmed/24919153) and [Tung et al. (2016)](http://giphy.com/gifs/barack-obama-president-they-tried-TEFplLVRDMWBi)]
		

#### Relevant papers 

* [Tung et al. (2016)](http://biorxiv.org/content/early/2016/07/08/062919) - Batch effects and the effective design of single-cell gene expression studies
* [Bacher and Kendziorski (2016)](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0927-y) - Design and computational analysis of single-cell RNA-sequencing experiments
* [Hicks et al. (2015)](http://biorxiv.org/content/early/2015/12/27/025528) - On the widespread and critical impact of systematic bias and batch effects in single-cell RNA-Seq data
* [Grun and van Oudenaarden (2015)](http://www.ncbi.nlm.nih.gov/pubmed/26544934) - Design and Analysis of Single-Cell Sequencing Experiments
* [Svensson et al. (2016)](http://biorxiv.org/content/early/2016/09/08/073692) - Power Analysis of Single Cell RNA?Sequencing Experiments

## Other resources

* [List of software packages for single-cell data analysis from Sean Davis](https://github.com/seandavi/awesome-single-cell/blob/master/README.md) - including RNA-Seq, ATAC-Seq, etc
* [List of software packages and papers for single-cell data analysis from Ming Tang](https://github.com/crazyhottommy/RNA-seq-analysis#single-cell-rna-seq)
* Notes/Slides from [Bioconductor 2016 Conference](http://bioconductor.org/help/course-materials/2016/BioC2016/) on analyzing scRNA-Seq data
	* [Slides from Sandrine Dudoit](http://bioconductor.org/help/course-materials/2016/BioC2016/InvitedTalks1/160624-Dudoit-scrnaseq.pdf) - discusses upcoming BioC packages including Scone, clusterExperiment and slingshot
	* [Workshop on single-cell differential expression and gene set enrichment with MAST from Andrew McDavid](http://bioconductor.org/help/course-materials/2016/BioC2016/ConcurrentWorkshops2/McDavid.html)
* [Bioconductor Workflow for basic analyses of scRNA-Seq data](http://bioconductor.org/help/workflows/simpleSingleCell/)
* [RNA-seq analysis resources from crazyhottommy](https://github.com/crazyhottommy/RNA-seq-analysis#single-cell-rna-seq) - Very broad list that includes some single cell RNA-seq packages and papers.
* [Analysis of single-cell RNA-seq data from Hemberg Lab](http://hemberg-lab.github.io/scRNA.seq.course/index.html) - discusses scater R/Bioconductor package, expression quantification with UMIs, expression QC, data visualizations, confounding, clustering, identifying 'important genes', ordering cells in pseudo-time and differential expression. 
* [Harvard STEM Cell Institute Single Cell Workshop 2015](http://hms-dbmi.github.io/scw/) - workshop on common computational analysis techniques for scRNA-seq data from differential expression to subpopulation identification and network analysis.


