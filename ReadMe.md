# bioconductorNotes

This is Github repo containing notes / usage examples for various R/Bioconductor packages.


## single cell RNA-Seq (scRNA-Seq)

#### basic processing and normalization 

* [scran](https://github.com/stephaniehicks/bioconductorNotes/blob/master/scRNASeq/basic_processing/scran.Rmd) - Basic analyses for single-cell RNA-Seq data. This R/Bioconductor package (1) estimates pool-based and spike-in based normalization size factors, (2) assigns cells to cell cycle phases, (3) detects highly variable and correlated genes across cells.


#### distance, clustering and tree structures 

* [sincell](https://github.com/stephaniehicks/bioconductorNotes/blob/master/scRNASeq/distance_clustering/sincell.Rmd) - Cell-state similarity and hierarchy for single-cell RNA-Seq data. This R/Bioconductor package (1) assesses cell-to-cell distance and (2) builds a cell-state hierarchy. Assessments of the robustness of the cell-state hierarchy are also provided.


#### cells ordered by time or space

* [monocle](https://github.com/stephaniehicks/bioconductorNotes/blob/master/scRNASeq/ordered_pseudotime/monocle.Rmd) 
* cellTree [bioC link](http://bioconductor.org/packages/release/bioc/html/cellTree.html) - Uses Latent Dirichlet Allocation to build tree structures connecting cells


#### Notes/Slides from [Bioconductor 2016 Conference](http://bioconductor.org/help/course-materials/2016/BioC2016/) on analyzing scRNA-Seq data

* [Slides from Sandrine Dudoit](http://bioconductor.org/help/course-materials/2016/BioC2016/InvitedTalks1/160624-Dudoit-scrnaseq.pdf) - discusses new BioC packages including Scone (for normalization and expression quantification), clusterExperiment (for resampling-based sequential ensemble clustering) and slingshot (cell lineage and pseudotime inference)
* [Workshop on single-cell differential expression and gene set enrichment with MAST from Andrew McDavid](http://bioconductor.org/help/course-materials/2016/BioC2016/ConcurrentWorkshops2/McDavid.html)


## Other

* [limma](https://github.com/stephaniehicks/bioconductorNotes/blob/master/limma.Rmd) - Analysis of gene expression data from microarrays or RNA-Seq platform technologies
* [minfi](https://github.com/stephaniehicks/bioconductorNotes/blob/master/minfi.Rmd) - Tools for analyzing and visualizing Illumina's 450k array data
* [GEOquery](https://github.com/stephaniehicks/bioconductorNotes/blob/master/GEOquery.Rmd) - Accesses data on [GEO](http://www.ncbi.nlm.nih.gov/geo/)
* [epivizr](https://github.com/stephaniehicks/bioconductorNotes/blob/master/epivizr.Rmd) - Interactive visualization of genomic data in a browser. Supports Bioconductor data structures such as `GenomicRanges` and `SummarizedExperiments`
