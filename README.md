sc-atacseq-explorer
===================

This single cell ATAC-Seq analysis pipeline is designed for intergative analysis of
dataset, initially processed by 10x Genomics [Cell Ranger ATAC][10xcellranger].
Jupyter Notebook format allows you to experiment with any types of parameters on the fly.

In addition to 10x Genomics results it offers:

* Flexible and clear data preprocessing and normalization methods
* Summary on different conditions in case of aggregated dataset
* Different types of clustering followed by heatmap and t-SNE visualizations in low dimensions space
* Top cluster markers visualization on heatmap / t-SNE plot
* Closest genes annotations for peaks and clusters
* Annotated cell-specific genes analysis
* Bigwig and BED files for clusters and markers ready-to-be-visualized in [JBR Genome Browser][jbr] by JetBrains Research
* Data preparation for [single cell explorer][sce] by Artyomov Lab, Washington University in St.Louis

Pipelines
---------
* [Seurat](https://www.biorxiv.org/content/biorxiv/early/2018/11/02/460147.full.pdf)
* [Cell Ranger ATAC](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/algorithms/overview) & [source code](https://github.com/10XGenomics/cellranger)
* [Scasat](https://academic.oup.com/nar/article/47/2/e10/5134327)
* [SnapATAC](https://www.biorxiv.org/content/10.1101/615179v2)
* [Scater](https://academic.oup.com/bioinformatics/article/33/8/1179/2907823)

Materials
---------
* [JetBrains Research BioLabs](https://research.jetbrains.org/groups/biolabs)

[10xcellranger]: https://www.10xgenomics.com/solutions/single-cell-atac/
[jbr]: https://research.jetbrains.org/groups/biolabs/tools/jbr-genome-browser
[sce]: https://artyomovlab.wustl.edu/shiny/single_cell_explorer