sc-atacseq-explorer
===================

This single cell ATAC-Seq analysis pipeline is designed for intergative analysis of
dataset, initially processed by 10x Genomics [Cell Ranger ATAC][10xcellranger].
Jupyter Notebook format allows you to experiment with any types of parameters on the fly.

It provides original 10x Genomics data re-analysis:

* Summary metrics visualization like cell count, total coverage
* Different conditions in case of several samples aggregated by `aggr` command.
* Automated BED files creation for each cluster / condition

In addition to 10x Genomics results it offers:

* Flexible and clear data preprocessing and normalization methods
* Summary on different conditions in case of aggregated dataset
* Different types of clustering followed by heatmap and t-SNE visualizations in low dimensions space
* MannWhitney U test based differential analysis
* Top cluster markers visualization of t-SNE plot
* Closest genes annotations for peaks and clusters
* Annotated cell-specific genes analysis
* BED files for clusters and markers ready-to-be-visualized in [JBR Genome Browser][jbr] by JetBrains Research
* Data preparation for [single cell explorer][sce] by Artyomov Lab, Washington University in St.Louis

Materials
---------
* [JetBrains Research BioLabs](https://research.jetbrains.org/groups/biolabs)
* [Cell Ranger source code](https://github.com/10XGenomics/cellranger)

[10xcellranger]: https://www.10xgenomics.com/solutions/single-cell-atac/
[jbr]: https://research.jetbrains.org/groups/biolabs/tools/jbr-genome-browser
[sce]: https://artyomovlab.wustl.edu/shiny/single_cell_explorer