[![JetBrains Research](https://jb.gg/badges/research.svg)](https://confluence.jetbrains.com/display/ALL/JetBrains+on+GitHub)

sc-atacseq-explorer
===================
This single cell ATAC-Seq analysis pipeline is designed for advanced analysis of dataset, 
produced by 10X Genomics [Cell Ranger ATAC][10xcellranger]. 
Aggregated datasets are also supported!

In addition to 10x Genomics results it offers:

* Capable to process aggregated data by 10X Genomics Cell Ranger ATAC.
* Summary on different conditions in case of aggregated dataset
* Flexible data processing with t-SNE/UMAP visualizations in low dimensions space
* User defined markers visualization as a heatmap
* Closest genes annotations for peaks and clusters
* Annotated markers analysis
* Bigwig and BED files for clusters and markers ready-to-be-visualized in [JBR Genome Browser][jbr]
* Data preparation for [Single Cell Explorer][sce]
* Save all the figures to ready for publication PDF format

10x Genomic Cell Ranger ATAC
----------------------------
* Launch batch cell ranger processing.\
  **NOTE**: we don't launch it in parallel because of martian framework used by Cell Ranger ATAC.

```
for SAMPLE in $(ls *.fastq.gz | sed -E 's/_S[0-9]_L00.*//g' | sort --unique); do
    cellranger-atac count --id=cra_${SAMPLE} --fastqs=${WORK_DIR} --sample ${SAMPLE} --reference ${REFERENCE}
done
```

* Create aggregation file `merged.csv`

```
library_id,fragments,cells
<id1>,<path1>/outs/fragments.tsv.gz,<path1>/outs/singlecell.csv
...
<idN>,<pathN>/outs/fragments.tsv.gz,<pathN>/outs/singlecell.csv
```

* Launch aggregation

```
cellranger-atac aggr --id=<id> --csv merged.csv --normalize=depth --reference=${REFERENCE}
```

Prerequisites
-------------
Conda environment `sc-atac-explorer` can be easily created for launching Jupyter Notebook:

    ```
    conda env create -f environment.yml
    conda activate sc-atac-explorer
    ```

Pipeline
--------
Launch jupyter notebook to proceed with the pipeline.

    ```
    conda activate sc-atac-explorer
    jupyter notebook
    ```

Other pipelines
---------------
* [Cell Ranger ATAC](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/algorithms/overview)
* [Seurat](https://www.biorxiv.org/content/biorxiv/early/2018/11/02/460147.full.pdf)
* [Scasat](https://academic.oup.com/nar/article/47/2/e10/5134327)
* [SnapATAC](https://www.biorxiv.org/content/10.1101/615179v2)
* [Scater](https://academic.oup.com/bioinformatics/article/33/8/1179/2907823)
* [Signac](https://satijalab.org/signac/articles/pbmc_vignette.html)
* [Benchmark](https://github.com/pinellolab/scATAC-benchmarking)


References
----------
* [Cell Ranger ATAC source code](https://github.com/10XGenomics/cellranger)
* [JetBrains Research BioLabs](https://research.jetbrains.org/groups/biolabs)

[10xcellranger]: https://www.10xgenomics.com/solutions/single-cell-atac/
[jbr]: https://research.jetbrains.org/groups/biolabs/tools/jbr-genome-browser
[sce]: https://artyomovlab.wustl.edu/shiny/single_cell_explorer
