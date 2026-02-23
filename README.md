# scitargets <img src="man/figures/logo.png" align="right" height="99"/>

The goal of scitargets is to help you analyse single-cell CITE-seq data by providing target factories and parametized report to deal with basic quality controls, HTO demultiplexing, and celltypes annotation.

## Installation

You can install the development version of scitargets like so:

```{r}
install.package("devtools")
# At the time I write this, seurat is on the version 5.3.0 on cran.
# I force seurat 5.3.1 or higher as the 5.3.0 version is bugged when using the azimuth package
devtools::install_github("satijalab/seurat", "fix/v.5.3.1", force=TRUE)
devtools::install_github("https://github.com/Pierre9344/scitargets", build_vignettes = T)
```

## Example

1.  In Rstudio, create a new project using scitargets project template (File -\> New Project -\> New Directory -\> Project template for scitargets).

2.  Copy your cellranger output inside **"data/cellranger_output/"** folder. Create a subfolder for each run.

3.  Open the \*\*"\_targets.R"\*\* script and modify it to run the pipeline on each run.

    -   Duplicate the **"tar_demultiplex_hto"** function to have one call per run and modify the argument (the run_id argument must correspond the subfolder in "data/cellranger_output/".

4.  (Optional), unquote the "**crew**" code in "\_targets.R" to allow a parallelization of the pipeline steps and accelerate the pipeline).

5.  Duplicate the "**.qmd**" reports (QMD folder) so that you have one document per run.

6.  Modify the document to fit your run:

    -   set the RUN_ID parameter (inside the yaml header) to the name of the run folder inside the **data/cellranger_output/** folder.

    -   replace all occurence of **"RUN_ID"** by **"\<your_run_id\>"** (same as the RUN_ID parameter in the yaml header).

7.  Modify the \*\*"\_quarto.yml"\*\* file to indicate each reports (qmd document) inside the website

8.  Run the pipeline using **"targets::tar_make()"**


```{r}
library(targets)
tar_make()
```

Check the vignettes for more details.

## Azimuth annotation

Build-in Azimuth annotation are currently only available for the (human) pbmc-ref of the SeuratData package.

It is also recommanded to use:

-  Seurat 5.3.1.1000
-  SeuratObject 5.2.0
-  Azimuth 0.5.0

Other versions may be compatible but not tested