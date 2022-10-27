# GWAS_pipeline
Within this pipeline, the following software is used:

• PLINK v1.90p 64-bit (2 Jun 2015) gcta-1.94.1

• R version 4.1.2

R packages used can be installed using conda:

library(plyr)  
library(ggplot2)  
library(ggpubr)  
library(broom)  
library(chron)

------------------------------------------------------------------------

-   This document describes what each script and step of the GWAS
    analysis pipeline does. For detailed commands, check out the
    **all_commands.sh** document.  
-   There are 2 main checkpoints in the pipeline:  

1.  **checkpoint 1**: After processing phenotypes, look at the
    percentage of variance explained by different covariates. If
    everything looks ok, proceed to the next step.
2.  **checkpoint 2**: If the SNP heritability estimates look ok, submit
    the next steps.

------------------------------------------------------------------------

The computational requirements for each step are based on a project with
N = 2,300, unpruned set of SNPs (\~6M SNPs) and 90 traits.

------------------------------------------------------------------------

1.  **Clone the git repo**

------------------------------------------------------------------------

2.  **Make all files in the code directory executable**  
    This is the directory that contains the GWAS scripts

``` bash
chmod +x /projects/ps-palmer/apurva/genetic_analysis/code

```

------------------------------------------------------------------------

3.  **Create directory structure**  
    This will create the directory structure for the project
    *p50_david_dietz_2020*

``` bash
/projects/ps-palmer/apurva/genetic_analysis/code/directory_structure.sh p50_david_dietz_2020
```

------------------------------------------------------------------------

4.  **Copy csv file with raw phenotypes and trait descriptions to
    data/raw_data**  

-   Ideally, we’d like to query the **gwas_phenotypes** and the
    **descriptions** tables from the database for GWAS analysis. We’ll
    need to upload the tables in the database regularly to make this
    possible. For now, we’ll copy these files to the
    **project_name/data/raw_data** folder.  
-   The project metadata \[PI name, project title etc\] should also be
    queried from the database. This information is in the
    **project_metadata** table in the sample_tracking schema. I’m
    storing this information in the **traits.RData** object for now.  
-   The **process_phenotypes.R** script will look for *raw_data.csv* and
    the *traits.RData* (data dictionary) in the
    **project_name/data/raw_data** folder.  

Here’s what the **process_phenotypes.R** script does:  
- removes white spaces in any trait or covariate values and sets them to
NA  
- quantile normalize the traits within each sex and center/project_name
value. (This is similar to using sex as a covariate)

There are 2 sets of covariates:  
- common covariates : coat color, cohort number. These are not specific
to an experiment.  
- covariates specific to an experiment : experiment box, age at
experiment. These should have the following naming convention: append
the experiment string in front of the covariate. For example, “loco_age”
will be used as a covariate for these locomotor traits: “loco_center”,
“loco_rear”, “loco_distance”, “loco_activity”

*output: data/residuals/residuals.RData*

**Before proceeding to the next steps, check the PVE explained by the
different covariates. These RData frames are saved in
pheno_processing_summary/covs folder**

------------------------------------------------------------------------

5.  **Subset genotype data and calculate GRM, SNP h2 estimates**  

-   *genotype_filters_on_vcf.sh* This will prepare the genotype data for
    GWAS. If we’re planning to get the PLINK binary file as an output of
    the genotyping pipeline, then we’ll only need to subset the file to
    the project rfids. This can be done using the *–keep* option in
    PLINK.  
-   *create_famid.R* This will create the sample id list. This list is
    an input to the –keep command for *genotype_filters_on_vcf.sh*
-   **Calculate 2 sets of GRMs**:  
-   *all_chromosomes_grm.sh*: GRM estimated from all SNPs on all
    chromosomes.  
-   *subtract_grm.sh*: calculated from SNPs on one chromosome. This GRM
    will be used for *–mlma-subtract-grm* option designed to parallelise
    the MLMA-LOCO analysis for large data set.  
    *More info*:
    <https://gcta.freeforums.net/thread/173/mixed-linear-model-association-analysis>   
-  *LZ_prep.sh* detailed steps for preparing the locuszoom backend
    SQLite database are outlined in this document:
    <https://www.dropbox.com/s/2ugc710e2sp9v98/locuszoom_standalone_instructions.pdf?dl=0>

``` bash
Rscript create_famid.R 

qsub -q condo -W depend=afterok:$geno_filters -l nodes=1:ppn=6 -l walltime=8:00:00  /projects/ps-palmer/apurva/genetic_analysis/code/extract_dosages.sh

geno_filters=`qsub -q condo -l nodes=1:ppn=6 -l walltime=3:00:00  /projects/ps-palmer/apurva/genetic_analysis/code/genotype_filters_on_vcf.sh`  

qsub -q condo -W depend=afterok:$geno_filters -l nodes=1:ppn=6 -l walltime=3:00:00  /projects/ps-palmer/apurva/genetic_analysis/code//subtract_grm.sh  

all_grm=`qsub -q condo -W depend=afterok:$geno_filters -l nodes=1:ppn=6 -l walltime=3:00:00  /projects/ps-palmer/apurva/genetic_analysis/code/all_chromosomes_grm.sh` 

qsub -q condo -W depend=afterok:$all_grm -l nodes=1:ppn=2 -l walltime=1:00:00 /projects/ps-palmer/apurva/genetic_analysis/code/snp_heritability.sh
```

------------------------------------------------------------------------
