#cd to project directory before executing all the steps
#if there is no qsub command in front of the script, it was not a computationally intensive step and I did not queue it on the cluster. This can cause unnecessary delays. 

#Step 1 
#Process phenotypes 

Rscript process_phenotypes.R


#########################################################################################


#check the PVE explained by covariates. These RData frames are saved in pheno_processing_summary/covs folder
#if everything looks ok, proceed to the next steps

#########################################################################################


#Step 2
#Prep genotype data, extract dosages for conditional analysis, estimate GRMs, prepare locuszoom standalone database for regional association plots
#These can be run as bash job dependencies:

Rscript create_famid.R 
qsub -q condo -W depend=afterok:$geno_filters -l nodes=1:ppn=6 -l walltime=8:00:00  /projects/ps-palmer/apurva/genetic_analysis/code/extract_dosages.sh
geno_filters=`qsub -q condo -l nodes=1:ppn=6 -l walltime=3:00:00  /projects/ps-palmer/apurva/genetic_analysis/code/genotype_filters_on_vcf.sh`  
qsub -q condo -W depend=afterok:$geno_filters -l nodes=1:ppn=6 -l walltime=3:00:00  /projects/ps-palmer/apurva/genetic_analysis/code//subtract_grm.sh  
all_grm=`qsub -q condo -W depend=afterok:$geno_filters -l nodes=1:ppn=6 -l walltime=3:00:00  /projects/ps-palmer/apurva/genetic_analysis/code/all_chromosomes_grm.sh` 
qsub -q condo -W depend=afterok:$all_grm -l nodes=1:ppn=2 -l walltime=1:00:00 /projects/ps-palmer/apurva/genetic_analysis/code/snp_heritability.sh


#locuszoom standalone database
#This is a series of bash and R commands. It also requires editing conf/m2zfast.conf, and changing the SQLITE_DB variable. For more information check out the locuszoom standalone instructions PDF. 
LZ_prep.sh


#########################################################################################


#check the SNP heritability estimates. If these are close to zero, it indicates some issue such as data being scrambled or some issue with how the traits were calculated. 
#if everything looks ok, proceed to the next steps

#########################################################################################


#Step 3

Rscript gwas_jobs_files.R

#depending on the number of traits you're analyzing, you'll need to wait before you submit all the array jobs because you'll hit the max jobs in the queue limit. if less than 75 traits, 
Rscript qsub_gwas_jobs.R


#calling qtls
run using qsub_calling_QTLs.R


#prep for conditional analysis

#conditional


#manhattan


#eqtl


#phewas


#locuszoom


#compress images


#knit report on the cluster


#clean up intermediate files from project directory to free up storage





