# An analysis of long-term mortality risks of air pollution in the UK Biobank

An example of survival analysis of long-term risks of air pollution on mortality using a large population-based cohort and time-varying exposure histories

------------------------------------------------------------------------

This repository stores the updated R code to reproduce the analysis presented in the article:

Vanoli J, et al. Long-term associations between time-varying exposure to ambient PM2.5 and mortality: an analysis of the UK Biobank. *Epidemiology*. 2025;36(1):1-10. DOI: 10.1097/EDE.0000000000001796 [[freely available here](http://www.ag-myresearch.com/2025_vanoli_epidemiol.html)]

The repo provides a series of R scripts that can be used to replicate the analysis described in the article (see section R scripts below). The main scripts included in the main folder *Rcode* can be used as a tutorial, and the user can run them to perform a simplified analysis illustrating the key steps using synthetic data. Additional scripts included in sub-folders are added to allow the complete replication of the study, although this would require access to the UKB database (see section Data below).

The work was supported by the Medical Research Council-UK (Grant ID: [MR/Y003330/1](https://gtr.ukri.org/projects?ref=MR%2FY003330%2F1)).

### R scripts

The series of R scripts are provided in the folder *Rcode* and its three sub-folders. They are organised as follows:

-   Main *Rcode* folder: it includes a series of scripts to be used as a tutorial to reproduce the main steps of the analysis, from setting the parameters, loading the synthetic datasets, performing the main and secondary analyses, to creating the main graphs. The results will be similar (but not identical) to the published article.
-   Sub-folder *origcode*: it includes the original scripts used for the analysis. Most of the code can be adapted and used with the synthetic datasets, although the complete replication would require the generation of the original datasets (see below).
-   Sub-folder *ukbcode*: it includes the original scripts that were used to derive a series of real datasets from the UKB database. These scripts cannot be run without having access to the UKB database, and they are included here for transparency and to allow the users to replicate the analysis in full. The real datasets used in the analysis were recreated as synthetic data and used for the tutorial (see below).
-   Sub-folder *synthcode*: it includes the scripts to generate the synthetic datasets (see the section Data below) using the real versions as a template. Similarly to the code in the two sub-folders above, these scripts require access to the original UKB database and can be run after those in *ukbcode* to generate the synthetic datasets.

### Data

The original data from the UKB cohort is available upon request (see [www.ukbiobank.ac.uk](https://www.ukbiobank.ac.uk)) but cannot be shared. For reproducibility purposes, a series of synthetic datasets are stored in Zenodo repository (<https://zenodo.org/records/13983170>), from which they can be downloaded manually or automatically using the R scripts. The datasets resemble the real data used in the analysis, and they were generated using the R package `synthpop` ([www.synthpop.org.uk](www.synthpop.org.uk)). Details on the data synthesis and codebooks for the datasets are provided in the repository.

The series of synthetic datasets are the following:

-   `synthbdcohortinfo`: basic cohort information regarding the follow-up period and birth/death dates for 502,360 participants.
-   `synthbdbasevar`: baseline variables, mostly collected at recruitment.
-   `synthpmdata`: annual average exposure to PM2.5 for each participant reconstructed using their residential history.
-   `synthoutdeath`: death records that occurred during the follow-up with date and ICD-10 code.

**Note**: the synthetic versions of the datasets resemble the real ones in several aspects, including the distribution and cross-correlation of variables as well as some of the effect associations. These can be used to replicate the analysis using the scripts provided in this repo, as well as for illustrative purposes in other tutorials. However, users are warned that **these datasets are fake** and must not be used for making inferences and addressing research hypotheses.
