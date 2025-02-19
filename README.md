# MicroPouch_NP

## Abstract
### background and Aims:
Chronic pouchitis is a common complication after ileal pouch-anal anastomosis (IPAA) with limited treatment options. In this proof-of-concept study, we aimed to investigate safety, efficacy, and microbiome changes associated with using fecal microbiota transplantation (FMT) from a donor with a normal functioning IPAA to induce remission in patients with chronic pouchitis.
### Methods:
The study was a case-series with a 4-week intervention period and 12-month follow-up. Eligible patients with chronic pouchitis were recruited from the Department of Gastrointestinal Surgery, Aalborg University Hospital, Denmark. Participants were treated with FMT derived from a donor with a normal functioning IPAA. Treatment was delivered daily by enema for two weeks followed by every second day for two weeks. Disease severity and quality of life (QoL) were accessed at inclusion and 30-day follow-up. Pouchitis Disease Activity Index (PDAI) <7 was considered equivalent to clinical remission. Fecal samples from participants, healthy donors, and IPAA donor were analyzed by shotgun metagenomic sequencing.
### Results:
Three patients with chronic pouchitis were included and completed the treatment protocol and follow-up visits. At the 30-day follow-up, all participants achieved clinical remission with reduced endoscopic inflammation. The median total PDAI score decreased from 8 (range 10–8) at baseline to 6 (range 6–5) at 30 days. Two participants reported improved QoL, while one reported no change. Few mild, self-limited adverse events were reported by all participants during treatment, with no serious events. Principal component analysis of fecal samples distinguished two clusters: healthy donors and the IPAA donor, with participant samples forming a separate cluster.
### Conclusion:
We observed that all participants achieved clinical remission with reduced endoscopic inflammation following a 4-week FMT intervention. Adverse events were mild and self-limited. Metagenomic analysis revealed distinct microbiome clusters between IPAA donor and recipients, both of which differed from those of healthy donors.

<!-- ## Cite -->

## Microbiome Analysis

This repository contains the code used for generating the microbiome analysis. All scripts for code generation can be found in `src` folder.

An external library is needed, which can be installed as:

`remotes::install_github("SebastianDall/mplibrary")`

Scripts are:
- `src/heatmaps.Rmd`: Produces heatmaps of the relative abundance of the top most abundant species.
- `src/alpha-beta_diversity.Rmd`: Produces alpha and beta diversity plots.
- `src/ordinationplot.Rmd`: Produces ordination plots.

Richness was defined as species with a relative abundance >0 and alpha diversity was calculated using the Shannon diversity index. Patient sample similarity to donors were calculated using both Sørensen coefficient on relative abundances and Bray-Curtis on Hellinger transformed relative abundances. The similarity would be measured as the median similarity to donor samples received and also as median similarity to a donor sample from each donor not received. For the placebo group, similarity was calculated as the median similarity to a donor sample from each donor used in the FMT group.


## Prerequisites
All data was analyzed using R (4.1.0) and RStudio.

## Data
Sequencing data is available at [PRJEB66493](https://www.ebi.ac.uk/ena/browser/view/PRJEB66493) & [PRJEB80556](https://www.ebi.ac.uk/ena/browser/view/PRJEB80556). To rerun the workflow install snakemake and run the snakemake pipeline as:

