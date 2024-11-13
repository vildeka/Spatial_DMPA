# Spatial transcriptomics unveils estrogen-modulated immune responses and structural alterations in the ectocervical mucosa of depot medroxyprogesterone acetate users

Vilde Kaldhusdal, Mathias Franzen Boger, Annelie Tjernlund, Adam D. Burgener, Frideborg Bradley, Julie Lajoie, Kenneth Omollo, Joshua Kimani, Keith Fowke, Paulo Czarnewski, and Kristina Broliden

![](./resources/Graphical%20abstract.png)

## Table of contents

-   [General info](#general-info)
-   [Dependencies](#dependencies)
-   [Data Availability Statment](#data-availability-statment)
-   [Repo description](#repo-description)
-   [Setup](#setup)

## General info {#general-info}

doi:

This repository contains the code related to the article "Spatial transcriptomics unveils estrogen-modulated immune responses and structural alterations in the ectocervical mucosa of depot medroxyprogesterone acetate users"

This project used multiple datasets:

-   Spatial transcriptomics data (10x Visium)

-   Transcriptomics data (bulk mRNA-SEQ)

## Dependencies

Project is created with:

-   R version: 4.1.2
-   RStudio version: 1.4.1106 (use version seperate from conda env.)
-   renv version: 0.15.2
-   Seurat version:

## Data Availability Statement

**Spattial transcriptomics count data** files and RDS object can be accessed in the Gene Expression Omnibus public repository, SuperSeries ID GSE217237. The raw transcriptomic sequencing data cannot be held in a public repository due to the sensitive nature of such personal data. Request for data access can be made to the Karolinska Institutet Research Data Office (contact via rdo\@ki.se), and access will be granted if the request meets the requirements of the data policy.

**Bulk transcriptomics count data** files can be accessed in the Gene Expression Omnibus public repository, SuperSeries ID GSE217237. The raw transcriptomic sequencing data cannot be held in a public repository due to the sensitive nature of such personal data. Request for data access can be made to the Karolinska Institutet Research Data Office (contact via rdo\@ki.se), and access will be granted if the request meets the requirements of the data policy.

## Repo description

-   **src**\
    contains all the analysis script including all preprocessing steps
-   **manuscript**\
    reproducible code for figures included in the manuscript
-   **md_files**\
    rendered versions of the Rmds in src and manuscript folder respectively
-   **bin**\
    R scripts with various helper functions, mostly visualization functions for the ST data
-   **data**\
    Empty folder to store data downloaded from GEO

<!-- -->

```         
project
│   README.md
│   renv.loc    
└───src
│   │   00_Preprocessing.Rmd
│   │   02_Analysis.Rmd
│   │   ...
|   └───md files (rendered versions of the Rmds)
│   |       │   00_Preprocessing.md
│   |       │   02_Analysis.md
│   |       │   ...
|   |
│   └───manuscript
│       │   Figure01.Rmd
│       │   Figure02.Rmd
│       │   ...
|       └───md files
│           │   Figure01.md
│           │   Figure02.md
│           │   ...
└───bin
│   │   file01.txt
│   │   file02.txt
└───data
│   │   file01.txt
│   │   file02.txt
│
```

## Setup

Recommended setup to run this project:

Conda + renv

1.  Clone the repo
2.  If not already installed download mini conda/conda
3.  In the terminal navigate to the project directory
4.  create a new enviroment:<br/> `conda env create -n Spatial_DMPA -f environment.yml`
5.  Activate the enviroment:<br/> `conda activate Spatial_DMPA`
6.  Open Rstudio:<br/> `rstudio& Spatial DMPA`
7.  Install all packages specified by the lockfile:<br/> `renv::restore()`
