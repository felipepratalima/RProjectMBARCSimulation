R Project MBARC Simulation
================
felipeprata
Thu Mar 8 18:41:42 2018

``` r
knitr::opts_chunk$set(echo = TRUE)
```

R Project MBARC Simulation is a R Project to develop multiple meta-omics (metagenomics shotgun, 16S amplicons and metatrasncriptomics) simulated dataset based on the MBARC-26 mock (<https://www.nature.com/articles/sdata201681/>). It was developed by Felipe Prata Lima (<http://lbi.usp.br/membros/>) and João Carlos Setubal (<http://www.iq.usp.br/setubal/>).

Following, we describe the steps we took to this simulation.

Preliminaries
-------------

-   We have been using the RStudio (<https://www.rstudio.com/>) to run our R code -- we recommend the same with the creation of a RStudio Project (<https://support.rstudio.com/hc/en-us/articles/200526207-Using-Projects>).
-   We have been using the Grinder (<https://sourceforge.net/projects/biogrinder/>) to run our reads simulation (single-end Illumina reads with 250 bp with +/- 50 bp standard deviation).
-   For 16S amplicons simulation we have been using the 515F-806R primers (<http://press.igsb.anl.gov/earthmicrobiome/protocols-and-standards/16s/>), which is designed to the ampliciation of the V4 region of the 16S rRNA gene.
-   Optionally, when CDSs files are not available we can predict then with Prodigal (<https://github.com/hyattpd/Prodigal>).

Load the required libraries
---------------------------

Load the following third-party libraries:

``` r
library(magrittr)
library(stringr)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(tools)
library(rentrez)
library(poweRlaw)
```

Load our libraries:

-   Helper functions to retrieve NCBI's information and data, and generate required files for the simulation:

``` r
source("../s4imosFunctions.R")
```

-   And helper methods to deal with NCBI's taxonomy (taxdump files):

``` r
## Install:
# devtools::install_github("felipepratalima/taxdumpr")
## Load:
library(taxdumpr)
## Instantiate an 
x <- Taxdumpr("../ncbi/nodes.dmp", "../ncbi/names.dmp")
```

Folder and programs settings
----------------------------

``` r
# sedwd("/path/to/RProjectMBARCSimulation/")

# PROJECT_FOLDER <- "/path/to/RProjectMBARCSimulation" %>% path.expand
# GENOMES_FOLDER <- "/path/to/RProjectMBARCSimulation/genomes" %>% path.expand
# CDSS_FOLDER <- "/path/to/RProjectMBARCSimulation/cdss" %>% path.expand
# PRODIGAL_LOCATION = "/path/to/prodigal_executable"

# SH_READS_NUMBER <- "1000000"
# MT_READS_NUMBER <- "1000000"
# AMP_READS_NUMBER <- "100000"
# AMP_CHIMERA_PERCENTUAL <- "10"

# READ_SIZE <- "250"
# READ_SIZE_STANDARD_DEVIATION <- "50"

## Amplicon settings
# AMPLICONS_SEARCH_SCRIPT_FILE <- "/path/to/s4imos_amplicon_search.pl"
# PRIMERS_FILE <- "/path/to/515f_806r.fasta" %>% path.expand

## Create folders
# dir.create(PROJECT_FOLDER, showWarnings = F)
# dir.create(GENOMES_FOLDER, showWarnings = F)
# dir.create(CDSS_FOLDER, showWarnings = F)
```

Prepare data
------------

Load the composition file:

``` r
composition.df <- read.delim("RProjectMBARCSimulation.csv")
```

That's the original data:

``` r
knitr::kable(composition.df)
```

| species                       |  numberOfCells|
|:------------------------------|--------------:|
| Clostridium perfringens       |             39|
| Clostridium thermocellum      |             15|
| Coraliomargarita akajimensis  |            144|
| Corynebacterium glutamicum    |             10|
| Desulfosporosinus acidiphilus |            409|
| Desulfosporosinus meridiei    |            561|
| Desulfotomaculum gibsoniae    |            535|
| Echinicola vietnamensis       |             41|
| Escherichia coli              |             16|
| Fervidobacterium pennivorans  |            672|
| Frateuria aurantia            |            317|
| Halovivax ruber               |            614|
| Hirschia baltica              |            400|
| Meiothermus silvanus          |            213|
| Natronobacterium gregoryi     |            569|
| Natronococcus occultus        |            933|
| Nocardiopsis dassonvillei     |              6|
| Olsenella uli                 |            304|
| Pseudomonas stutzeri          |            164|
| Salmonella bongori            |             31|
| Salmonella enterica           |             40|
| Segniliparus rotundus         |            149|
| Spirochaeta smaragdinae       |            467|
| Streptococcus pyogenes        |             16|
| Terriglobus roseus            |            155|
| Thermobacillus composti       |              7|

Adjust the relative abundance:

``` r
composition.df$relativeAbundance <- composition.df$numberOfCells / sum(composition.df$numberOfCells) * 100
```

Specify a search term which is used to locate the genome in NCBI's RefSeq or Assembly (we use rentrez to locate the genomes). Our pattern is to locate a genome by the taxonomy id, using the string txid\[TAXONOMY\_ID\_NUMBER\]. For example, by txid54736, where 54736 is the taxonomy id for the Salmonella bongori species, we will be able to retrieve information and data for Salmonella bongori serovar 66:z41:- str. SA19983605, a strain of this species with genome sequence stored in NCBI's genome (<https://www.ncbi.nlm.nih.gov/genome/?term=txid1243617>).

``` r
composition.df$searchTerm <- taxdumpr::getStandardTaxonomyIdsByNames(x, composition.df$species) %>% paste0("txid", .)
```

Remove duplications (different rows with same search term) when present. When it happens, relative abundances are merged (it is not the case here):

``` r
composition.df <- removeDuplicated(composition.df)
```

Check the composition (the sum of relative abundances should be equals to 100%):

``` r
composition.df <- checkComposition(composition.df)
```

And that's our base simulation data:

``` r
knitr::kable(composition.df)
```

| searchTerm |  relativeAbundance| downloadStatus |
|:-----------|------------------:|:---------------|
| txid1502   |          0.5712612| FALSE          |
| txid1515   |          0.2197158| FALSE          |
| txid395922 |          2.1092720| FALSE          |
| txid1718   |          0.1464772| FALSE          |
| txid885581 |          5.9909184| FALSE          |
| txid79209  |          8.2173722| FALSE          |
| txid102134 |          7.8365314| FALSE          |
| txid390884 |          0.6005566| FALSE          |
| txid562    |          0.2343636| FALSE          |
| txid93466  |          9.8432694| FALSE          |
| txid81475  |          4.6433280| FALSE          |
| txid387341 |          8.9937015| FALSE          |
| txid2724   |          5.8590889| FALSE          |
| txid52022  |          3.1199648| FALSE          |
| txid44930  |          8.3345540| FALSE          |
| txid29288  |         13.6663249| FALSE          |
| txid2014   |          0.0878863| FALSE          |
| txid133926 |          4.4529076| FALSE          |
| txid316    |          2.4022265| FALSE          |
| txid54736  |          0.4540794| FALSE          |
| txid28901  |          0.5859089| FALSE          |
| txid286802 |          2.1825106| FALSE          |
| txid55206  |          6.8404863| FALSE          |
| txid1314   |          0.2343636| FALSE          |
| txid392734 |          2.2703970| FALSE          |
| txid377615 |          0.1025341| FALSE          |

Obtain Genomes and CDS's
------------------------

The follow piece of code, retrieves the genomes and CDS's sequences which will be used for our simulation. We could avoid the FOR statement and simply use the getGenomesAndCdss function, but to have more control (because of some issues we had while running the code in the RStudio, e.g. instability with some programming statements) we choose to run this code this way.

``` r
# aux.df <- NULL
# for (i in 1:nrow(composition.df)) {
#   if (aux.df %>% is.null) {
#     aux.df <- getGenomesAndCdss(GENOMES_FOLDER, CDSS_FOLDER, composition.df[i,], FALSE)
#   } else {
#     aux.df <- rbind(
#       aux.df,
#       getGenomesAndCdss(GENOMES_FOLDER, CDSS_FOLDER, composition.df[i,], FALSE)
#     )
#   }
# }
```

Once this piece code requires a very long time to run, we saved its results (on our original running) and loaded it to help us with this documentation:

``` r
aux.df <- readRDS("aux.df.RDS")

composition.df <- aux.df

knitr::kable(composition.df)
```

| searchTerm |  relativeAbundance| downloadStatus | assemblyname                            | taxid   | organism                                                                          | speciestaxid | speciesname                     | ftppath\_refseq                                                                                                  | ftppath\_genbank                                                                                                 | usedFtpPathForGenome                                                                                             | usedFtpPathForCdss                                                                                               |
|:-----------|------------------:|:---------------|:----------------------------------------|:--------|:----------------------------------------------------------------------------------|:-------------|:--------------------------------|:-----------------------------------------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------|
| txid1502   |          0.5712612| TRUE           | ASM1328v1                               | 195103  | Clostridium perfringens ATCC 13124 (firmicutes)                                   | 1502         | Clostridium perfringens         | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/285/GCF_000013285.1_ASM1328v1>                               | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/013/285/GCA_000013285.1_ASM1328v1>                               | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/013/285/GCA_000013285.1_ASM1328v1>                               | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/013/285/GCA_000013285.1_ASM1328v1>                               |
| txid1515   |          0.2197158| TRUE           | ASM1586v1                               | 203119  | Ruminiclostridium thermocellum ATCC 27405 (firmicutes)                            | 1515         | Ruminiclostridium thermocellum  | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/015/865/GCF_000015865.1_ASM1586v1>                               | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/015/865/GCA_000015865.1_ASM1586v1>                               | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/015/865/GCA_000015865.1_ASM1586v1>                               | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/015/865/GCA_000015865.1_ASM1586v1>                               |
| txid395922 |          2.1092720| TRUE           | ASM2590v1                               | 583355  | Coraliomargarita akajimensis DSM 45221 (verrucomicrobia)                          | 395922       | Coraliomargarita akajimensis    | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/025/905/GCF_000025905.1_ASM2590v1>                               | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/025/905/GCA_000025905.1_ASM2590v1>                               | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/025/905/GCA_000025905.1_ASM2590v1>                               | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/025/905/GCA_000025905.1_ASM2590v1>                               |
| txid1718   |          0.1464772| TRUE           | ASM1132v1                               | 196627  | Corynebacterium glutamicum ATCC 13032 (high GC Gram+)                             | 1718         | Corynebacterium glutamicum      | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/011/325/GCF_000011325.1_ASM1132v1>                               | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/011/325/GCA_000011325.1_ASM1132v1>                               | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/011/325/GCA_000011325.1_ASM1132v1>                               | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/011/325/GCA_000011325.1_ASM1132v1>                               |
| txid885581 |          5.9909184| TRUE           | ASM25511v3                              | 646529  | Desulfosporosinus acidiphilus SJ4 (firmicutes)                                    | 885581       | Desulfosporosinus acidiphilus   | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/255/115/GCF_000255115.2_ASM25511v3>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/255/115/GCA_000255115.3_ASM25511v3>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/255/115/GCA_000255115.3_ASM25511v3>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/255/115/GCA_000255115.3_ASM25511v3>                              |
| txid79209  |          8.2173722| TRUE           | ASM23138v3                              | 768704  | Desulfosporosinus meridiei DSM 13257 (firmicutes)                                 | 79209        | Desulfosporosinus meridiei      | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/231/385/GCF_000231385.2_ASM23138v3>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/231/385/GCA_000231385.3_ASM23138v3>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/231/385/GCA_000231385.3_ASM23138v3>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/231/385/GCA_000231385.3_ASM23138v3>                              |
| txid102134 |          7.8365314| TRUE           | ASM23371v3                              | 767817  | Desulfotomaculum gibsoniae DSM 7213 (firmicutes)                                  | 102134       | Desulfotomaculum gibsoniae      | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/233/715/GCF_000233715.2_ASM23371v3>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/233/715/GCA_000233715.3_ASM23371v3>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/233/715/GCA_000233715.3_ASM23371v3>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/233/715/GCA_000233715.3_ASM23371v3>                              |
| txid390884 |          0.6005566| TRUE           | ASM32570v1                              | 926556  | Echinicola vietnamensis DSM 17526 (CFB group bacteria)                            | 390884       | Echinicola vietnamensis         | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/325/705/GCF_000325705.1_ASM32570v1>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/325/705/GCA_000325705.1_ASM32570v1>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/325/705/GCA_000325705.1_ASM32570v1>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/325/705/GCA_000325705.1_ASM32570v1>                              |
| txid562    |          0.2343636| TRUE           | ASM2632v1                               | 585056  | Escherichia coli UMN026 (E. coli)                                                 | 562          | Escherichia coli                | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/026/325/GCF_000026325.1_ASM2632v1>                               | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/026/325/GCA_000026325.2_ASM2632v2>                               | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/026/325/GCA_000026325.2_ASM2632v2>                               | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/026/325/GCA_000026325.2_ASM2632v2>                               |
| txid93466  |          9.8432694| TRUE           | ASM23540v3                              | 771875  | Fervidobacterium pennivorans DSM 9078 (thermotogales)                             | 93466        | Fervidobacterium pennivorans    | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/235/405/GCF_000235405.2_ASM23540v3>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/235/405/GCA_000235405.3_ASM23540v3>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/235/405/GCA_000235405.3_ASM23540v3>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/235/405/GCA_000235405.3_ASM23540v3>                              |
| txid81475  |          4.6433280| TRUE           | ASM24225v3                              | 767434  | Frateuria aurantia DSM 6220 (g-proteobacteria)                                    | 81475        | Frateuria aurantia              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/242/255/GCF_000242255.2_ASM24225v3>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/242/255/GCA_000242255.3_ASM24225v3>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/242/255/GCA_000242255.3_ASM24225v3>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/242/255/GCA_000242255.3_ASM24225v3>                              |
| txid387341 |          8.9937015| TRUE           | ASM32852v1                              | 797302  | Halovivax ruber XH-70 (euryarchaeotes)                                            | 387341       | Halovivax ruber                 | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/328/525/GCF_000328525.1_ASM32852v1>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/328/525/GCA_000328525.1_ASM32852v1>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/328/525/GCA_000328525.1_ASM32852v1>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/328/525/GCA_000328525.1_ASM32852v1>                              |
| txid2724   |          5.8590889| TRUE           | ASM2378v1                               | 582402  | Hirschia baltica ATCC 49814 (a-proteobacteria)                                    | 2724         | Hirschia baltica                | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/023/785/GCF_000023785.1_ASM2378v1>                               | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/023/785/GCA_000023785.1_ASM2378v1>                               | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/023/785/GCA_000023785.1_ASM2378v1>                               | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/023/785/GCA_000023785.1_ASM2378v1>                               |
| txid52022  |          3.1199648| TRUE           | ASM9212v1                               | 526227  | Meiothermus silvanus DSM 9946 (bacteria)                                          | 52022        | Meiothermus silvanus            | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/092/125/GCF_000092125.1_ASM9212v1>                               | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/092/125/GCA_000092125.1_ASM9212v1>                               | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/092/125/GCA_000092125.1_ASM9212v1>                               | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/092/125/GCA_000092125.1_ASM9212v1>                               |
| txid44930  |          8.3345540| TRUE           | IMG-taxon 2693429900 annotated assembly | 44930   | Natronobacterium gregoryi (euryarchaeotes)                                        | 44930        | Natronobacterium gregoryi       | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/114/025/GCF_900114025.1_IMG-taxon_2693429900_annotated_assembly> | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/114/025/GCA_900114025.1_IMG-taxon_2693429900_annotated_assembly> | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/114/025/GCA_900114025.1_IMG-taxon_2693429900_annotated_assembly> | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/114/025/GCA_900114025.1_IMG-taxon_2693429900_annotated_assembly> |
| txid29288  |         13.6663249| TRUE           | ASM32868v1                              | 694430  | Natronococcus occultus SP4 (euryarchaeotes)                                       | 29288        | Natronococcus occultus          | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/328/685/GCF_000328685.1_ASM32868v1>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/328/685/GCA_000328685.1_ASM32868v1>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/328/685/GCA_000328685.1_ASM32868v1>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/328/685/GCA_000328685.1_ASM32868v1>                              |
| txid2014   |          0.0878863| TRUE           | ASM9298v1                               | 446468  | Nocardiopsis dassonvillei subsp. dassonvillei DSM 43111 (high GC Gram+)           | 2014         | Nocardiopsis dassonvillei       | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/092/985/GCF_000092985.1_ASM9298v1>                               | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/092/985/GCA_000092985.1_ASM9298v1>                               | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/092/985/GCA_000092985.1_ASM9298v1>                               | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/092/985/GCA_000092985.1_ASM9298v1>                               |
| txid133926 |          4.4529076| TRUE           | ASM14384v1                              | 633147  | Olsenella uli DSM 7084 (actinobacteria)                                           | 133926       | Olsenella uli                   | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/143/845/GCF_000143845.1_ASM14384v1>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/143/845/GCA_000143845.1_ASM14384v1>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/143/845/GCA_000143845.1_ASM14384v1>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/143/845/GCA_000143845.1_ASM14384v1>                              |
| txid316    |          2.4022265| TRUE           | ASM21960v1                              | 316     | Pseudomonas stutzeri (g-proteobacteria)                                           | 316          | Pseudomonas stutzeri            | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/219/605/GCF_000219605.1_ASM21960v1>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/219/605/GCA_000219605.1_ASM21960v1>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/219/605/GCA_000219605.1_ASM21960v1>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/219/605/GCA_000219605.1_ASM21960v1>                              |
| txid54736  |          0.4540794| TRUE           | ASM221192v1                             | 1243617 | Salmonella bongori serovar 66:z41:- str. SA19983605 (enterobacteria)              | 54736        | Salmonella bongori              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/211/925/GCF_002211925.1_ASM221192v1>                             | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/211/925/GCA_002211925.1_ASM221192v1>                             | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/211/925/GCA_002211925.1_ASM221192v1>                             | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/211/925/GCA_002211925.1_ASM221192v1>                             |
| txid28901  |          0.5859089| TRUE           | ASM694v2                                | 99287   | Salmonella enterica subsp. enterica serovar Typhimurium str. LT2 (enterobacteria) | 28901        | Salmonella enterica             | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2>                                | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/006/945/GCA_000006945.2_ASM694v2>                                | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/006/945/GCA_000006945.2_ASM694v2>                                | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/006/945/GCA_000006945.2_ASM694v2>                                |
| txid286802 |          2.1825106| TRUE           | ASM9282v1                               | 640132  | Segniliparus rotundus DSM 44985 (high GC Gram+)                                   | 286802       | Segniliparus rotundus           | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/092/825/GCF_000092825.1_ASM9282v1>                               | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/092/825/GCA_000092825.1_ASM9282v1>                               | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/092/825/GCA_000092825.1_ASM9282v1>                               | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/092/825/GCA_000092825.1_ASM9282v1>                               |
| txid55206  |          6.8404863| TRUE           | ASM14398v1                              | 573413  | Sediminispirochaeta smaragdinae DSM 11293 (spirochetes)                           | 55206        | Sediminispirochaeta smaragdinae | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/143/985/GCF_000143985.1_ASM14398v1>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/143/985/GCA_000143985.1_ASM14398v1>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/143/985/GCA_000143985.1_ASM14398v1>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/143/985/GCA_000143985.1_ASM14398v1>                              |
| txid1314   |          0.2343636| TRUE           | ASM678v2                                | 160490  | Streptococcus pyogenes M1 GAS (firmicutes)                                        | 1314         | Streptococcus pyogenes          | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/785/GCF_000006785.2_ASM678v2>                                | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/006/785/GCA_000006785.2_ASM678v2>                                | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/006/785/GCA_000006785.2_ASM678v2>                                | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/006/785/GCA_000006785.2_ASM678v2>                                |
| txid392734 |          2.2703970| TRUE           | IMG-taxon 2690315654 annotated assembly | 392734  | Terriglobus roseus (bacteria)                                                     | 392734       | Terriglobus roseus              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/102/185/GCF_900102185.1_IMG-taxon_2690315654_annotated_assembly> | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/102/185/GCA_900102185.1_IMG-taxon_2690315654_annotated_assembly> | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/102/185/GCA_900102185.1_IMG-taxon_2690315654_annotated_assembly> | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/102/185/GCA_900102185.1_IMG-taxon_2690315654_annotated_assembly> |
| txid377615 |          0.1025341| TRUE           | ASM22770v3                              | 717605  | Thermobacillus composti KWC4 (firmicutes)                                         | 377615       | Thermobacillus composti         | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/227/705/GCF_000227705.2_ASM22770v3>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/227/705/GCA_000227705.3_ASM22770v3>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/227/705/GCA_000227705.3_ASM22770v3>                              | <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/227/705/GCA_000227705.3_ASM22770v3>                              |

The downloadStatus column indicates that we got both Genomes and CDS's when is TRUE. The usedFtpPathForGenome and usedFtpPathForCdss columns the base URL used to download these files. The taxid and organism columns give more details about the organism genome, including the strain. A problem with download.file function, is that if the downloaded file does not exist it create an empty file. We check this as follows:

``` r
composition.df$genomeCheck <- checkSequencesFiles("genomes", composition.df)
composition.df$cdsCheck <- checkSequencesFiles("cdss", composition.df)
```

If we have cases with genomes found but CDS's not found, we can use the genome to predict the CDS's using Prodigal (not the case here):

``` r
# composition.df <- predictMissingCdss(GENOMES_FOLDER, CDSS_FOLDER, composition.df, PRODIGAL_LOCATION)
```

Then we check the files again (not the case here):

``` r
# composition.df$genomeCheck <- checkSequencesFiles(GENOMES_FOLDER, composition.df)
# composition.df$cdsCheck <- checkSequencesFiles(CDSS_FOLDER, composition.df)
# composition.df$downloadStatus <- composition.df$genomeCheck & composition.df$cdsCheck
# composition.df <- composition.df %>% filter(downloadStatus)
```

Once it is not the case here, we just indicate that we did not use Prodigal for predict CDS's (a metadata necessary information):

``` r
composition.df$prodigalCds <- FALSE
```

Simulation with grinder
-----------------------

Create the references files for Grinder simulation:

``` r
# createGenomesReferenceFileForGrinder(GENOMES_FOLDER, PROJECT_FOLDER)
# createCdssReferenceFileForGrinder(CDSS_FOLDER, PROJECT_FOLDER)
```

Add amplicons information (verify primer compatibility and recompute relative abundance based on this):

``` r
# composition.df <- addAmpliconsInformation(PROJECT_FOLDER, PRIMERS_FILE, AMPLICONS_SEARCH_SCRIPT_FILE, composition.df)
```

Export reference metadata:

``` r
# exportReferenceMetada(PROJECT_FOLDER, composition.df)
```

Export ranks:

``` r
# exportRanks(PROJECT_FOLDER, composition.df)
```

Generated grinder command for simulate metagenomimcs shotgun:

``` r
cat(readLines('grinderScriptForSh.sh'), sep = '\n')
```

    ## Warning in readLines("grinderScriptForSh.sh"): incomplete final line found
    ## on 'grinderScriptForSh.sh'

    ## grinder -reference_file /work/metazoo/users/felipeprata/RProjectMBARCSimulation/genomes.fna -total_reads 1000000 -read_dist 250 normal 50 -length_bias 1 -random_seed 1 -qual_levels 30 10 -fastq_output 1 -output_dir /work/metazoo/users/felipeprata/RProjectMBARCSimulation/sh -md poly4 3e-3 3.3e-8 -mr 80 20 -af /work/metazoo/users/felipeprata/RProjectMBARCSimulation/sh-mt-ranks.txt

Generated grinder command for simulate metatranscriptomics:

``` r
cat(readLines('grinderScriptForMt.sh'), sep = '\n')
```

    ## Warning in readLines("grinderScriptForMt.sh"): incomplete final line found
    ## on 'grinderScriptForMt.sh'

    ## grinder -reference_file /work/metazoo/users/felipeprata/RProjectMBARCSimulation/cdss.fna -total_reads 1000000 -read_dist 250 normal 50 -length_bias 1 -random_seed 1 -qual_levels 30 10 -fastq_output 1 -output_dir /work/metazoo/users/felipeprata/RProjectMBARCSimulation/mt -md poly4 3e-3 3.3e-8 -mr 80 20 -af /work/metazoo/users/felipeprata/RProjectMBARCSimulation/sh-mt-ranks.txt

Generated grinder command for simulate 16S amplicons metagenomics:

``` r
cat(readLines('grinderScriptForAmp.sh'), sep = '\n')
```

    ## Warning in readLines("grinderScriptForAmp.sh"): incomplete final line found
    ## on 'grinderScriptForAmp.sh'

    ## grinder -reference_file /work/metazoo/users/felipeprata/RProjectMBARCSimulation/genomes.fna -total_reads 100000 -read_dist 250 normal 50 -copy_bias 1 -random_seed 1 -qual_levels 30 10 -fastq_output 1 -output_dir /work/metazoo/users/felipeprata/RProjectMBARCSimulation/amp -md poly4 3e-3 3.3e-8 -mr 80 20 -af /work/metazoo/users/felipeprata/RProjectMBARCSimulation/amp-ranks.txt -length_bias 0 -unidirectional 1 -chimera_perc 10 -fr /work/metazoo/users/felipeprata/RProjectMBARCSimulation/primers.fna
