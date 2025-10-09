
# NWFSC Eulachon GSI Project

Primary Contacts: Maddie Betts (madison.betts@noaa.gov), Krista Nichols (krista.nichols@noaa.gov), Mia Nahom (mia.nahom@noaa.gov)

## Background

Eulachon (*Thaleichthys pacificus*, Osmeridae) are anadromous smelts that are endemic to the North Pacific Ocean, with a range that stretches from northern California to the southeastern Bering Sea along the coast of Alaska. Eulachon are an important prey species for birds, marine mammals, and fishes, and culturally important for First Nation and American Indigenous peoples. Recent studies of the population structure of eulachon have identified a number of distinct populations. The southern population of eulachon (ranging from the Nass River in BC, Canada, to the Mad River in CA, USA) has faced sharp population declines in recent decades and was listed as threatened under the Endangered Species Act in 2010. 

## Objective

In collaboration with Fisheries and Oceans Canada (DFO), the primary objective of this project was to identify what genetic stock eulachon bycatch samples caught by the West Coast Ocean Shrimp Fishery belong to, and to what extent the fishery is impacting the threatened Southern population. To do this, we utilized high-throughput sequencing technologies and the genetic baseline designed by Sutherland et al. (2012), which included 521 variant single nucleotide polymorphisms (SNPs), to  assign each of ~1600 eulachon bycatch samples to a specific genetic stock.

## Methods

### 1. Sample Selection and Preparation for Sequencing

We selected 1,632 eulachon bycatch fin clip samples from a total of several thousand samples collected between 2013 and 2025. Samples were preserved in 95% ethanol. Samples were extracted using the DNeasy Blood and Tissue kit on a QIAcube high-throughput automated extraction machine. DNA was eluted in 200 uL EB and stored at -20 degrees C. 

### 2. GTseq 

Samples were sequenced on an Illumina MiSeq using the [Genotyping-by-Thousands](https://pubmed.ncbi.nlm.nih.gov/25476721/) (GTseq) method developed by Nate Campbell. We utilized and optimized a panel of loci originally developed by DFO for sequencing - the final panel used contained 475 loci.

### 3. GTscore

.fastq files from completed runs were demultiplexed and then analyzed using the GTscore pipeline originally developed by Garrett McKinney (see his github [here](https://github.com/gjmckinney/GTscore)). 

### 4. GSI Analysis

Samples were identified with a genetic stock ID using additional bioinformatics pipelines. 

## Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.